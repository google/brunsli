// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <brunsli/brunsli_decode.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <string>
#include <vector>

#include <brotli/decode.h>
#include "../common/constants.h"
#include "../common/context.h"
#include "../common/lehmer_code.h"
#include "../common/platform.h"
#include "../common/predict.h"
#include "../common/quant_matrix.h"
#include "./ans_decode.h"
#include "./arith_decode.h"
#include "./bit_reader.h"
#include "./brunsli_input.h"
#include "./context_map_decode.h"
#include <brunsli/jpeg_data_writer.h>
#include "./state.h"

namespace brunsli {

using ::brunsli::internal::dec::BlockI32;
using ::brunsli::internal::dec::ComponentMeta;
using ::brunsli::internal::dec::Stage;
using ::brunsli::internal::dec::State;

using ::brunsli::internal::dec::PrepareMeta;
using ::brunsli::internal::dec::UpdateSubsamplingDerivatives;

static const int kNumDirectCodes = 8;
static const int kCoeffAlphabetSize = kNumDirectCodes + 10;

static const int32_t kKnownSectionTags =
    (1u << kBrunsliSignatureTag) | (1u << kBrunsliHeaderTag) |
    (1u << kBrunsliMetaDataTag) | (1u << kBrunsliJPEGInternalsTag) |
    (1u << kBrunsliQuantDataTag) | (1u << kBrunsliHistogramDataTag) |
    (1u << kBrunsliDCDataTag) | (1u << kBrunsliACDataTag) |
    (1u << kBrunsliOriginalJpgTag);

static const int32_t kKnownHeaderVarintTags =
    (1u << kBrunsliHeaderWidthTag) | (1u << kBrunsliHeaderHeightTag) |
    (1u << kBrunsliHeaderVersionCompTag) | (1u << kBrunsliHeaderSubsamplingTag);

bool IsBrunsli(const uint8_t* data, const size_t len) {
  static const uint8_t kSignature[6] = {
      /* marker */ 0x0A,
      /* length */ 0x04, 0x42, 0xD2, 0xD5, 0x4E};
  static const size_t kSignatureLen = sizeof(kSignature);
  if (len < kSignatureLen) return false;
  return (memcmp(kSignature, data, kSignatureLen) == 0);
}

// Returns ceil(a/b).
inline int DivCeil(int a, int b) { return (a + b - 1) / b; }

// Decodes a number in the range [0..255], by reading 1 - 11 bits.
inline uint32_t DecodeVarLenUint8(BrunsliBitReader* br) {
  if (BrunsliBitReaderReadBits(br, 1)) {
    uint32_t nbits = BrunsliBitReaderReadBits(br, 3);
    if (nbits == 0) {
      return 1u;
    } else {
      return BrunsliBitReaderReadBits(br, nbits) + (1u << nbits);
    }
  }
  return 0;
}

uint32_t DecodeVarint(BrunsliBitReader* br, size_t max_bits) {
  uint32_t result = 0;
  for (size_t b = 0; b < max_bits; ++b) {
    if (b + 1 != max_bits && !BrunsliBitReaderReadBits(br, 1)) {
      break;
    }
    result |= BrunsliBitReaderReadBits(br, 1) << b;
  }
  return result;
}

size_t DecodeLimitedVarint(BrunsliBitReader* br, int nbits, int max_symbols) {
  size_t bits = 0;
  size_t shift = 0;
  for (size_t b = 0; b < max_symbols && BrunsliBitReaderReadBits(br, 1); ++b) {
    const size_t next_bits = BrunsliBitReaderReadBits(br, nbits);
    bits |= next_bits << shift;
    shift += nbits;
  }
  return bits;
}

std::string GenerateApp0Marker(uint8_t app0_status) {
  static const uint8_t kStaticApp0Data[17] = {
      0xe0, 0x00, 0x10, 'J',  'F',  'I',  'F',  0x00, 0x01,
      0x01, 0x00, 0x00, 0x01, 0x00, 0x01, 0x00, 0x00};
  std::string app0_marker(reinterpret_cast<const char*>(kStaticApp0Data),
                          sizeof(kStaticApp0Data));
  app0_marker[9] = app0_status & 1u ? 2 : 1;
  app0_status >>= 1u;
  app0_marker[10] = app0_status & 0x3u;
  app0_status >>= 2u;
  uint16_t x_dens = kApp0Densities[app0_status];
  app0_marker[11] = app0_marker[13] = x_dens >> 8u;
  app0_marker[12] = app0_marker[14] = x_dens & 0xFFu;
  return app0_marker;
}

std::string GenerateAppMarker(uint8_t marker, uint8_t code) {
  std::string s;
  if (marker == 0x80) {
    s = std::string(reinterpret_cast<const char*>(AppData_0xe2), 3161);
    s[84] = code;
  } else if (marker == 0x81) {
    s = std::string(reinterpret_cast<const char*>(AppData_0xec), 18);
    s[15] = code;
  } else {
    BRUNSLI_DCHECK(marker == 0x82);
    s = std::string(reinterpret_cast<const char*>(AppData_0xee), 15);
    s[10] = code;
  }
  return s;
}

// TODO: avoid string input.
bool AddMetaData(const std::string& metadata, JPEGData* jpg) {
  size_t pos = 0;
  size_t short_marker_count = 0;
  while (pos < metadata.size()) {
    const uint8_t marker = static_cast<uint8_t>(metadata[pos++]);
    if (marker == 0xD9) {
      jpg->tail_data = metadata.substr(pos);
      break;
    } else if (marker < 0x40) {
      if (++short_marker_count > kBrunsliShortMarkerLimit) {
        return false;
      }
      jpg->app_data.push_back(GenerateApp0Marker(marker));
    } else if (marker >= 0x80 && marker <= 0x82) {
      if (++short_marker_count > kBrunsliShortMarkerLimit) {
        return false;
      }
      if (pos >= metadata.size()) return false;
      const uint8_t code = metadata[pos++];
      jpg->app_data.push_back(GenerateAppMarker(marker, code));
    } else {
      if (pos + 1 >= metadata.size()) return false;
      const uint8_t hi = metadata[pos];
      const uint8_t lo = metadata[pos + 1];
      const size_t marker_len = (hi << 8u) + lo;
      if (marker == 0xFE) {
        jpg->com_data.push_back(metadata.substr(pos, marker_len));
      } else if ((marker >> 4u) == 0x0E) {
        jpg->app_data.push_back(metadata.substr(pos - 1, marker_len + 1));
      } else {
        return false;
      }
      pos += marker_len;
    }
  }
  return true;
}

bool DecodeQuantTables(BrunsliBitReader* br, JPEGData* jpg) {
  if (!BrunsliBitReaderReadMoreInput(br)) {
    return false;
  }
  bool have_internals_data = !jpg->quant.empty();
  size_t num_quant_tables = BrunsliBitReaderReadBits(br, 2) + 1;
  if (jpg->quant.size() != num_quant_tables) {
    return false;
  }
  for (size_t i = 0; i < num_quant_tables; ++i) {
    JPEGQuantTable* q = &jpg->quant[i];
    int data_precision = 0;
    if (!BrunsliBitReaderReadBits(br, 1)) {
      const size_t short_code = BrunsliBitReaderReadBits(br, 3);
      for (size_t k = 0; k < kDCTBlockSize; ++k) {
        q->values[k] = kStockQuantizationTables[(i > 0) ? 1 : 0][short_code][k];
      }
    } else {
      const uint32_t q_factor = BrunsliBitReaderReadBits(br, 6);
      uint8_t predictor[kDCTBlockSize];
      FillQuantMatrix(i > 0, q_factor, predictor);
      int delta = 0;
      for (size_t k = 0; k < kDCTBlockSize; ++k) {
        if (!BrunsliBitReaderReadMoreInput(br)) {
          return false;
        }
        if (BrunsliBitReaderReadBits(br, 1)) {
          const int sign = BrunsliBitReaderReadBits(br, 1);
          const int diff = DecodeVarint(br, 16) + 1;
          if (sign) {
            delta -= diff;
          } else {
            delta += diff;
          }
        }
        const int j = kJPEGNaturalOrder[k];
        const int quant_value = predictor[j] + delta;
        q->values[j] = quant_value;
        if (quant_value <= 0) {
          return false;
        }
        if (quant_value >= 256) {
          data_precision = 1;
        }
        if (quant_value >= 65536) {
          return false;
        }
      }
    }
    if (!have_internals_data) {
      q->precision = data_precision;
      q->is_last = true;
      q->index = i;
    } else {
      if (q->precision != data_precision) return false;
    }
  }
  for (size_t i = 0; i < jpg->components.size(); ++i) {
    JPEGComponent* c = &jpg->components[i];
    c->quant_idx = BrunsliBitReaderReadBits(br, 2);
    if (c->quant_idx >= jpg->quant.size()) {
      return false;
    }
  }
  return true;
}

bool DecodeHuffmanCode(BrunsliBitReader* br, JPEGHuffmanCode* huff,
                       bool is_known_last) {
  if (!BrunsliBitReaderReadMoreInput(br)) {
    return false;
  }
  huff->slot_id = BrunsliBitReaderReadBits(br, 2);
  int is_dc_table = (BrunsliBitReaderReadBits(br, 1) == 0);
  huff->slot_id += is_dc_table ? 0 : 0x10;
  huff->is_last = is_known_last || BrunsliBitReaderReadBits(br, 1);
  huff->counts[0] = 0;
  int found_match = BrunsliBitReaderReadBits(br, 1);
  if (found_match) {
    if (is_dc_table) {
      int huff_table_idx = BrunsliBitReaderReadBits(br, 1);
      memcpy(&huff->counts[1], kStockDCHuffmanCodeCounts[huff_table_idx],
             sizeof(kStockDCHuffmanCodeCounts[0]));
      memcpy(&huff->values[0], kStockDCHuffmanCodeValues[huff_table_idx],
             sizeof(kStockDCHuffmanCodeValues[0]));
    } else {
      int huff_table_idx = BrunsliBitReaderReadBits(br, 1);
      memcpy(&huff->counts[1], kStockACHuffmanCodeCounts[huff_table_idx],
             sizeof(kStockACHuffmanCodeCounts[0]));
      memcpy(&huff->values[0], kStockACHuffmanCodeValues[huff_table_idx],
             sizeof(kStockACHuffmanCodeValues[0]));
    }
    return true;
  }
  int total_count = 0;
  int space = 1u << kJpegHuffmanMaxBitLength;
  int max_len = BrunsliBitReaderReadBits(br, 4) + 1;
  int max_count = is_dc_table ? kJpegDCAlphabetSize : kJpegHuffmanAlphabetSize;
  space -= 1u << (kJpegHuffmanMaxBitLength - max_len);
  for (int i = 1; i <= max_len; ++i) {
    size_t shift = kJpegHuffmanMaxBitLength - i;
    int count_limit = std::min(max_count - total_count, space >> shift);
    if (count_limit > 0) {
      int nbits = Log2FloorNonZero(count_limit) + 1;
      int count = BrunsliBitReaderReadBits(br, nbits);
      if (count > count_limit) {
        return false;
      }
      huff->counts[i] = count;
      total_count += count;
      space -= count * (1u << shift);
    }
  }
  ++huff->counts[max_len];

  PermutationCoder p(
      is_dc_table
          ? std::vector<uint8_t>(kDefaultDCValues, std::end(kDefaultDCValues))
          : std::vector<uint8_t>(kDefaultACValues, std::end(kDefaultACValues)));
  for (int i = 0; i < total_count; ++i) {
    if (!BrunsliBitReaderReadMoreInput(br)) {
      return false;
    }

    const int nbits = p.num_bits();
    const int code = DecodeLimitedVarint(br, 2, (nbits + 1) >> 1u);
    const int value = p.Remove(code);
    if (value < 0) {
      return false;
    }
    huff->values[i] = value;
  }
  huff->values[total_count] = kJpegHuffmanAlphabetSize;
  return true;
}

bool DecodeScanInfo(BrunsliBitReader* br, JPEGScanInfo* si) {
  if (!BrunsliBitReaderReadMoreInput(br)) {
    return false;
  }
  si->Ss = BrunsliBitReaderReadBits(br, 6);
  si->Se = BrunsliBitReaderReadBits(br, 6);
  si->Ah = BrunsliBitReaderReadBits(br, 4);
  si->Al = BrunsliBitReaderReadBits(br, 4);
  si->components.resize(BrunsliBitReaderReadBits(br, 2) + 1);
  for (size_t i = 0; i < si->components.size(); ++i) {
    si->components[i].comp_idx = BrunsliBitReaderReadBits(br, 2);
    si->components[i].dc_tbl_idx = BrunsliBitReaderReadBits(br, 2);
    si->components[i].ac_tbl_idx = BrunsliBitReaderReadBits(br, 2);
  }
  int last_block_idx = -1;
  while (BrunsliBitReaderReadBits(br, 1)) {
    if (!BrunsliBitReaderReadMoreInput(br)) {
      return false;
    }
    int block_idx = last_block_idx + DecodeVarint(br, 28) + 1;
    si->reset_points.insert(block_idx);
    last_block_idx = block_idx;
    if (last_block_idx > (1u << 30u)) {
      // At most 8K x 8K x num_channels blocks are expected. That is, typically,
      // 1.5 * 2^27. 2^30 should be sufficient for any sane image.
      return false;
    }
  }
  last_block_idx = 0;
  int last_num = 0;
  while (BrunsliBitReaderReadBits(br, 1)) {
    if (!BrunsliBitReaderReadMoreInput(br)) {
      return false;
    }
    int block_idx = last_block_idx + DecodeVarint(br, 28);
    if (block_idx > last_block_idx) {
      if (last_num > 0) {
        JPEGScanInfo::ExtraZeroRunInfo info;
        info.block_idx = last_block_idx;
        info.num_extra_zero_runs = last_num;
        si->extra_zero_runs.push_back(info);
        last_num = 0;
      }
    }
    ++last_num;
    last_block_idx = block_idx;
    if (last_block_idx > (1u << 30u)) {
      // At most 8K x 8K x num_channels blocks are expected. That is, typically,
      // 1.5 * 2^27. 2^30 should be sufficient for any sane image.
      return false;
    }
  }
  if (last_num > 0) {
    JPEGScanInfo::ExtraZeroRunInfo info;
    info.block_idx = last_block_idx;
    info.num_extra_zero_runs = last_num;
    si->extra_zero_runs.push_back(info);
  }
  return true;
}

bool DecodeAuxData(BrunsliBitReader* br, JPEGData* jpg) {
  bool have_dri = false;
  int num_scans = 0;
  size_t dht_count = 0;
  uint8_t marker;
  do {
    if (!BrunsliBitReaderReadMoreInput(br)) {
      return false;
    }
    marker = 0xc0 + BrunsliBitReaderReadBits(br, 6);
    jpg->marker_order.push_back(marker);
    if (marker == 0xc4) ++dht_count;
    if (marker == 0xdd) have_dri = true;
    if (marker == 0xda) ++num_scans;
  } while (marker != 0xd9);
  if (have_dri) {
    jpg->restart_interval = BrunsliBitReaderReadBits(br, 16);
  }

  size_t terminal_huffman_code_count = 0;
  for (int i = 0;; ++i) {
    const bool is_known_last = BrunsliBitReaderReadBits(br, 1);
    JPEGHuffmanCode huff;
    if (!DecodeHuffmanCode(br, &huff, is_known_last)) {
      return false;
    }
    if (huff.is_last) {
      terminal_huffman_code_count++;
    }
    jpg->huffman_code.push_back(huff);
    if (is_known_last) break;
    if (i == kMaxDHTMarkers) {
      // Too many Huffman codes for a valid bit-stream. Normally, a jpeg file
      // can have any arbitrary number of DHT, DQT, etc. But i prefer we force a
      // reasonable lower bound instead of open door to likely forged BRN input.
      return false;
    }
  }
  if (dht_count != terminal_huffman_code_count) {
    BRUNSLI_LOG_ERROR() << "Invalid number of DHT markers" << BRUNSLI_ENDL();
    return false;
  }

  jpg->scan_info.resize(num_scans);
  for (int i = 0; i < num_scans; ++i) {
    if (!DecodeScanInfo(br, &jpg->scan_info[i])) {
      return false;
    }
  }
  int num_quant_tables = BrunsliBitReaderReadBits(br, 2) + 1;
  jpg->quant.resize(num_quant_tables);
  for (int i = 0; i < num_quant_tables; ++i) {
    JPEGQuantTable* q = &jpg->quant[i];
    q->index = BrunsliBitReaderReadBits(br, 2);
    q->is_last = (i == num_quant_tables - 1) || BrunsliBitReaderReadBits(br, 1);
    q->precision = BrunsliBitReaderReadBits(br, 4);
    if (q->precision > 1) {
      BRUNSLI_LOG_ERROR() << "Invalid quantization table precision: "
                          << q->precision << BRUNSLI_ENDL();
      return false;
    }
    // note that q->values[] are initialized to invalid 0 values.
  }
  int comp_ids = BrunsliBitReaderReadBits(br, 2);
  static const size_t kMinRequiredComponents[4] = {
      3 /* Ids123*/, 1 /* IdsGray */, 3 /* IdsRGB */, 0 /* IdsCustom */
  };
  if (jpg->components.size() < kMinRequiredComponents[comp_ids]) {
    BRUNSLI_LOG_ERROR() << "Insufficient number of components for ColorId #"
                        << comp_ids << BRUNSLI_ENDL();
    return false;
  }
  if (comp_ids == kComponentIds123) {
    jpg->components[0].id = 1;
    jpg->components[1].id = 2;
    jpg->components[2].id = 3;
  } else if (comp_ids == kComponentIdsGray) {
    jpg->components[0].id = 1;
  } else if (comp_ids == kComponentIdsRGB) {
    jpg->components[0].id = 'R';
    jpg->components[1].id = 'G';
    jpg->components[2].id = 'B';
  } else {
    BRUNSLI_DCHECK(comp_ids == kComponentIdsCustom);
    for (size_t i = 0; i < jpg->components.size(); ++i) {
      jpg->components[i].id = BrunsliBitReaderReadBits(br, 8);
    }
  }

  // security: limit is 32b for n_size
  size_t n_size = DecodeLimitedVarint(br, 8, 4);
  jpg->has_zero_padding_bit = (n_size > 0);
  if (n_size > 0) {
    if (n_size > PaddingBitsLimit(*jpg)) {
      BRUNSLI_LOG_ERROR() << "Suspicious number of padding bits " << n_size
                          << BRUNSLI_ENDL();
      return false;
    }
    jpg->padding_bits.resize(n_size);
    for (size_t i = 0; i < n_size; ++i) {
      if (!BrunsliBitReaderReadMoreInput(br)) {
        return false;
      }
      jpg->padding_bits[i] = BrunsliBitReaderReadBits(br, 1);
    }
  }
  return true;
}

bool DecodeCoeffOrder(int* order, BrunsliInput* in) {
  int lehmer[kDCTBlockSize] = {0};
  static const int kSpan = 16;
  for (int i = 0; i < kDCTBlockSize; i += kSpan) {
    if (!in->ReadBits(1)) continue;  // span is all-zero
    const int start = (i > 0) ? i : 1;
    const int end = i + kSpan;
    for (int j = start; j < end; ++j) {
      int v = 0;
      while (v <= kDCTBlockSize) {
        const int bits = in->ReadBits(3);
        v += bits;
        if (bits < 7) break;
      }
      if (v > kDCTBlockSize) return false;
      lehmer[j] = v;
    }
  }
  int end = kDCTBlockSize - 1;
  while (end >= 1 && lehmer[end] == 0) {
    --end;
  }
  if (lehmer[end] == 1) return false;
  for (int i = 1; i <= end; ++i) {
    if (lehmer[i] == 0) return false;
    --lehmer[i];
  }
  if (!DecodeLehmerCode(lehmer, kDCTBlockSize, order)) return false;
  for (int k = 0; k < kDCTBlockSize; ++k) {
    order[k] = kJPEGNaturalOrder[order[k]];
  }
  return true;
}

int DecodeNumNonzeros(Prob* const p, BinaryArithmeticDecoder* ac,
                      BrunsliInput* in) {
  const size_t kMaxBits = 6;
  int val = 1;
  for (size_t b = 0; b < kMaxBits; ++b) {
    const int bit = ac->ReadBit(p[val - 1].get_proba(), in);
    p[val - 1].Add(bit);
    val = 2 * val + bit;
  }
  return val - (1u << kMaxBits);
}

bool DecodeDC(State* state, BrunsliInput* in) {
  const std::vector<ComponentMeta>& meta = state->meta;
  const size_t num_components = meta.size();
  const int mcu_rows = meta[0].height_in_blocks / meta[0].v_samp;

  std::vector<ComponentStateDC> comps(num_components);
  for (int c = 0; c < num_components; ++c) {
    comps[c].SetWidth(meta[c].width_in_blocks);
  }

  BinaryArithmeticDecoder ac;
  ANSDecoder ans;
  ans.Init(in);
  in->InitBitReader();
  ac.Init(in);

  // We decode DC components in the following interleaved manner:
  //   v_samp[0] rows from component 0
  //   v_samp[1] rows from component 1
  //   v_samp[2] rows from component 2
  //   v_samp[3] rows from component 3 (if present)
  //
  // E.g. in a YUV420 image, we decode 2 rows of DC components from Y and then
  // 1 row of DC components from U and 1 row of DC components from V.
  for (int mcu_y = 0; mcu_y < mcu_rows; ++mcu_y) {
    for (size_t i = 0; i < num_components; ++i) {
      ComponentStateDC* c = &comps[i];
      const ComponentMeta& m = meta[i];
      const uint8_t* context_map = state->context_map + i * kNumAvrgContexts;
      const size_t ac_stride = m.ac_stride;
      const size_t b_stride = m.b_stride;
      const int width = m.width_in_blocks;
      int y = mcu_y * m.v_samp;
      int* const prev_sgn = &c->prev_sign[1];
      int* const prev_abs = &c->prev_abs_coeff[2];
      for (int iy = 0; iy < m.v_samp; ++iy, ++y) {
        coeff_t* coeffs = m.ac_coeffs + y * ac_stride;
        uint8_t* block_state = m.block_state + y * b_stride;
        for (int x = 0; x < width; ++x) {
          const int is_empty_ctx =
              IsEmptyBlockContext(&c->prev_is_nonempty[1], x);
          Prob* const is_empty_p = &c->is_empty_block_prob[is_empty_ctx];
          const bool is_empty_block = !ac.ReadBit(is_empty_p->get_proba(), in);
          is_empty_p->Add(!is_empty_block);
          c->prev_is_nonempty[x + 1] = !is_empty_block;
          *block_state = is_empty_block;
          int abs_val = 0;
          int sign = 0;
          if (!is_empty_block) {
            int is_zero = 0;
            Prob* const p_is_zero = &c->is_zero_prob;
            is_zero = ac.ReadBit(p_is_zero->get_proba(), in);
            p_is_zero->Add(is_zero);
            if (!is_zero) {
              const int avg_ctx = WeightedAverageContextDC(prev_abs, x);
              const int sign_ctx = prev_sgn[x] * 3 + prev_sgn[x - 1];
              Prob* const sign_p = &c->sign_prob[sign_ctx];
              sign = ac.ReadBit(sign_p->get_proba(), in);
              sign_p->Add(sign);
              const int entropy_ix = context_map[avg_ctx];
              int code = ans.ReadSymbol(state->entropy_codes[entropy_ix], in);
              if (code < kNumDirectCodes) {
                abs_val = code + 1;
              } else {
                const size_t nbits = code - kNumDirectCodes;
                Prob* const p_first_extra_bit = &c->first_extra_bit_prob[nbits];
                const int first_extra_bit =
                    ac.ReadBit(p_first_extra_bit->get_proba(), in);
                p_first_extra_bit->Add(first_extra_bit);
                int extra_bits_val = first_extra_bit << nbits;
                if (nbits > 0) {
                  extra_bits_val |= in->ReadBits(nbits);
                }
                abs_val = kNumDirectCodes - 1 + (2 << nbits) + extra_bits_val;
              }
            }
          }
          prev_abs[x] = abs_val;
          prev_sgn[x] = abs_val ? sign + 1 : 0;
          coeffs[0] = ((1 - 2 * sign) * abs_val +
                       PredictWithAdaptiveMedian(coeffs, x, y, ac_stride));
          ++block_state;
          coeffs += kDCTBlockSize;
        }
      }
    }
  }
  if (!ans.CheckCRC()) return false;
  return (in->error_ == 0);
}

bool DecodeAC(State* state, BrunsliInput* in) {
  const std::vector<ComponentMeta>& meta = state->meta;
  const size_t num_components = meta.size();
  const int mcu_rows = meta[0].height_in_blocks / meta[0].v_samp;

  std::vector<ComponentState> comps(num_components);
  for (size_t c = 0; c < num_components; ++c) {
    comps[c].SetWidth(meta[c].width_in_blocks);
    ComputeACPredictMultipliers(&meta[c].quant[0], comps[c].mult_row,
                                comps[c].mult_col);
  }

  BinaryArithmeticDecoder ac;
  ANSDecoder ans;
  ans.Init(in);
  in->InitBitReader();
  ac.Init(in);

  for (int i = 0; i < num_components; ++i) {
    if (!DecodeCoeffOrder(comps[i].order, in)) return false;
  }

  for (int mcu_y = 0; mcu_y < mcu_rows; ++mcu_y) {
    for (int i = 0; i < num_components; ++i) {
      ComponentState& c = comps[i];
      const ComponentMeta& m = meta[i];
      const uint8_t* context_map =
          state->context_map + m.context_offset * kNumAvrgContexts;
      const int context_bits = m.context_bits;
      const int width = m.width_in_blocks;
      const size_t ac_stride = m.ac_stride;
      const size_t b_stride = m.b_stride;
      int y = mcu_y * m.v_samp;
      int prev_row_delta = (1 - 2 * (y & 1u)) * (width + 3) * kDCTBlockSize;
      for (int iy = 0; iy < m.v_samp; ++iy, ++y) {
        coeff_t* coeffs = m.ac_coeffs + y * ac_stride;
        const coeff_t* prev_row_coeffs = coeffs - ac_stride;
        const coeff_t* prev_col_coeffs = coeffs - kDCTBlockSize;
        const uint8_t* block_state = m.block_state + y * b_stride;

        int* prev_sgn = &c.prev_sign[kDCTBlockSize];
        int* prev_abs =
            &c.prev_abs_coeff[((y & 1u) * (width + 3) + 2) * kDCTBlockSize];
        for (int x = 0; x < width; ++x) {
          const bool is_empty_block = *block_state;
          int last_nz = 0;
          if (!is_empty_block) {
            const int nonzero_ctx =
                NumNonzerosContext(&c.prev_num_nonzeros[1], x, y);
            last_nz =
                DecodeNumNonzeros(c.num_nonzero_prob[nonzero_ctx], &ac, in);
          }
          for (int k = kDCTBlockSize - 1; k > last_nz; --k) {
            prev_sgn[k] = 0;
            prev_abs[k] = 0;
          }
          int num_nonzeros = 0;
          for (int k = last_nz; k >= 1; --k) {
            int is_zero = 0;
            if (k < last_nz) {
              const int bucket = kNonzeroBuckets[num_nonzeros - 1];
              const int is_zero_ctx = bucket * kDCTBlockSize + k;
              Prob* const p = &c.is_zero_prob[is_zero_ctx];
              is_zero = ac.ReadBit(p->get_proba(), in);
              p->Add(is_zero);
            }
            int abs_val = 0;
            int sign = 1;
            const int k_nat = c.order[k];
            if (!is_zero) {
              int avg_ctx = 0;
              int sign_ctx = kMaxAverageContext;
              if (k_nat < 8) {
                if (y > 0) {
                  const int ctx = ACPredictContextRow(prev_row_coeffs + k_nat,
                                                      coeffs + k_nat,
                                                      &c.mult_col[k_nat * 8]);
                  avg_ctx = std::abs(ctx);
                  sign_ctx += ctx;
                }
              } else if ((k_nat & 7u) == 0) {
                if (x > 0) {
                  const int ctx =
                      ACPredictContextCol(prev_col_coeffs + k_nat,
                                          coeffs + k_nat, &c.mult_row[k_nat]);
                  avg_ctx = std::abs(ctx);
                  sign_ctx += ctx;
                }
              } else {
                avg_ctx = WeightedAverageContext(prev_abs + k, prev_row_delta);
                sign_ctx = prev_sgn[k] * 3 + prev_sgn[k - kDCTBlockSize];
              }
              sign_ctx = sign_ctx * kDCTBlockSize + k;
              Prob* const sign_p = &c.sign_prob[sign_ctx];
              sign = ac.ReadBit(sign_p->get_proba(), in);
              sign_p->Add(sign);
              prev_sgn[k] = sign + 1;
              sign = 1 - 2 * sign;
              const int z_dens_ctx =
                  ZeroDensityContext(num_nonzeros, k, context_bits);
              const int histo_ix = z_dens_ctx * kNumAvrgContexts + avg_ctx;
              const int entropy_ix = context_map[histo_ix];
              int code = ans.ReadSymbol(state->entropy_codes[entropy_ix], in);
              if (code < kNumDirectCodes) {
                abs_val = code + 1;
              } else {
                int nbits = code - kNumDirectCodes;
                Prob* p = &c.first_extra_bit_prob[k * 10 + nbits];
                int first_extra_bit = ac.ReadBit(p->get_proba(), in);
                p->Add(first_extra_bit);
                int extra_bits_val = first_extra_bit << nbits;
                if (nbits > 0) {
                  extra_bits_val |= in->ReadBits(nbits);
                }
                abs_val = kNumDirectCodes - 1 + (2u << nbits) + extra_bits_val;
              }
              ++num_nonzeros;
            } else {
              prev_sgn[k] = 0;
            }
            int coeff = sign * abs_val;
            coeffs[k_nat] = coeff;
            prev_abs[k] = abs_val;
          }
          c.prev_num_nonzeros[x + 1] = num_nonzeros;
          ++block_state;
          coeffs += kDCTBlockSize;
          prev_sgn += kDCTBlockSize;
          prev_abs += kDCTBlockSize;
          prev_row_coeffs += kDCTBlockSize;
          prev_col_coeffs += kDCTBlockSize;
        }
        prev_row_delta *= -1;
      }
    }
  }
  if (!ans.CheckCRC()) return false;
  return (in->error_ == 0);
}

static bool CheckCanRead(State* state, size_t required) {
  if (state->len < required) return false;
  return (state->pos <= state->len - required);
}

static bool CheckCanReadByte(State* state) { return (state->pos < state->len); }

static uint8_t ReadByte(State* state) { return state->data[state->pos++]; }

static bool DecodeBase128(State* state, size_t* val) {
  *val = 0;
  uint64_t b = 0x80;
  for (size_t i = 0; (i < 9) && (b & 0x80u); ++i) {
    if (!CheckCanReadByte(state)) return false;
    b = ReadByte(state);
    *val |= (b & 0x7Fu) << (i * 7);
  }
  return ((b & 0x80u) == 0);
}

static bool DecodeDataLength(State* state, size_t* data_len) {
  if (!DecodeBase128(state, data_len)) return false;
  return CheckCanRead(state, *data_len);
}

static Stage Fail(State* state, BrunsliStatus result) {
  state->result = result;
  return Stage::ERROR;
}

Stage VerifySignature(State* state) {
  if (!CheckCanRead(state, kBrunsliSignatureSize)) {
    return Fail(state, BRUNSLI_NOT_ENOUGH_DATA);
  }
  const bool is_signature_ok =
      (memcmp(state->data + state->pos, kBrunsliSignature,
              kBrunsliSignatureSize) != 0);
  state->pos += kBrunsliSignatureSize;
  state->tags_met |= 1u << kBrunsliSignatureTag;
  if (is_signature_ok) return Fail(state, BRUNSLI_INVALID_BRN);
  return Stage::HEADER;
}

// Parses the brunsli header starting at data[*pos] and fills in *jpg.
// Sets *pos to the position after the header.
// Returns BRUNSLI_OK, unless the data is not valid brunsli byte stream
// or is truncated.
Stage DecodeHeader(State* state, JPEGData* jpg) {
  if (!CheckCanReadByte(state)) return Fail(state, BRUNSLI_INVALID_BRN);
  const uint8_t section_tag = ReadByte(state);
  if (section_tag != SectionMarker(kBrunsliHeaderTag)) {
    return Fail(state, BRUNSLI_INVALID_BRN);
  }
  state->tags_met |= 1u << kBrunsliHeaderTag;

  size_t marker_len = 0;
  if (!DecodeDataLength(state, &marker_len)) {
    return Fail(state, BRUNSLI_INVALID_BRN);
  }
  const size_t marker_end = state->pos + marker_len;
  int32_t tags_met = 0;
  size_t varint_values[16] = {0};

  while (state->pos < marker_end) {
    const uint8_t marker = ReadByte(state);
    const size_t tag = marker >> 3u;
    if (tag == 0 || tag > 15) return Fail(state, BRUNSLI_INVALID_BRN);
    const int32_t tag_bit = 1u << tag;
    const size_t wiring_type = marker & 0x7u;
    if (wiring_type != kBrunsliWiringTypeVarint &&
        wiring_type != kBrunsliWiringTypeLengthDelimited) {
      return Fail(state, BRUNSLI_INVALID_BRN);
    }
    if (tags_met & tag_bit) return Fail(state, BRUNSLI_INVALID_BRN);
    tags_met |= tag_bit;

    const bool is_varint = (wiring_type == kBrunsliWiringTypeVarint);
    if (kKnownHeaderVarintTags & tag_bit) {
      if (!is_varint) return Fail(state, BRUNSLI_INVALID_BRN);
    }

    size_t value = 0;
    if (!DecodeBase128(state, &value)) return Fail(state, BRUNSLI_INVALID_BRN);
    if (is_varint) {
      varint_values[tag] = value;
    } else {
      if (!CheckCanRead(state, value)) return Fail(state, BRUNSLI_INVALID_BRN);
      // Skip section; current version of Brunsli does not use any.
      state->pos += value;
    }
  }
  if (state->pos != marker_end) return Fail(state, BRUNSLI_INVALID_BRN);

  const bool has_version = tags_met & (1 << kBrunsliHeaderVersionCompTag);
  if (!has_version) return Fail(state, BRUNSLI_INVALID_BRN);
  const size_t version_and_comp_count =
      varint_values[kBrunsliHeaderVersionCompTag];

  const int version = version_and_comp_count >> 2;
  jpg->version = version;

  if (version == 1) {  // fallback mode
    // TODO: do we need this?
    jpg->width = 0;
    jpg->height = 0;
    return Stage::FALLBACK;
  }

  if (version != 0) {  // unknown mode
    return Fail(state, BRUNSLI_INVALID_BRN);
  }

  // Otherwise version == 0 i.e. regular brunsli.

  // Do not allow "original_jpg" for regular Brunsli files.
  state->tags_met |= 1 << kBrunsliOriginalJpgTag;

  const bool has_width = tags_met & (1 << kBrunsliHeaderWidthTag);
  if (!has_width) return Fail(state, BRUNSLI_INVALID_BRN);
  const size_t width = varint_values[kBrunsliHeaderWidthTag];
  const bool has_height = tags_met & (1 << kBrunsliHeaderHeightTag);
  if (!has_height) return Fail(state, BRUNSLI_INVALID_BRN);
  const size_t height = varint_values[kBrunsliHeaderHeightTag];

  if (width == 0 || height == 0) return Fail(state, BRUNSLI_INVALID_BRN);
  if (width > kMaxDimPixels || height > kMaxDimPixels) {
    return Fail(state, BRUNSLI_INVALID_BRN);
  }
  jpg->width = width;
  jpg->height = height;

  const size_t num_components = (version_and_comp_count & 3) + 1;
  jpg->components.resize(num_components);

  const bool has_subsampling = tags_met & (1 << kBrunsliHeaderSubsamplingTag);
  if (!has_subsampling) return Fail(state, BRUNSLI_INVALID_BRN);
  size_t subsampling_code = varint_values[kBrunsliHeaderSubsamplingTag];

  for (size_t i = 0; i < jpg->components.size(); ++i) {
    JPEGComponent* c = &jpg->components[i];
    c->v_samp_factor = (subsampling_code & 0xF) + 1;
    subsampling_code >>= 4;
    c->h_samp_factor = (subsampling_code & 0xF) + 1;
    subsampling_code >>= 4;
    if (c->v_samp_factor > kBrunsliMaxSampling) {
      return Fail(state, BRUNSLI_INVALID_BRN);
    }
    if (c->h_samp_factor > kBrunsliMaxSampling) {
      return Fail(state, BRUNSLI_INVALID_BRN);
    }
  }
  if (!UpdateSubsamplingDerivatives(jpg)) {
    return Fail(state, BRUNSLI_INVALID_BRN);
  }

  PrepareMeta(jpg, state);

  return Stage::SECTION;
}

// Decompress Brotli stream and discard output.
// Returns true IFF stream is valid, and contains exactly |expected_output_size|
// bytes.
bool ValidateBrotliStream(const uint8_t* data, const size_t len,
                          const size_t expected_output_size) {
  BrotliDecoderState* s =
      BrotliDecoderCreateInstance(nullptr, nullptr, nullptr);
  if (s == nullptr) {
    return false;
  }

  size_t available_in = len;
  const uint8_t* next_in = data;
  size_t available_out = 0;
  BrotliDecoderResult result;
  bool sane = true;
  size_t remaining_expected_output = expected_output_size;

  while (true) {
    result = BrotliDecoderDecompressStream(s, &available_in, &next_in,
                                           &available_out, nullptr, nullptr);
    size_t output_size = 0;
    BrotliDecoderTakeOutput(s, &output_size);
    if (remaining_expected_output < output_size) {
      sane = false;
      break;
    } else {
      remaining_expected_output -= output_size;
    }
    if (result == BROTLI_DECODER_RESULT_SUCCESS) break;
    if (result == BROTLI_DECODER_RESULT_NEEDS_MORE_INPUT ||
        result == BROTLI_DECODER_RESULT_ERROR) {
      sane = false;
      break;
    }
  }
  BrotliDecoderDestroyInstance(s);

  if (available_in != 0) {
    sane = false;
  } else if (remaining_expected_output != 0) {
    sane = false;
  }

  return sane;
}

static bool DecodeMetaDataSection(State* state, JPEGData* jpg) {
  if (state->pos == state->section_end) return true;
  if (state->pos + 1 == state->section_end) {
    return AddMetaData(std::string(1, ReadByte(state)), jpg);
  }

  size_t metadata_size = 0;
  if (!DecodeBase128(state, &metadata_size)) return false;

  const uint8_t* compressed_data = state->data + state->pos;
  if (state->section_end <= state->pos) return false;
  const size_t compressed_size = state->section_end - state->pos;

  // Make additional check if compressed data is suspicious,
  // i.e. expected output is larger than 1GiB, or compression ratio is larger
  // than 4K.
  // This will protect from broken streams that would require allocating
  // gigantic chunk of memory.
  // TODO: make AddMetaData more stream-friendly; in this case temporary
  //               "metadata" string does not have to be allocated at all.
  bool is_suspicious = (metadata_size >= (1u << 30)) ||
                       ((metadata_size >> 12) > compressed_size);
  if (is_suspicious) {
    bool is_valid_brotli_stream =
        ValidateBrotliStream(compressed_data, compressed_size, metadata_size);
    if (!is_valid_brotli_stream) {
      return false;
    }
  }

  std::string metadata(metadata_size, 0);
  BrotliDecoderResult result =
      BrotliDecoderDecompress(compressed_size, compressed_data, &metadata_size,
                              reinterpret_cast<uint8_t*>(&metadata[0]));
  if (result != BROTLI_DECODER_RESULT_SUCCESS) {
    return false;
  }
  if (!AddMetaData(metadata, jpg)) {
    return false;
  }

  state->pos += compressed_size;
  return true;
}

static bool DecodeJPEGInternalsSection(State* state, JPEGData* jpg) {
  if (state->pos == state->section_end) return false;

  // TODO: merge BitReader into State
  size_t section_len = state->section_end - state->pos;
  BrunsliBitReader br;
  BrunsliBitReaderInit(&br, state->data + state->pos, section_len);
  if (!DecodeAuxData(&br, jpg)) {
    return false;
  }
  int tail_length = BrunsliBitReaderJumpToByteBoundary(&br);
  // TODO: will become no-op / DCHECK after BitReader modernization.
  if (tail_length < 0) return false;
  BRUNSLI_DCHECK(static_cast<size_t>(tail_length) < section_len);
  size_t consumed = section_len - tail_length;
  state->pos += consumed;

  for (size_t i = 0; i < jpg->marker_order.size(); ++i) {
    if (jpg->marker_order[i] != 0xff) {
      continue;
    }
    size_t data_size = 0;
    if (!DecodeDataLength(state, &data_size)) {
      return false;
    }
    jpg->inter_marker_data.emplace_back(
        reinterpret_cast<const char*>(state->data + state->pos), data_size);
    state->pos += data_size;
  }
  return true;
}

static bool DecodeQuantDataSection(State* state, JPEGData* jpg) {
  if (state->pos == state->section_end) return false;

  // TODO: merge BitReader into State
  size_t section_len = state->section_end - state->pos;
  BrunsliBitReader br;
  BrunsliBitReaderInit(&br, state->data + state->pos, section_len);
  if (!DecodeQuantTables(&br, jpg)) {
    return false;
  }

  int tail_length = BrunsliBitReaderJumpToByteBoundary(&br);
  if (tail_length != 0) return false;
  state->pos += section_len;
  return true;
}

static bool DecodeHistogramDataSection(State* state, JPEGData* jpg) {
  if (state->pos == state->section_end) return false;

  size_t num_components = jpg->components.size();
  BRUNSLI_DCHECK(num_components != 0);

  std::vector<ComponentMeta>& meta = state->meta;

  // TODO: merge BitReader into State
  size_t section_len = state->section_end - state->pos;
  BrunsliBitReader br;
  BrunsliBitReaderInit(&br, state->data + state->pos, section_len);

  size_t num_contexts = num_components;
  for (size_t i = 0; i < num_components; ++i) {
    int scheme = BrunsliBitReaderReadBits(&br, 3);
    if (scheme >= kNumSchemes) return false;
    meta[i].context_bits = scheme;
    meta[i].context_offset = num_contexts;
    num_contexts += kNumNonzeroContextSkip[scheme];
  }
  state->num_contexts = num_contexts;

  state->num_histograms = DecodeVarLenUint8(&br) + 1;

  if (!state->shallow_histograms) {
    state->context_map_.resize(state->num_contexts * kNumAvrgContexts);
    if (!DecodeContextMap(state->num_histograms, state->context_map_.size(),
                          state->context_map_.data(), &br)) {
      return false;
    }
    state->context_map = state->context_map_.data();

    state->entropy_codes_.resize(state->num_histograms);
    for (size_t i = 0; i < state->num_histograms; ++i) {
      if (!state->entropy_codes_[i].ReadFromBitStream(kCoeffAlphabetSize,
                                                      &br)) {
        return false;
      }
    }
    state->entropy_codes = state->entropy_codes_.data();

    int tail_length = BrunsliBitReaderJumpToByteBoundary(&br);
    if (tail_length != 0) return false;
  }

  state->pos += section_len;
  return true;
}

static bool DecodeDCDataSection(State* state) {
  // TODO: merge BrunsliInput into State
  size_t section_len = state->section_end - state->pos;
  BrunsliInput in(state->data + state->pos, section_len);

  if (!DecodeDC(state, &in)) return false;

  if (in.len_ != in.pos_) return false;
  state->pos += section_len;
  return true;
}

bool DecodeACDataSection(State* state) {
  // TODO: merge BrunsliInput into State
  size_t section_len = state->section_end - state->pos;
  BrunsliInput in(state->data + state->pos, section_len);

  if (!DecodeAC(state, &in)) return false;

  if (in.len_ != in.pos_) return false;
  state->pos += section_len;
  return true;
}

static Stage DecodeOriginalJpg(State* state, JPEGData* jpg) {
  if (!CheckCanReadByte(state)) return Fail(state, BRUNSLI_INVALID_BRN);
  const uint8_t tag = ReadByte(state);
  if (tag != SectionMarker(kBrunsliOriginalJpgTag)) {
    return Fail(state, BRUNSLI_INVALID_BRN);
  }
  size_t jpg_len = 0;
  if (!DecodeDataLength(state, &jpg_len)) {
    return Fail(state, BRUNSLI_INVALID_BRN);
  }
  jpg->original_jpg = state->data + state->pos;
  jpg->original_jpg_size = jpg_len;

  state->pos += jpg_len;
  return Stage::DONE;
}

static Stage ParseSection(State* state) {
  if (!CheckCanReadByte(state)) {
    // End of input.
    const bool has_ac = state->tags_met & (1 << kBrunsliACDataTag);
    if (!has_ac) {
      return Fail(state, BRUNSLI_NOT_ENOUGH_DATA);
    }
    return Stage::DONE;
  }
  const uint8_t marker = ReadByte(state);

  const size_t tag = marker >> 3;
  if (tag == 0 || tag > 15) return Fail(state, BRUNSLI_INVALID_BRN);
  const int32_t tag_bit = 1u << tag;

  const size_t wiring_type = marker & 0x7u;
  const bool is_section = (wiring_type == kBrunsliWiringTypeLengthDelimited);
  if (wiring_type != kBrunsliWiringTypeVarint &&
      wiring_type != kBrunsliWiringTypeLengthDelimited) {
    return Fail(state, BRUNSLI_INVALID_BRN);
  }

  if (state->tags_met & tag_bit) {
    BRUNSLI_LOG_ERROR() << "Duplicate marker " << std::hex
                        << static_cast<int>(marker) << BRUNSLI_ENDL();
    return Fail(state, BRUNSLI_INVALID_BRN);
  }
  state->tags_met |= tag_bit;

  const bool is_known_section_tag = kKnownSectionTags & tag_bit;
  if (!is_section) {
    if (is_known_section_tag) {
      return Fail(state, BRUNSLI_INVALID_BRN);
    } else {
      // No known varint tags on top level.
      size_t dummy;
      if (!DecodeBase128(state, &dummy)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      return Stage::SECTION;
    }
  }

  size_t marker_len = 0;
  if (!DecodeDataLength(state, &marker_len)) {
    return Fail(state, BRUNSLI_INVALID_BRN);
  }

  state->tag = tag;
  state->section_end = state->pos + marker_len;
  return Stage::SECTION_BODY;
}

static bool HasSection(State* state, int32_t tag) {
  return state->tags_met & (1 << tag);
}

static Stage ProcessSection(State* state, JPEGData* jpg) {
  const int32_t tag_bit = 1 << state->tag;
  const bool is_known_section_tag = kKnownSectionTags & tag_bit;
  if (!is_known_section_tag) {
    // Skip section content.
    state->pos = state->section_end;
    return Stage::SECTION;
  }

  if (state->skip_tags & tag_bit) {
    // Skip section content.
    state->pos = state->section_end;
    return Stage::SECTION;
  }

  switch (state->tag) {
    case kBrunsliMetaDataTag:
      if (!DecodeMetaDataSection(state, jpg)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      break;

    case kBrunsliJPEGInternalsTag:
      if (!DecodeJPEGInternalsSection(state, jpg)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      break;

    case kBrunsliQuantDataTag:
      if (!HasSection(state, kBrunsliJPEGInternalsTag)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      if (!DecodeQuantDataSection(state, jpg)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      break;

    case kBrunsliHistogramDataTag:
      if (!HasSection(state, kBrunsliJPEGInternalsTag)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      if (!DecodeHistogramDataSection(state, jpg)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      break;

    case kBrunsliDCDataTag:
      if (!HasSection(state, kBrunsliHistogramDataTag)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      if (!HasSection(state, kBrunsliQuantDataTag)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      internal::dec::WarmupMeta(jpg, state);
      if (!DecodeDCDataSection(state)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      break;

    case kBrunsliACDataTag:
      if (!HasSection(state, kBrunsliDCDataTag)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      internal::dec::WarmupMeta(jpg, state);
      if (!DecodeACDataSection(state)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      break;

    default:
      /* Unreachable */
      return Fail(state, BRUNSLI_INVALID_BRN);
  }

  if (state->pos != state->section_end) {
    return Fail(state, BRUNSLI_INVALID_BRN);
  }

  return Stage::SECTION;
}

namespace internal {
namespace dec {

bool UpdateSubsamplingDerivatives(JPEGData* jpg) {
  for (size_t i = 0; i < jpg->components.size(); ++i) {
    JPEGComponent* c = &jpg->components[i];
    jpg->max_h_samp_factor = std::max(jpg->max_h_samp_factor, c->h_samp_factor);
    jpg->max_v_samp_factor = std::max(jpg->max_v_samp_factor, c->v_samp_factor);
  }
  jpg->MCU_rows = DivCeil(jpg->height, jpg->max_v_samp_factor * 8);
  jpg->MCU_cols = DivCeil(jpg->width, jpg->max_h_samp_factor * 8);
  for (size_t i = 0; i < jpg->components.size(); ++i) {
    JPEGComponent* c = &jpg->components[i];
    c->width_in_blocks = jpg->MCU_cols * c->h_samp_factor;
    c->height_in_blocks = jpg->MCU_rows * c->v_samp_factor;
    int64_t num_blocks =
        static_cast<int64_t>(c->width_in_blocks) * c->height_in_blocks;
    if (num_blocks > kBrunsliMaxNumBlocks) {
      return false;
    }
    c->num_blocks = num_blocks;
  }
  return true;
}

void PrepareMeta(const JPEGData* jpg, State* state) {
  size_t num_components = jpg->components.size();
  state->block_state_.resize(num_components);
  std::vector<ComponentMeta>& meta = state->meta;
  meta.resize(num_components);
  for (size_t i = 0; i < num_components; ++i) {
    const JPEGComponent& c = jpg->components[i];
    ComponentMeta& m = meta[i];
    m.h_samp = c.h_samp_factor;
    m.v_samp = c.v_samp_factor;
    m.width_in_blocks = jpg->MCU_cols * m.h_samp;
    m.height_in_blocks = jpg->MCU_rows * m.v_samp;
  }
}

void WarmupMeta(JPEGData* jpg, State* state) {
  std::vector<ComponentMeta>& meta = state->meta;
  const size_t num_components = meta.size();

  if (!state->is_storage_allocated) {
    state->is_storage_allocated = true;
    for (size_t i = 0; i < num_components; ++i) {
      size_t num_blocks = meta[i].width_in_blocks * meta[i].height_in_blocks;
      jpg->components[i].coeffs.resize(num_blocks * kDCTBlockSize);
      state->block_state_[i].resize(num_blocks);
      meta[i].block_state = state->block_state_[i].data();
    }
  }

  if (!state->is_meta_warm) {
    state->is_meta_warm = true;
    for (size_t c = 0; c < num_components; ++c) {
      ComponentMeta& m = meta[c];
      const JPEGQuantTable& q = jpg->quant[jpg->components[c].quant_idx];
      m.ac_coeffs = jpg->components[c].coeffs.data();
      m.ac_stride = m.width_in_blocks * kDCTBlockSize;
      m.b_stride = m.width_in_blocks;
      memcpy(m.quant.data(), q.values.data(),
             kDCTBlockSize * sizeof(m.quant[0]));
    }
  }
}

BrunsliStatus ProcessJpeg(State* state, JPEGData* jpg) {
  while (true) {
    switch (state->stage) {
      case Stage::SIGNATURE:
        state->stage = VerifySignature(state);
        break;

      case Stage::HEADER:
        state->stage = DecodeHeader(state, jpg);
        break;

      case Stage::FALLBACK:
        state->stage = DecodeOriginalJpg(state, jpg);
        break;

      case Stage::SECTION:
        state->stage = ParseSection(state);
        break;

      case Stage::SECTION_BODY:
        state->stage = ProcessSection(state, jpg);
        break;

      case Stage::DONE:
        // It is expected that there is no garbage after the valid brunsli
        // stream.
        if (state->pos != state->len) return BRUNSLI_INVALID_BRN;
        return BRUNSLI_OK;

      case Stage::ERROR:
        return state->result;

      default:
        /* Unreachable */
        return BRUNSLI_DECOMPRESSION_ERROR;
    }
  }
}

}  // namespace dec
}  // namespace internal

BrunsliStatus BrunsliDecodeJpeg(const uint8_t* data, const size_t len,
                                JPEGData* jpg) {
  if (!data) return BRUNSLI_INVALID_PARAM;

  State state;
  state.data = data;
  state.len = len;

  return internal::dec::ProcessJpeg(&state, jpg);
}

size_t BrunsliEstimateDecoderPeakMemoryUsage(const uint8_t* data,
                                             const size_t len) {
  if (!data) return BRUNSLI_INVALID_PARAM;

  State state;
  state.data = data;
  state.len = len;
  state.skip_tags = ~(1u << kBrunsliHistogramDataTag);
  state.shallow_histograms = true;

  JPEGData jpg;
  BrunsliStatus status = internal::dec::ProcessJpeg(&state, &jpg);

  if (status != BRUNSLI_OK) return 0;

  size_t out_size = 2 * len;
  size_t total_num_blocks = 0;
  size_t component_state_size = 0;
  for (size_t i = 0; i < jpg.components.size(); ++i) {
    const JPEGComponent& c = jpg.components[i];
    total_num_blocks += c.num_blocks;
    component_state_size += ComponentState::SizeInBytes(c.width_in_blocks);
  }
  size_t jpeg_data_size = total_num_blocks * kDCTBlockSize * sizeof(coeff_t);
  size_t context_map_size =
      state.num_contexts * kNumAvrgContexts * sizeof(int32_t);
  size_t histogram_size = state.num_histograms * sizeof(ANSDecodingData);
  size_t decode_peak = context_map_size + histogram_size + component_state_size;
  size_t jpeg_writer_size = (1 << 17) + (1 << 16) * sizeof(int32_t);
  return (out_size + jpeg_data_size + std::max(decode_peak, jpeg_writer_size));
}

}  // namespace brunsli
