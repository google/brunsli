// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./brunsli_decode.h"

#include <algorithm>
#include <cstdlib>
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
#include "./brunsli_aux_data.h"
#include "./brunsli_input.h"
#include "./context_map_decode.h"
#include "./jpeg_data_writer.h"

namespace brunsli {

static const int kNumDirectCodes = 8;
static const int kCoeffAlphabetSize = kNumDirectCodes + 10;

// Returns ceil(a/b).
inline int DivCeil(int a, int b) {
  return (a + b - 1) / b;
}

// Decodes a number in the range [0..255], by reading 1 - 11 bits.
inline int DecodeVarLenUint8(BrunsliBitReader* br) {
  if (BrunsliBitReaderReadBits(br, 1)) {
    int nbits = static_cast<int>(BrunsliBitReaderReadBits(br, 3));
    if (nbits == 0) {
      return 1;
    } else {
      return static_cast<int>(BrunsliBitReaderReadBits(br, nbits)) +
             (1 << nbits);
    }
  }
  return 0;
}

int DecodeVarint(BrunsliBitReader* br, int max_bits) {
  int n = 0;
  for (int b = 0; b < max_bits; ++b) {
    if (b + 1 != max_bits && !BrunsliBitReaderReadBits(br, 1)) {
      break;
    }
    n |= BrunsliBitReaderReadBits(br, 1) << b;
  }
  return n;
}

size_t DecodeLimitedVarint(BrunsliBitReader* br, int nbits, int max_symbols) {
  size_t bits = 0;
  int shift = 0;
  for (int b = 0; b < max_symbols && BrunsliBitReaderReadBits(br, 1); ++b) {
    const size_t next_bits =
    static_cast<size_t>(BrunsliBitReaderReadBits(br, nbits));
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
  app0_marker[9] = app0_status & 1 ? 2 : 1;
  app0_status >>= 1;
  app0_marker[10] = app0_status & 0x3;
  app0_status >>= 2;
  int x_dens = kApp0Densities[app0_status];
  app0_marker[11] = app0_marker[13] = x_dens >> 8;
  app0_marker[12] = app0_marker[14] = x_dens & 0xff;
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

bool AddMetaData(const std::string& metadata, JPEGData* jpg) {
  size_t pos = 0;
  size_t short_marker_count = 0;
  while (pos < metadata.size()) {
    const uint8_t marker = metadata[pos++];
    if (marker == 0xd9) {
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
      const size_t marker_len = (hi << 8) + lo;
      if (marker == 0xfe) {
        jpg->com_data.push_back(metadata.substr(pos, marker_len));
      } else if ((marker >> 4) == 0xe) {
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
      const int qfactor = BrunsliBitReaderReadBits(br, 6);
      uint8_t predictor[kDCTBlockSize];
      FillQuantMatrix(i > 0, qfactor, predictor);
      int delta = 0;
      for (int k = 0; k < kDCTBlockSize; ++k) {
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

bool DecodeHuffmanCode(BrunsliBitReader* br,
                       JPEGHuffmanCode* huff,
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
  int space = 1 << kJpegHuffmanMaxBitLength;
  int max_len = BrunsliBitReaderReadBits(br, 4) + 1;
  int max_count = is_dc_table ? kJpegDCAlphabetSize : kJpegHuffmanAlphabetSize;
  space -= (1 << (kJpegHuffmanMaxBitLength - max_len));
  for (int i = 1; i <= max_len; ++i) {
    int count_limit = std::min(max_count - total_count,
                               space >> (kJpegHuffmanMaxBitLength - i));
    if (count_limit > 0) {
      int nbits = Log2FloorNonZero(count_limit) + 1;
      int count = BrunsliBitReaderReadBits(br, nbits);
      if (count > count_limit) {
        return false;
      }
      huff->counts[i] = count;
      total_count += count;
      space -= count * (1 << (kJpegHuffmanMaxBitLength - i));
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
    const int code = DecodeLimitedVarint(br, 2, (nbits + 1) >> 1);
    const int value = p.Remove(code);
    if (value < 0) {
      return false;
    }
    huff->values[i] = value;
  }
  huff->values[total_count] = kJpegHuffmanAlphabetSize;
  return true;
}

bool DecodeScanInfo(BrunsliBitReader* br,
                    JPEGScanInfo* si) {
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
    if (last_block_idx > (1 << 30)) {
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
    if (last_block_idx > (1 << 30)) {
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
  for (int i = 0; ; ++i) {
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
      // Too many Huffman codes for a valid bitstream. Normally, a jpeg file can
      // have any arbitrary number of DHT, DQT, etc. But i prefer we force a
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
  for (int i = 0 ; i < num_quant_tables; ++i) {
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

  // security: limit is 32b for nsize
  size_t nsize = DecodeLimitedVarint(br, 8, 4);
  jpg->has_zero_padding_bit = (nsize > 0);
  if (nsize > 0) {
    if (nsize > PaddingBitsLimit(*jpg)) {
      BRUNSLI_LOG_ERROR() << "Suspicious number of padding bits " << nsize
                          << BRUNSLI_ENDL();
      return false;
    }
    jpg->padding_bits.resize(nsize);
    for (size_t i = 0; i < nsize; ++i) {
      if (!BrunsliBitReaderReadMoreInput(br)) {
        return false;
      }
      jpg->padding_bits[i] = BrunsliBitReaderReadBits(br, 1);
    }
  }
  return true;
}

bool DecodeCoeffOrder(int* order, BrunsliInput* in) {
  int lehmer[kDCTBlockSize] = { 0 };
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
  const int kMaxBits = 6;
  int val = 1;
  for (int b = 0; b < kMaxBits; ++b) {
    const int bit = ac->ReadBit(p[val - 1].get_proba(), in);
    p[val - 1].Add(bit);
    val = 2 * val + bit;
  }
  return val - (1 << kMaxBits);
}

bool DecodeDC(int mcu_cols,
              int mcu_rows,
              int num_components,
              const int h_samp[kMaxComponents],
              const int v_samp[kMaxComponents],
              const std::vector<uint8_t>& context_map,
              const std::vector<ANSDecodingData>& entropy_codes,
              coeff_t* all_coeffs[kMaxComponents],
              std::vector<bool>* const block_state,
              BrunsliInput* in) {
  std::vector<ComponentStateDC> comps(num_components);
  int total_num_blocks = 0;
  for (int i = 0; i < num_components; ++i) {
    comps[i].SetWidth(mcu_cols * h_samp[i]);
    total_num_blocks += mcu_cols * mcu_rows * h_samp[i] * v_samp[i];
  }
  block_state->resize(total_num_blocks);

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
  int block_ipos = 0;
  for (int mcu_y = 0; mcu_y < mcu_rows; ++mcu_y) {
    for (int i = 0; i < num_components; ++i) {
      ComponentStateDC* const c = &comps[i];
      const uint8_t* const cur_context_map = &context_map[i * kNumAvrgContexts];
      const int width = c->width;
      int y = mcu_y * v_samp[i];
      int block_ix = y * width;
      coeff_t* coeffs = &all_coeffs[i][block_ix * kDCTBlockSize];
      int* const prev_sgn = &c->prev_sign[1];
      int* const prev_abs = &c->prev_abs_coeff[2];
      for (int iy = 0; iy < v_samp[i]; ++iy, ++y) {
        for (int x = 0; x < width; ++x) {
          const int is_empty_ctx =
              IsEmptyBlockContext(&c->prev_is_nonempty[1], x);
          Prob* const is_empty_p = &c->is_empty_block_prob[is_empty_ctx];
          const bool is_empty_block = !ac.ReadBit(is_empty_p->get_proba(), in);
          is_empty_p->Add(!is_empty_block);
          c->prev_is_nonempty[x + 1] = !is_empty_block;
          (*block_state)[block_ipos] = is_empty_block;
          int absval = 0;
          int sign = 0;
          if (!is_empty_block) {
            int is_zero = 0;
            Prob* const p = &c->is_zero_prob;
            is_zero = ac.ReadBit(p->get_proba(), in);
            p->Add(is_zero);
            if (!is_zero) {
              const int avrg_ctx = WeightedAverageContextDC(prev_abs, x);
              const int sign_ctx = prev_sgn[x] * 3 + prev_sgn[x - 1];
              Prob* const sign_p = &c->sign_prob[sign_ctx];
              sign = ac.ReadBit(sign_p->get_proba(), in);
              sign_p->Add(sign);
              const int entropy_ix = cur_context_map[avrg_ctx];
              int code = ans.ReadSymbol(entropy_codes[entropy_ix], in);
              if (code < kNumDirectCodes) {
                absval = code + 1;
              } else {
                const int nbits = code - kNumDirectCodes;
                Prob* const p = &c->first_extra_bit_prob[nbits];
                const int first_extra_bit = ac.ReadBit(p->get_proba(), in);
                p->Add(first_extra_bit);
                int extra_bits_val = first_extra_bit << nbits;
                if (nbits > 0) {
                  extra_bits_val |= in->ReadBits(nbits);
                }
                absval = kNumDirectCodes - 1 + (2 << nbits) + extra_bits_val;
              }
            }
          }
          prev_abs[x] = absval;
          prev_sgn[x] = absval ? sign + 1 : 0;
          coeffs[0] = ((1 - 2 * sign) * absval +
                       PredictWithAdaptiveMedian(coeffs, x, y, width));
          ++block_ipos;
          ++block_ix;
          coeffs += kDCTBlockSize;
        }
      }
    }
  }
  if (!ans.CheckCRC()) return false;
  if (in->error_) return false;
  return true;
}

bool DecodeAC(const int mcu_cols,
              const int mcu_rows,
              const int num_components,
              const int h_samp[kMaxComponents],
              const int v_samp[kMaxComponents],
              const int all_quant[kMaxComponents][kDCTBlockSize],
              const std::vector<int> context_bits,
              const std::vector<uint8_t>& context_map,
              const std::vector<ANSDecodingData>& entropy_codes,
              const std::vector<bool> block_state,
              coeff_t* all_coeffs[kMaxComponents],
              BrunsliInput* in) {
  int num_contexts = num_components;
  std::vector<ComponentState> comps(num_components);
  for (int i = 0; i < num_components; ++i) {
    comps[i].SetWidth(mcu_cols * h_samp[i]);
    comps[i].context_offset = num_contexts * kNumAvrgContexts;
    num_contexts += kNumNonzeroContextSkip[context_bits[i]];
    ComputeACPredictMultipliers(&all_quant[i][0],
                                &comps[i].mult_row[0],
                                &comps[i].mult_col[0]);
  }

  BinaryArithmeticDecoder ac;
  ANSDecoder ans;
  ans.Init(in);
  in->InitBitReader();
  ac.Init(in);

  for (int i = 0; i < num_components; ++i) {
    if (!DecodeCoeffOrder(&comps[i].order[0], in)) {
      return false;
    }
  }

  int block_ipos = 0;
  for (int mcu_y = 0; mcu_y < mcu_rows; ++mcu_y) {
    for (int i = 0; i < num_components; ++i) {
      ComponentState* const c = &comps[i];
      const uint8_t* const cur_context_map = &context_map[c->context_offset];
      const int cur_ctx_bits = context_bits[i];
      const int width = c->width;
      int y = mcu_y * v_samp[i];
      int block_ix = y * width;
      coeff_t* coeffs = &all_coeffs[i][block_ix * kDCTBlockSize];
      const coeff_t* prev_row_coeffs =
          &all_coeffs[i][(block_ix - width) * kDCTBlockSize];
      const coeff_t* prev_col_coeffs =
          &all_coeffs[i][(block_ix - 1) * kDCTBlockSize];
      int prev_row_delta =
          (1 - 2 * (y & 1)) * (width + 3) * kDCTBlockSize;
      for (int iy = 0; iy < v_samp[i]; ++iy, ++y) {
        int* prev_sgn = &c->prev_sign[kDCTBlockSize];
        int* prev_abs =
            &c->prev_abs_coeff[((y & 1) * (width + 3) + 2) * kDCTBlockSize];
        for (int x = 0; x < width; ++x) {
          const bool is_empty_block = block_state[block_ipos];
          int last_nz = 0;
          if (!is_empty_block) {
            const int nzero_ctx =
                NumNonzerosContext(&c->prev_num_nonzeros[1], x, y);
            last_nz =
                DecodeNumNonzeros(c->num_nonzero_prob[nzero_ctx], &ac, in);
          }
          for (int k = kDCTBlockSize - 1; k > last_nz; --k) {
            prev_sgn[k] = 0;
            prev_abs[k] = 0;
          }
          int num_nzeros = 0;
          for (int k = last_nz; k >= 1; --k) {
            int is_zero = 0;
            if (k < last_nz) {
              const int bucket = kNonzeroBuckets[num_nzeros - 1];
              const int is_zero_ctx = bucket * kDCTBlockSize + k;
              Prob* const p = &c->is_zero_prob[is_zero_ctx];
              is_zero = ac.ReadBit(p->get_proba(), in);
              p->Add(is_zero);
            }
            int absval = 0;
            int sign = 1;
            const int k_nat = c->order[k];
            if (!is_zero) {
              int avrg_ctx = 0;
              int sign_ctx = kMaxAverageContext;
              if (k_nat < 8) {
                if (y > 0) {
                  const int ctx = ACPredictContextRow(prev_row_coeffs + k_nat,
                                                      coeffs + k_nat,
                                                      &c->mult_col[k_nat * 8]);
                  avrg_ctx = std::abs(ctx);
                  sign_ctx += ctx;
                }
              } else if ((k_nat & 7) == 0) {
                if (x > 0) {
                  const int ctx = ACPredictContextCol(prev_col_coeffs + k_nat,
                                                      coeffs + k_nat,
                                                      &c->mult_row[k_nat]);
                  avrg_ctx = std::abs(ctx);
                  sign_ctx += ctx;
                }
              } else {
                avrg_ctx = WeightedAverageContext(prev_abs + k,
                                                  prev_row_delta);
                sign_ctx = prev_sgn[k] * 3 + prev_sgn[k - kDCTBlockSize];
              }
              sign_ctx = sign_ctx * kDCTBlockSize + k;
              Prob* const sign_p = &c->sign_prob[sign_ctx];
              sign = ac.ReadBit(sign_p->get_proba(), in);
              sign_p->Add(sign);
              prev_sgn[k] = sign + 1;
              sign = 1 - 2 * sign;
              const int zdens_ctx =
                  ZeroDensityContext(num_nzeros, k, cur_ctx_bits);
              const int histo_ix = zdens_ctx * kNumAvrgContexts + avrg_ctx;
              const int entropy_ix = cur_context_map[histo_ix];
              int code = ans.ReadSymbol(entropy_codes[entropy_ix], in);
              if (code < kNumDirectCodes) {
                absval = code + 1;
              } else {
                int nbits = code - kNumDirectCodes;
                Prob* p = &c->first_extra_bit_prob[k * 10 + nbits];
                int first_extra_bit = ac.ReadBit(p->get_proba(), in);
                p->Add(first_extra_bit);
                int extra_bits_val = first_extra_bit << nbits;
                if (nbits > 0) {
                  extra_bits_val |= in->ReadBits(nbits);
                }
                absval = kNumDirectCodes - 1 + (2 << nbits) + extra_bits_val;
              }
              ++num_nzeros;
            } else {
              prev_sgn[k] = 0;
            }
            int coeff = sign * absval;
            coeffs[k_nat] = coeff;
            prev_abs[k] = absval;
          }
          c->prev_num_nonzeros[x + 1] = num_nzeros;
          ++block_ipos;
          ++block_ix;
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
  if (in->error_) return false;
  return true;
}

struct JPEGDecodingState {
  std::vector<int> context_bits;
  std::vector<uint8_t> context_map;
  std::vector<ANSDecodingData> entropy_codes;
  std::vector<bool> block_state;
};

bool DecodeBase128(const uint8_t* data, const size_t len,
                   size_t* pos, size_t* val) {
  int shift = 0;
  uint64_t b;
  *val = 0;
  do {
    if (*pos >= len || shift > 57) {
      return false;
    }
    b = data[(*pos)++];
    *val |= (b & 0x7f) << shift;
    shift += 7;
  } while (b & 0x80);
  return true;
}

bool DecodeDataLength(const uint8_t* data, const size_t len,
                      size_t* pos, size_t* data_len) {
  if (!DecodeBase128(data, len, pos, data_len)) {
    return false;
  }
  return *data_len <= len && *pos <= len - *data_len;
}

// Returns true if there is a brunsli signature starting at data[*pos].
// Sets *pos to the position after the signature.
BrunsliStatus VerifySignature(const uint8_t* data, const size_t len,
                              size_t* pos) {
  if (data == NULL || len < kBrunsliSignatureSize) {
    return BRUNSLI_NOT_ENOUGH_DATA;
  }
  if (*pos > len - kBrunsliSignatureSize ||
      memcmp(&data[*pos], kBrunsliSignature, kBrunsliSignatureSize) != 0) {
    return BRUNSLI_INVALID_BRN;
  }
  *pos += kBrunsliSignatureSize;
  return BRUNSLI_OK;
}

// Parses the brunsli header starting at data[*pos] and fills in *jpg.
// Sets *pos to the position after the header.
// Returns BRUNSLI_OK, unless the data is not valid brunsli byte stream
// or is truncated.
BrunsliStatus DecodeHeader(const uint8_t* data, const size_t len, size_t* pos,
                           JPEGData* jpg) {
  if (*pos >= len || data[(*pos)++] != SectionMarker(kBrunsliHeaderTag)) {
    return BRUNSLI_INVALID_BRN;
  }
  size_t marker_len = 0;
  if (!DecodeDataLength(data, len, pos, &marker_len)) {
    return BRUNSLI_INVALID_BRN;
  }
  size_t marker_end = *pos + marker_len;
  bool tags_met[16] = { false };
  bool known_varint_tags[16] = {false};
  for (size_t tag :
       {kBrunsliHeaderWidthTag, kBrunsliHeaderHeightTag,
        kBrunsliHeaderVersionCompTag, kBrunsliHeaderSubsamplingTag}) {
    known_varint_tags[tag] = true;
  }
  size_t varint_values[16] = {0};

  while (*pos < marker_end) {
    const uint8_t marker = data[(*pos)++];
    const size_t tag = marker >> 3;
    if (tag == 0 || tag > 15) return BRUNSLI_INVALID_BRN;
    const size_t wiring_type = marker & 0x7;
    if (wiring_type != kBrunsliWiringTypeVarint &&
        wiring_type != kBrunsliWiringTypeLengthDelimited) {
      return BRUNSLI_INVALID_BRN;
    }

    if (tags_met[tag]) {
      return BRUNSLI_INVALID_BRN;
    }
    tags_met[tag] = true;

    const bool is_varint = (wiring_type == kBrunsliWiringTypeVarint);
    if (known_varint_tags[tag] && !is_varint) return BRUNSLI_INVALID_BRN;

    size_t value = 0;
    if (!DecodeBase128(data, len, pos, &value)) return BRUNSLI_INVALID_BRN;
    if (is_varint) {
      varint_values[tag] = value;
    } else {
      // Skip section; current version of Brunsli does not use any.
      if (value > len || *pos > len - value) return BRUNSLI_INVALID_BRN;
      *pos += value;
    }
  }
  if (*pos != marker_end) {
    return BRUNSLI_INVALID_BRN;
  }

  if (!tags_met[kBrunsliHeaderVersionCompTag]) return BRUNSLI_INVALID_BRN;
  const size_t version_and_comp_count =
      varint_values[kBrunsliHeaderVersionCompTag];

  const int version = version_and_comp_count >> 2;
  jpg->version = version;

  if (version == 0) {  // regular brunsli
    if (!tags_met[kBrunsliHeaderWidthTag]) return BRUNSLI_INVALID_BRN;
    const size_t width = varint_values[kBrunsliHeaderWidthTag];
    if (!tags_met[kBrunsliHeaderHeightTag]) return BRUNSLI_INVALID_BRN;
    const size_t height = varint_values[kBrunsliHeaderHeightTag];

    if (width == 0 || height == 0) return BRUNSLI_INVALID_BRN;
    if (width > kMaxDimPixels || height > kMaxDimPixels) {
      return BRUNSLI_INVALID_BRN;
    }
    jpg->width = width;
    jpg->height = height;

    const size_t comp_count = (version_and_comp_count & 3) + 1;
    jpg->components.resize(comp_count);

    if (!tags_met[kBrunsliHeaderSubsamplingTag]) return BRUNSLI_INVALID_BRN;
    size_t subsampling_code = varint_values[kBrunsliHeaderSubsamplingTag];

    for (size_t i = 0; i < jpg->components.size(); ++i) {
      JPEGComponent* c = &jpg->components[i];
      c->v_samp_factor = (subsampling_code & 0xF) + 1;
      if (c->v_samp_factor > kBrunsliMaxSampling) return BRUNSLI_INVALID_BRN;
      subsampling_code >>= 4;
      c->h_samp_factor = (subsampling_code & 0xF) + 1;
      if (c->h_samp_factor > kBrunsliMaxSampling) return BRUNSLI_INVALID_BRN;
      subsampling_code >>= 4;
    }
    if (!UpdateSubsamplingDerivatives(jpg)) return BRUNSLI_INVALID_BRN;
  } else if (version == 1) {  // fallback mode
    // TODO: do we need this?
    jpg->width = 0;
    jpg->height = 0;
  } else {
    return BRUNSLI_INVALID_BRN;
  }

  return BRUNSLI_OK;
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
    result = BrotliDecoderDecompressStream(s,
        &available_in, &next_in, &available_out, nullptr, nullptr);
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

bool DecodeMetaDataSection(const uint8_t* data, const size_t len,
                           JPEGData* jpg) {
  if (len == 0) {
    return true;
  } else if (len == 1) {
    return AddMetaData(std::string(1, data[0]), jpg);
  }

  size_t pos = 0;
  size_t metadata_size = 0;
  if (!DecodeBase128(data, len, &pos, &metadata_size)) {
    return false;
  }
  size_t compressed_metadata_size = len - pos;

  // Make additional check if compressed data is suspicious,
  // i.e. expected output is larger than 1GiB, or compression ratio is larger
  // than 4K.
  // This will protect from broken streams that would require allocating
  // gigantic chunk of memory.
  // TODO: make AddMetaData more stream-friendly; in this case temporary
  //               "metadata" string does not have to be allocated at all.
  bool is_suspicious = (metadata_size >= (1u << 30)) ||
                       ((metadata_size >> 12) > compressed_metadata_size);
  if (is_suspicious) {
    bool is_valid_brotli_stream = ValidateBrotliStream(
        &data[pos], compressed_metadata_size, metadata_size);
    if (!is_valid_brotli_stream) {
      return false;
    }
  }

  std::string metadata(metadata_size, 0);
  BrotliDecoderResult result = BrotliDecoderDecompress(
      compressed_metadata_size, &data[pos], &metadata_size,
      reinterpret_cast<uint8_t*>(&metadata[0]));
  if (result != BROTLI_DECODER_RESULT_SUCCESS) {
    return false;
  }
  if (!AddMetaData(metadata, jpg)) {
    return false;
  }
  return true;
}

bool DecodeJPEGInternalsSection(const uint8_t* data, const size_t len,
                                JPEGData* jpg) {
  if (len == 0) {
    return false;
  }
  BrunsliBitReader br;
  BrunsliBitReaderInit(&br, data, len);
  if (!DecodeAuxData(&br, jpg)) {
    return false;
  }
  int tail_length = BrunsliBitReaderJumpToByteBoundary(&br);
  // TODO: will become no-op / DCHECK after BitReader modernization.
  if (tail_length < 0) return false;
  BRUNSLI_DCHECK(static_cast<size_t>(tail_length) < len);
  size_t pos = len - tail_length;
  for (size_t i = 0; i < jpg->marker_order.size(); ++i) {
    if (jpg->marker_order[i] != 0xff) {
      continue;
    }
    size_t data_size = 0;
    if (!DecodeDataLength(data, len, &pos, &data_size)) {
      return false;
    }
    jpg->inter_marker_data.push_back(
        std::string(reinterpret_cast<const char*>(&data[pos]), data_size));
    pos += data_size;
  }
  return (pos == len);
}

bool DecodeQuantDataSection(const uint8_t* data, const size_t len,
                            JPEGData* jpg) {
  if (len == 0) {
    return false;
  }
  BrunsliBitReader br;
  BrunsliBitReaderInit(&br, data, len);
  if (!DecodeQuantTables(&br, jpg)) {
    return false;
  }

  return (BrunsliBitReaderJumpToByteBoundary(&br) == 0);
}

bool DecodeHistogramDataSection(const uint8_t* data, const size_t len,
                                BrunsliReadMode mode,
                                JPEGDecodingState* s,
                                JPEGData* jpg,
                                BrunsliAuxData* aux) {
  if (jpg->components.empty()) {
    // Histogram data can not be decoded without knowing the number of
    // components from the header.
    return false;
  }
  if (len == 0) {
    return false;
  }
  BrunsliBitReader br;
  BrunsliBitReaderInit(&br, data, len);
  int num_contexts = jpg->components.size();
  s->context_bits.resize(jpg->components.size());
  for (size_t i = 0; i < jpg->components.size(); ++i) {
    int scheme = BrunsliBitReaderReadBits(&br, 3);
    if (scheme >= kNumSchemes) {
      return false;
    }
    s->context_bits[i] = scheme;
    num_contexts += kNumNonzeroContextSkip[scheme];
  }
  int num_histograms = DecodeVarLenUint8(&br) + 1;
  if (mode == BRUNSLI_READ_SIZES) {
    if (aux != nullptr) {
      aux->num_brunsli_contexts = num_contexts;
      aux->num_brunsli_histograms = num_histograms;
    }
    return true;
  }
  s->context_map.resize(num_contexts * kNumAvrgContexts);
  if (!DecodeContextMap(
          num_histograms, s->context_map.size(), &s->context_map[0], &br)) {
    return false;
  }

  s->entropy_codes.resize(num_histograms);
  for (int i = 0; i < num_histograms; ++i) {
    if (!s->entropy_codes[i].ReadFromBitStream(kCoeffAlphabetSize, &br)) {
      return false;
    }
  }

  return (BrunsliBitReaderJumpToByteBoundary(&br) == 0);
}

bool DecodeDCDataSection(const uint8_t* data, const size_t len,
                         JPEGDecodingState* s,
                         JPEGData* jpg) {
  if (jpg->width == 0 ||
      jpg->height == 0 ||
      jpg->MCU_rows == 0 ||
      jpg->MCU_cols == 0 ||
      jpg->components.empty() ||
      jpg->quant.empty() ||
      s->context_map.empty()) {
    // DC data can not be decoded without knowing the width and the height
    // and the number of components from the header, the quantization tables
    // from the quant data section and the context map from the histogram data
    // section.
    return false;
  }
  coeff_t* coeffs[kMaxComponents] = { NULL };
  int h_samp[kMaxComponents] = { 0 };
  int v_samp[kMaxComponents] = { 0 };
  for (size_t i = 0; i < jpg->components.size(); ++i) {
    JPEGComponent* c = &jpg->components[i];
    coeffs[i] = &c->coeffs[0];
    h_samp[i] = c->h_samp_factor;
    v_samp[i] = c->v_samp_factor;
  }
  BrunsliInput in(data, len);
  if (!DecodeDC(jpg->MCU_cols, jpg->MCU_rows, jpg->components.size(),
                h_samp, v_samp, s->context_map, s->entropy_codes, coeffs,
                &s->block_state, &in)) {
    return false;
  }
  return (in.len_ == in.pos_);
}

bool DecodeACDataSection(const uint8_t* data, const size_t len,
                         JPEGDecodingState* s,
                         JPEGData* jpg) {
  if (jpg->width == 0 ||
      jpg->height == 0 ||
      jpg->MCU_rows == 0 ||
      jpg->MCU_cols == 0 ||
      jpg->components.empty() ||
      jpg->quant.empty() ||
      s->block_state.empty() ||
      s->context_map.empty()) {
    // AC data can not be decoded without knowing the width and the height
    // and the number of components from the header, the quantization tables
    // from the quant data section, the context map from the histogram data
    // section and the block states from the DC data section.
    return false;
  }
  coeff_t* coeffs[kMaxComponents] = { NULL };
  int quant[kMaxComponents][kDCTBlockSize];
  int h_samp[kMaxComponents] = { 0 };
  int v_samp[kMaxComponents] = { 0 };
  for (size_t i = 0; i < jpg->components.size(); ++i) {
    JPEGComponent* c = &jpg->components[i];
    if (c->quant_idx >= jpg->quant.size()) {
      return false;
    }
    const JPEGQuantTable& q = jpg->quant[c->quant_idx];
    coeffs[i] = &c->coeffs[0];
    memcpy(&quant[i][0], &q.values[0], kDCTBlockSize * sizeof(quant[0][0]));
    for (int k = 0; k < kDCTBlockSize; ++k) {
      if (quant[i][k] == 0) {
        return false;
      }
    }
    h_samp[i] = c->h_samp_factor;
    v_samp[i] = c->v_samp_factor;
  }
  BrunsliInput in(data, len);
  if (!DecodeAC(jpg->MCU_cols, jpg->MCU_rows, jpg->components.size(),
                h_samp, v_samp, quant, s->context_bits, s->context_map,
                s->entropy_codes, s->block_state, coeffs, &in)) {
    return false;
  }
  return true;
}

static BrunsliStatus DecodeOriginalJpg(const uint8_t* data, const size_t len,
                                       size_t* pos, JPEGData* jpg) {
  if (*pos >= len) {
    return BRUNSLI_INVALID_BRN;
  }
  if (data[(*pos)++] != SectionMarker(kBrunsliOriginalJpgTag)) {
    return BRUNSLI_INVALID_BRN;
  }
  size_t jpg_len = 0;
  if (!DecodeDataLength(data, len, pos, &jpg_len)) {
    return BRUNSLI_INVALID_BRN;
  }
  jpg->original_jpg = data + *pos;
  jpg->original_jpg_size = jpg_len;
  *pos += jpg_len;
  if (*pos != len) {
    return BRUNSLI_INVALID_BRN;
  }
  return BRUNSLI_OK;
}

BrunsliStatus BrunsliDecodeJpeg(const uint8_t* data, const size_t len,
                                BrunsliReadMode mode,
                                JPEGData* jpg,
                                BrunsliAuxData* aux) {
  size_t pos = 0;
  BrunsliStatus status;
  status = VerifySignature(data, len, &pos);
  if (status != BRUNSLI_OK) return status;

  status = DecodeHeader(data, len, &pos, jpg);
  if (status != BRUNSLI_OK || mode == BRUNSLI_READ_HEADER) return status;

  bool tags_met[16] = { false };
  // Signature and header can appear only at the start of the data.
  for (size_t tag : {kBrunsliSignatureTag, kBrunsliHeaderTag}) {
    tags_met[tag] = true;
  }
  bool known_section_tags[16] = { false };
  for (size_t tag :
       {kBrunsliSignatureTag, kBrunsliHeaderTag, kBrunsliMetaDataTag,
        kBrunsliJPEGInternalsTag, kBrunsliQuantDataTag,
        kBrunsliHistogramDataTag, kBrunsliDCDataTag, kBrunsliACDataTag,
        kBrunsliOriginalJpgTag}) {
    known_section_tags[tag] = true;
  }

  if (jpg->version == 1) {
    return DecodeOriginalJpg(data, len, &pos, jpg);
  }
  // Do not allow "original_jpg" for regular Brunsli files.
  tags_met[kBrunsliOriginalJpgTag] = true;

  // Allocate coefficient buffers.
  if (mode == BRUNSLI_READ_ALL) {
    for (size_t i = 0; i < jpg->components.size(); ++i) {
      JPEGComponent* c = &jpg->components[i];
      c->coeffs.resize(c->num_blocks * kDCTBlockSize);
    }
  }

  JPEGDecodingState s;
  while (pos < len) {
    const uint8_t marker = data[pos++];
    const size_t tag = marker >> 3;
    if (tag == 0 || tag > 15) return BRUNSLI_INVALID_BRN;
    const size_t wiring_type = marker & 0x7;
    if (wiring_type != kBrunsliWiringTypeVarint &&
        wiring_type != kBrunsliWiringTypeLengthDelimited) {
      return BRUNSLI_INVALID_BRN;
    }

    if (tags_met[tag]) {
      BRUNSLI_LOG_ERROR() << "Duplicate marker " << std::hex
                          << static_cast<int>(marker) << BRUNSLI_ENDL();
      return BRUNSLI_INVALID_BRN;
    }
    tags_met[tag] = true;

    size_t value = 0;
    if (!DecodeDataLength(data, len, &pos, &value)) {
      return BRUNSLI_INVALID_BRN;
    }

    const bool is_section = (wiring_type == kBrunsliWiringTypeLengthDelimited);
    if (known_section_tags[tag] && !is_section) return BRUNSLI_INVALID_BRN;

    // No varint tags on top level.
    if (!is_section) continue;
    const size_t marker_len = value;

    if (mode == BRUNSLI_READ_SIZES && tag != kBrunsliHistogramDataTag) {
      pos += marker_len;
      continue;
    }

    bool ok = false;
    switch (tag) {
      case kBrunsliMetaDataTag:
        ok = DecodeMetaDataSection(&data[pos], marker_len, jpg);
        break;
      case kBrunsliJPEGInternalsTag:
        ok = DecodeJPEGInternalsSection(&data[pos], marker_len, jpg);
        break;
      case kBrunsliQuantDataTag:
        ok = DecodeQuantDataSection(&data[pos], marker_len, jpg);
        break;
      case kBrunsliHistogramDataTag:
        ok = DecodeHistogramDataSection(&data[pos], marker_len, mode, &s, jpg,
                                        aux);
        break;
      case kBrunsliDCDataTag:
        ok = tags_met[kBrunsliQuantDataTag] &&
             DecodeDCDataSection(&data[pos], marker_len, &s, jpg);
        break;
      case kBrunsliACDataTag:
        ok = tags_met[kBrunsliDCDataTag] &&
             DecodeACDataSection(&data[pos], marker_len, &s, jpg);
        break;
      default:
        // Skip unrecognized sections.
        ok = true;
        break;
    }
    if (!ok) {
      return BRUNSLI_INVALID_BRN;
    }
    if (mode == BRUNSLI_READ_SIZES && tag == kBrunsliHistogramDataTag) {
      return BRUNSLI_OK;
    }
    pos += marker_len;
  }
  // It is expected that there is no garbage after the valid brunsli stream.
  if (pos != len) {
    return BRUNSLI_INVALID_BRN;
  }
  return BRUNSLI_OK;
}

}  // namespace brunsli
