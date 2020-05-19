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
#include <brunsli/jpeg_data.h>
#include "../common/lehmer_code.h"
#include "../common/platform.h"
#include "../common/predict.h"
#include "../common/quant_matrix.h"
#include <brunsli/status.h>
#include <brunsli/types.h>
#include "./ans_decode.h"
#include "./arith_decode.h"
#include "./bit_reader.h"
#include "./brunsli_input.h"
#include "./context_map_decode.h"
#include <brunsli/jpeg_data_writer.h>
#include "./state.h"
#include "./state_internal.h"

namespace brunsli {

using ::brunsli::internal::dec::AcDcState;
using ::brunsli::internal::dec::BlockI32;
using ::brunsli::internal::dec::ComponentMeta;
using ::brunsli::internal::dec::FallbackState;
using ::brunsli::internal::dec::HeaderState;
using ::brunsli::internal::dec::InternalState;
using ::brunsli::internal::dec::MetadataDecompressionStage;
using ::brunsli::internal::dec::MetadataState;
using ::brunsli::internal::dec::PrepareMeta;
using ::brunsli::internal::dec::SectionHeaderState;
using ::brunsli::internal::dec::SectionState;
using ::brunsli::internal::dec::Stage;
using ::brunsli::internal::dec::State;
using ::brunsli::internal::dec::UpdateSubsamplingDerivatives;

static const int kNumDirectCodes = 8;
static const int kCoeffAlphabetSize = kNumDirectCodes + 10;

static const uint32_t kKnownSectionTags =
    (1u << kBrunsliSignatureTag) | (1u << kBrunsliHeaderTag) |
    (1u << kBrunsliMetaDataTag) | (1u << kBrunsliJPEGInternalsTag) |
    (1u << kBrunsliQuantDataTag) | (1u << kBrunsliHistogramDataTag) |
    (1u << kBrunsliDCDataTag) | (1u << kBrunsliACDataTag) |
    (1u << kBrunsliOriginalJpgTag);

static const uint32_t kKnownHeaderVarintTags =
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
  if (BrunsliBitReaderRead(br, 1)) {
    uint32_t nbits = BrunsliBitReaderRead(br, 3);
    if (nbits == 0) {
      return 1u;
    } else {
      return BrunsliBitReaderRead(br, nbits) + (1u << nbits);
    }
  }
  return 0;
}

uint32_t DecodeVarint(BrunsliBitReader* br, size_t max_bits) {
  uint32_t result = 0;
  for (size_t b = 0; b < max_bits; ++b) {
    if (b + 1 != max_bits && !BrunsliBitReaderRead(br, 1)) {
      break;
    }
    result |= BrunsliBitReaderRead(br, 1) << b;
  }
  return result;
}

template<size_t kChunkSize>
size_t DecodeLimitedVarint(BrunsliBitReader* br, size_t max_symbols) {
  size_t bits = 0;
  size_t shift = 0;
  for (size_t b = 0; b < max_symbols && BrunsliBitReaderRead(br, 1); ++b) {
    const size_t next_bits = BrunsliBitReaderRead(br, kChunkSize);
    bits |= next_bits << shift;
    shift += kChunkSize;
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

bool ProcessMetaData(const uint8_t* data, size_t len, MetadataState* state,
                     JPEGData* jpg) {
  size_t pos = 0;
  while (pos < len) {
    switch (state->stage) {
      case MetadataState::READ_MARKER: {
        state->marker = static_cast<uint8_t>(data[pos++]);
        if (state->marker == 0xD9) {
          jpg->tail_data = std::string();
          state->stage = MetadataState::READ_TAIL;
          continue;
        } else if (state->marker < 0x40) {
          state->short_marker_count++;
          if (state->short_marker_count > kBrunsliShortMarkerLimit) {
            return false;
          }
          jpg->app_data.push_back(GenerateApp0Marker(state->marker));
          continue;
        } else if (state->marker >= 0x80 && state->marker <= 0x82) {
          state->short_marker_count++;
          if (state->short_marker_count > kBrunsliShortMarkerLimit) {
            return false;
          }
          state->stage = MetadataState::READ_CODE;
          continue;
        }
        // Otherwise - mutlibyte sequence.
        if ((state->marker != 0xFE) && ((state->marker >> 4u) != 0x0E)) {
          return false;
        }
        state->stage = MetadataState::READ_LENGTH_HI;
        continue;
      }

      case MetadataState::READ_TAIL: {
        jpg->tail_data.append(data + pos, data + len);
        pos = len;
        continue;
      }

      case MetadataState::READ_CODE: {
        const uint8_t code = data[pos++];
        jpg->app_data.push_back(GenerateAppMarker(state->marker, code));
        state->stage = MetadataState::READ_MARKER;
        continue;
      }

      case MetadataState::READ_LENGTH_HI: {
        state->length_hi = data[pos++];
        state->stage = MetadataState::READ_LENGTH_LO;
        continue;
      }

      case MetadataState::READ_LENGTH_LO: {
        const uint8_t lo = data[pos++];
        size_t marker_len = (state->length_hi << 8u) + lo;
        if (marker_len < 2) return false;
        state->remaining_multibyte_length = marker_len - 2;
        uint8_t head[3] = {state->marker, state->length_hi, lo};
        std::vector<std::string>* dest =
            (state->marker == 0xFE) ? &jpg->com_data : &jpg->app_data;
        dest->emplace_back(head, head + 3);
        state->multibyte_sink = &dest->back();
        // Turn state machine to default state in case there is no payload in
        // multibyte sequence. This is important when such a sequence concludes
        // the input.
        state->stage = (state->remaining_multibyte_length > 0)
                           ? MetadataState::READ_MULTIBYTE
                           : MetadataState::READ_MARKER;
        continue;
      }

      case MetadataState::READ_MULTIBYTE: {
        size_t chunk_size =
            std::min(state->remaining_multibyte_length, len - pos);
        state->multibyte_sink->append(data + pos, data + pos + chunk_size);
        state->remaining_multibyte_length -= chunk_size;
        pos += chunk_size;
        if (state->remaining_multibyte_length == 0) {
          state->stage = MetadataState::READ_MARKER;
        }
        continue;
      }

      default: return false;
    }
  }
  return true;
}

bool DecodeQuantTables(BrunsliBitReader* br, JPEGData* jpg) {
  bool have_internals_data = !jpg->quant.empty();
  size_t num_quant_tables = BrunsliBitReaderRead(br, 2) + 1;
  if (jpg->quant.size() != num_quant_tables) {
    return false;
  }
  for (size_t i = 0; i < num_quant_tables; ++i) {
    JPEGQuantTable* q = &jpg->quant[i];
    int data_precision = 0;
    if (!BrunsliBitReaderRead(br, 1)) {
      const size_t short_code = BrunsliBitReaderRead(br, 3);
      for (size_t k = 0; k < kDCTBlockSize; ++k) {
        q->values[k] = kStockQuantizationTables[(i > 0) ? 1 : 0][short_code][k];
      }
    } else {
      const uint32_t q_factor = BrunsliBitReaderRead(br, 6);
      uint8_t predictor[kDCTBlockSize];
      FillQuantMatrix(i > 0, q_factor, predictor);
      int delta = 0;
      for (size_t k = 0; k < kDCTBlockSize; ++k) {
        if (BrunsliBitReaderRead(br, 1)) {
          const int sign = BrunsliBitReaderRead(br, 1);
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
    c->quant_idx = BrunsliBitReaderRead(br, 2);
    if (c->quant_idx >= jpg->quant.size()) {
      return false;
    }
  }
  return BrunsliBitReaderIsHealthy(br);
}

bool DecodeHuffmanCode(BrunsliBitReader* br, JPEGHuffmanCode* huff,
                       bool is_known_last) {
  huff->slot_id = BrunsliBitReaderRead(br, 2);
  int is_dc_table = (BrunsliBitReaderRead(br, 1) == 0);
  huff->slot_id += is_dc_table ? 0 : 0x10;
  huff->is_last = is_known_last || BrunsliBitReaderRead(br, 1);
  huff->counts[0] = 0;
  int found_match = BrunsliBitReaderRead(br, 1);
  if (found_match) {
    if (is_dc_table) {
      int huff_table_idx = BrunsliBitReaderRead(br, 1);
      memcpy(&huff->counts[1], kStockDCHuffmanCodeCounts[huff_table_idx],
             sizeof(kStockDCHuffmanCodeCounts[0]));
      memcpy(&huff->values[0], kStockDCHuffmanCodeValues[huff_table_idx],
             sizeof(kStockDCHuffmanCodeValues[0]));
    } else {
      int huff_table_idx = BrunsliBitReaderRead(br, 1);
      memcpy(&huff->counts[1], kStockACHuffmanCodeCounts[huff_table_idx],
             sizeof(kStockACHuffmanCodeCounts[0]));
      memcpy(&huff->values[0], kStockACHuffmanCodeValues[huff_table_idx],
             sizeof(kStockACHuffmanCodeValues[0]));
    }
    return BrunsliBitReaderIsHealthy(br);
  }
  int total_count = 0;
  int space = 1u << kJpegHuffmanMaxBitLength;
  int max_len = BrunsliBitReaderRead(br, 4) + 1;
  int max_count = is_dc_table ? kJpegDCAlphabetSize : kJpegHuffmanAlphabetSize;
  space -= 1u << (kJpegHuffmanMaxBitLength - max_len);
  for (int i = 1; i <= max_len; ++i) {
    size_t shift = kJpegHuffmanMaxBitLength - i;
    int count_limit = std::min(max_count - total_count, space >> shift);
    if (count_limit > 0) {
      int nbits = Log2FloorNonZero(count_limit) + 1;
      int count = BrunsliBitReaderRead(br, nbits);
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
    const int nbits = p.num_bits();
    const size_t code = DecodeLimitedVarint<2>(br, (nbits + 1) >> 1u);
    uint8_t value;
    if (!p.Remove(code, &value)) {
      return false;
    }
    huff->values[i] = value;
  }
  huff->values[total_count] = kJpegHuffmanAlphabetSize;
  return BrunsliBitReaderIsHealthy(br);
}

bool DecodeScanInfo(BrunsliBitReader* br, JPEGScanInfo* si) {
  si->Ss = BrunsliBitReaderRead(br, 6);
  si->Se = BrunsliBitReaderRead(br, 6);
  si->Ah = BrunsliBitReaderRead(br, 4);
  si->Al = BrunsliBitReaderRead(br, 4);
  si->components.resize(BrunsliBitReaderRead(br, 2) + 1);
  for (size_t i = 0; i < si->components.size(); ++i) {
    si->components[i].comp_idx = BrunsliBitReaderRead(br, 2);
    si->components[i].dc_tbl_idx = BrunsliBitReaderRead(br, 2);
    si->components[i].ac_tbl_idx = BrunsliBitReaderRead(br, 2);
  }
  int last_block_idx = -1;
  while (BrunsliBitReaderRead(br, 1)) {
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
  while (BrunsliBitReaderRead(br, 1)) {
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
  return BrunsliBitReaderIsHealthy(br);
}

bool DecodeAuxData(BrunsliBitReader* br, JPEGData* jpg) {
  bool have_dri = false;
  int num_scans = 0;
  size_t dht_count = 0;
  uint8_t marker;
  do {
    if (!BrunsliBitReaderIsHealthy(br)) return false;
    marker = 0xc0 + BrunsliBitReaderRead(br, 6);
    jpg->marker_order.push_back(marker);
    if (marker == 0xc4) ++dht_count;
    if (marker == 0xdd) have_dri = true;
    if (marker == 0xda) ++num_scans;
  } while (marker != 0xd9);
  if (have_dri) {
    jpg->restart_interval = BrunsliBitReaderRead(br, 16);
  }

  size_t terminal_huffman_code_count = 0;
  for (int i = 0;; ++i) {
    const bool is_known_last = BrunsliBitReaderRead(br, 1);
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
  int num_quant_tables = BrunsliBitReaderRead(br, 2) + 1;
  jpg->quant.resize(num_quant_tables);
  for (int i = 0; i < num_quant_tables; ++i) {
    JPEGQuantTable* q = &jpg->quant[i];
    q->index = BrunsliBitReaderRead(br, 2);
    q->is_last = (i == num_quant_tables - 1) || BrunsliBitReaderRead(br, 1);
    q->precision = BrunsliBitReaderRead(br, 4);
    if (q->precision > 1) {
      BRUNSLI_LOG_ERROR() << "Invalid quantization table precision: "
                          << q->precision << BRUNSLI_ENDL();
      return false;
    }
    // note that q->values[] are initialized to invalid 0 values.
  }
  int comp_ids = BrunsliBitReaderRead(br, 2);
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
      jpg->components[i].id = BrunsliBitReaderRead(br, 8);
    }
  }

  // security: limit is 32b for n_size
  size_t n_size = DecodeLimitedVarint<8>(br, 4);
  jpg->has_zero_padding_bit = (n_size > 0);
  if (n_size > 0) {
    if (n_size > PaddingBitsLimit(*jpg)) {
      BRUNSLI_LOG_ERROR() << "Suspicious number of padding bits " << n_size
                          << BRUNSLI_ENDL();
      return false;
    }
    jpg->padding_bits.resize(n_size);
    for (size_t i = 0; i < n_size; ++i) {
      jpg->padding_bits[i] = BrunsliBitReaderRead(br, 1);
    }
  }
  return BrunsliBitReaderIsHealthy(br);
}

static bool BRUNSLI_NOINLINE DecodeCoeffOrder(uint32_t* order, BitSource* br,
                                              WordSource* in) {
  uint32_t lehmer[kDCTBlockSize] = {0};
  static const int kSpan = 16;
  for (int i = 0; i < kDCTBlockSize; i += kSpan) {
    if (!br->ReadBits(1, in)) continue;  // span is all-zero
    const int start = (i > 0) ? i : 1;
    const int end = i + kSpan;
    for (int j = start; j < end; ++j) {
      uint32_t v = 0;
      while (v <= kDCTBlockSize) {
        const uint32_t bits = br->ReadBits(3, in);
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

static size_t DecodeNumNonzeros(Prob* p, BinaryArithmeticDecoder* ac,
                                WordSource* in) {
  // To simplity BST navigation, we use 1-based indexing.
  Prob* bst = p - 1;
  size_t ctx = 1;

  for (size_t b = 0; b < kNumNonZeroBits; ++b) {
    const int bit = ac->ReadBit(bst[ctx].get_proba(), in);
    bst[ctx].Add(bit);
    ctx = 2 * ctx + bit;
  }

  // Leaf index in the level corresponds to the resuling value.
  size_t val = ctx - (1u << kNumNonZeroBits);
  BRUNSLI_DCHECK(val <= kNumNonZeroTreeSize);
  return val;
}

void EnsureSubdecodersInitialized(State* state, WordSource* in) {
  InternalState& s = *state->internal;
  if (!s.subdecoders_initialized) {
    s.ans_decoder.Init(in);
    s.bit_reader.Init(in);
    s.arith_decoder.Init(in);
    s.subdecoders_initialized = true;
  }
}

bool FinalizeSubdecoders(State* state) {
  InternalState& s = *state->internal;
  if (!s.ans_decoder.CheckCRC()) return false;
  if (!s.bit_reader.Finish()) return false;
  s.subdecoders_initialized = false;
  return true;
}

bool DecodeDC(State* state, WordSource* in) {
  const std::vector<ComponentMeta>& meta = state->meta;
  const size_t num_components = meta.size();
  const int mcu_rows = meta[0].height_in_blocks / meta[0].v_samp;
  InternalState& s = *state->internal;
  AcDcState& ac_dc_state = s.ac_dc;

  std::vector<ComponentStateDC>& comps = ac_dc_state.dc;
  if (comps.empty()) {
    comps.resize(num_components);
    for (size_t c = 0; c < num_components; ++c) {
      comps[c].SetWidth(meta[c].width_in_blocks);
    }
  }

  EnsureSubdecodersInitialized(state, in);
  ANSDecoder ans = s.ans_decoder;
  BitSource br = s.bit_reader;
  BinaryArithmeticDecoder ac = s.arith_decoder;

  // We decode DC components in the following interleaved manner:
  //   v_samp[0] rows from component 0
  //   v_samp[1] rows from component 1
  //   v_samp[2] rows from component 2
  //   v_samp[3] rows from component 3 (if present)
  //
  // E.g. in a YUV420 image, we decode 2 rows of DC components from Y and then
  // 1 row of DC components from U and 1 row of DC components from V.
  for (int mcu_y = ac_dc_state.next_mcu_y; mcu_y < mcu_rows; ++mcu_y) {
    for (size_t i = ac_dc_state.next_component; i < num_components; ++i) {
      ComponentStateDC* c = &comps[i];
      const ComponentMeta& m = meta[i];
      const uint8_t* context_map = state->context_map + i * kNumAvrgContexts;
      const size_t ac_stride = m.ac_stride;
      const size_t b_stride = m.b_stride;
      const int width = m.width_in_blocks;
      int y = mcu_y * m.v_samp + ac_dc_state.next_iy;
      int* const prev_sgn = &c->prev_sign[1];
      int* const prev_abs = &c->prev_abs_coeff[2];
      for (int iy = ac_dc_state.next_iy; iy < m.v_samp; ++iy, ++y) {
        coeff_t* coeffs =
            m.ac_coeffs + y * ac_stride + ac_dc_state.next_x * kDCTBlockSize;
        uint8_t* block_state =
            m.block_state + y * b_stride + ac_dc_state.next_x;
        for (int x = ac_dc_state.next_x; x < width; ++x) {
          // Each iteration uses up to 4 words.
          const int is_empty_ctx =
              IsEmptyBlockContext(&c->prev_is_nonempty[1], x);
          Prob* BRUNSLI_RESTRICT is_empty_p =
              &c->is_empty_block_prob[is_empty_ctx];
          const bool is_empty_block = !ac.ReadBit(is_empty_p->get_proba(), in);
          is_empty_p->Add(!is_empty_block);
          c->prev_is_nonempty[x + 1] = !is_empty_block;
          *block_state = is_empty_block;
          int abs_val = 0;
          int sign = 0;
          if (!is_empty_block) {
            Prob* BRUNSLI_RESTRICT p_is_zero = &c->is_zero_prob;
            int is_zero = ac.ReadBit(p_is_zero->get_proba(), in);
            p_is_zero->Add(is_zero);
            if (!is_zero) {
              const int avg_ctx = WeightedAverageContextDC(prev_abs, x);
              const int sign_ctx = prev_sgn[x] * 3 + prev_sgn[x - 1];
              Prob* BRUNSLI_RESTRICT sign_p = &c->sign_prob[sign_ctx];
              sign = ac.ReadBit(sign_p->get_proba(), in);
              sign_p->Add(sign);
              const int entropy_ix = context_map[avg_ctx];
              int code = ans.ReadSymbol(state->entropy_codes[entropy_ix], in);
              if (code < kNumDirectCodes) {
                abs_val = code + 1;
              } else {
                const size_t nbits = code - kNumDirectCodes;
                Prob* BRUNSLI_RESTRICT p_first_extra_bit =
                    &c->first_extra_bit_prob[nbits];
                int first_extra_bit =
                    ac.ReadBit(p_first_extra_bit->get_proba(), in);
                p_first_extra_bit->Add(first_extra_bit);
                int extra_bits_val = first_extra_bit << nbits;
                if (nbits > 0) {
                  extra_bits_val |= br.ReadBits(nbits, in);
                }
                abs_val = kNumDirectCodes - 1 + (2 << nbits) + extra_bits_val;
              }
            }
          }
          prev_abs[x] = abs_val;
          prev_sgn[x] = abs_val ? sign + 1 : 0;
          coeffs[0] = ((1 - 2 * sign) * abs_val +
                       PredictWithAdaptiveMedian(coeffs, x, y, ac_stride));
          block_state++;
          coeffs += kDCTBlockSize;
        }
        ac_dc_state.next_x = 0;
      }
      ac_dc_state.next_iy = 0;
    }
    ac_dc_state.next_component = 0;
  }
  // ac_dc_state.next_mcu_y = 0;

  // Prepare for AC decoding.
  ac_dc_state.next_mcu_y = 0;
  ac_dc_state.next_component = 0;
  ac_dc_state.next_iy = 0;
  ac_dc_state.next_x = 0;

  comps.clear();
  comps.shrink_to_fit();

  s.ans_decoder = ans;
  s.bit_reader = br;
  s.arith_decoder = ac;
  if (!FinalizeSubdecoders(state)) return false;

  return (in->error_ == 0);
}

static void BRUNSLI_NOINLINE DecodeEmptyAcBlock(
    int* BRUNSLI_RESTRICT prev_sgn, int* BRUNSLI_RESTRICT prev_abs) {
  for (int k = 1; k < kDCTBlockSize; ++k) {
    prev_sgn[k] = 0;
    prev_abs[k] = 0;
  }
}

/** All the necessary things for decoding AC block. */
struct AcBlockCookie {
  int x;
  int y;
  uint8_t* BRUNSLI_RESTRICT prev_num_nonzeros;
  int* BRUNSLI_RESTRICT prev_sgn;
  int* BRUNSLI_RESTRICT prev_abs;
  Prob* BRUNSLI_RESTRICT num_nonzero_prob;
  BinaryArithmeticDecoder* BRUNSLI_RESTRICT ac;
  WordSource* BRUNSLI_RESTRICT in;
  ANSDecoder* BRUNSLI_RESTRICT ans;
  BitSource* BRUNSLI_RESTRICT br;
  coeff_t* BRUNSLI_RESTRICT coeffs;
  const coeff_t* BRUNSLI_RESTRICT prev_row_coeffs;
  const coeff_t* BRUNSLI_RESTRICT prev_col_coeffs;
  Prob* BRUNSLI_RESTRICT is_zero_prob;
  const uint32_t* BRUNSLI_RESTRICT order;
  const int* BRUNSLI_RESTRICT mult_col;
  const int* BRUNSLI_RESTRICT mult_row;
  int prev_row_delta;
  Prob* BRUNSLI_RESTRICT sign_prob;
  int context_bits;
  const uint8_t* BRUNSLI_RESTRICT context_map;
  const ANSDecodingData* BRUNSLI_RESTRICT entropy_codes;
  Prob* BRUNSLI_RESTRICT first_extra_bit_prob;
};

static size_t BRUNSLI_NOINLINE DecodeAcBlock(const AcBlockCookie& cookie) {
  AcBlockCookie c = cookie;
  ANSDecoder ans = *c.ans;
  BitSource br = *c.br;
  BinaryArithmeticDecoder ac = *c.ac;

  // Each iteration uses up to 179 words.
  size_t num_nonzeros = 0;

  const uint8_t nonzero_ctx = NumNonzerosContext(c.prev_num_nonzeros, c.x, c.y);
  size_t last_nz = DecodeNumNonzeros(
      c.num_nonzero_prob + kNumNonZeroTreeSize * nonzero_ctx, &ac, c.in);
  for (size_t k = last_nz + 1; k < kDCTBlockSize; ++k) {
    c.prev_sgn[k] = 0;
    c.prev_abs[k] = 0;
  }
  for (size_t k = last_nz; k > 0; --k) {
    int is_zero = 0;
    if (k < last_nz) {
      const int bucket = kNonzeroBuckets[num_nonzeros - 1];
      const int is_zero_ctx = bucket * kDCTBlockSize + k;
      Prob& p = c.is_zero_prob[is_zero_ctx];
      is_zero = ac.ReadBit(p.get_proba(), c.in);
      p.Add(is_zero);
    }
    int abs_val = 0;
    int sign = 1;
    const int k_nat = c.order[k];
    if (!is_zero) {
      int avg_ctx = 0;
      int sign_ctx = kMaxAverageContext;
      if (k_nat < 8) {
        if (c.y > 0) {
          ACPredictContextRow(c.prev_row_coeffs + k_nat, c.coeffs + k_nat,
                              c.mult_col + k_nat * 8, &avg_ctx, &sign_ctx);
        }
      } else if ((k_nat & 7u) == 0) {
        if (c.x > 0) {
          ACPredictContextCol(c.prev_col_coeffs + k_nat, c.coeffs + k_nat,
                              c.mult_row + k_nat, &avg_ctx, &sign_ctx);
        }
      } else {
        avg_ctx = WeightedAverageContext(c.prev_abs + k, c.prev_row_delta);
        sign_ctx =
            c.prev_sgn[k] * 3 + c.prev_sgn[static_cast<int>(k) - kDCTBlockSize];
      }
      sign_ctx = sign_ctx * kDCTBlockSize + k;
      Prob& sign_p = c.sign_prob[sign_ctx];
      sign = ac.ReadBit(sign_p.get_proba(), c.in);
      sign_p.Add(sign);
      c.prev_sgn[k] = sign + 1;
      sign = 1 - 2 * sign;
      const int z_dens_ctx =
          ZeroDensityContext(num_nonzeros, k, c.context_bits);
      const int histo_ix = z_dens_ctx * kNumAvrgContexts + avg_ctx;
      const int entropy_ix = c.context_map[histo_ix];
      int code = ans.ReadSymbol(c.entropy_codes[entropy_ix], c.in);
      if (code < kNumDirectCodes) {
        abs_val = code + 1;
      } else {
        int nbits = code - kNumDirectCodes;
        Prob& p = c.first_extra_bit_prob[k * 10 + nbits];
        int first_extra_bit = ac.ReadBit(p.get_proba(), c.in);
        p.Add(first_extra_bit);
        int extra_bits_val = first_extra_bit << nbits;
        if (nbits > 0) {
          extra_bits_val |= br.ReadBits(nbits, c.in);
        }
        abs_val = kNumDirectCodes - 1 + (2u << nbits) + extra_bits_val;
      }
      ++num_nonzeros;
    } else {
      c.prev_sgn[k] = 0;
    }
    int coeff = sign * abs_val;
    c.coeffs[k_nat] = coeff;
    c.prev_abs[k] = abs_val;
  }

  *c.ans = ans;
  *c.br = br;
  *c.ac = ac;

  return num_nonzeros;
}

bool DecodeAC(State* state, WordSource* in) {
  const std::vector<ComponentMeta>& meta = state->meta;
  const size_t num_components = meta.size();
  const int mcu_rows = meta[0].height_in_blocks / meta[0].v_samp;
  InternalState& s = *state->internal;
  AcDcState& ac_dc_state = s.ac_dc;

  std::vector<ComponentState>& comps = ac_dc_state.ac;
  if (comps.empty()) {
    comps.resize(num_components);
    for (size_t c = 0; c < num_components; ++c) {
      comps[c].SetWidth(meta[c].width_in_blocks);
      ComputeACPredictMultipliers(&meta[c].quant[0], comps[c].mult_row,
                                  comps[c].mult_col);
    }
  }

  EnsureSubdecodersInitialized(state, in);

  if (!ac_dc_state.ac_coeffs_order_decoded) {
    while (ac_dc_state.next_component < num_components) {
      // Uses up to 121 word.
      if (!DecodeCoeffOrder(comps[ac_dc_state.next_component].order,
                            &s.bit_reader, in)) {
        return false;
      }
      ac_dc_state.next_component++;
    }
    ac_dc_state.next_component = 0;
    ac_dc_state.ac_coeffs_order_decoded = true;
  }

  AcBlockCookie c;
  c.ac = &s.arith_decoder;
  c.in = in;
  c.ans = &s.ans_decoder;
  c.br = &s.bit_reader;
  c.entropy_codes = state->entropy_codes;

  for (int mcu_y = ac_dc_state.next_mcu_y; mcu_y < mcu_rows; ++mcu_y) {
    for (size_t i = ac_dc_state.next_component; i < num_components; ++i) {
      ComponentState& cst = comps[i];
      c.prev_num_nonzeros = cst.prev_num_nonzeros.data();
      c.num_nonzero_prob = cst.num_nonzero_prob;
      c.is_zero_prob = cst.is_zero_prob.data();
      c.order = cst.order;
      c.mult_col = cst.mult_col;
      c.mult_row = cst.mult_row;
      c.sign_prob = cst.sign_prob.data();
      c.first_extra_bit_prob = cst.first_extra_bit_prob.data();
      const ComponentMeta& m = meta[i];
      c.context_map = state->context_map + m.context_offset * kNumAvrgContexts;
      c.context_bits = m.context_bits;
      const int width = m.width_in_blocks;
      const size_t ac_stride = m.ac_stride;
      const size_t b_stride = m.b_stride;
      const int next_iy = ac_dc_state.next_iy;
      c.y = mcu_y * m.v_samp + next_iy;
      c.prev_row_delta = (1 - 2 * (c.y & 1u)) * (width + 3) * kDCTBlockSize;
      for (int iy = next_iy; iy < m.v_samp; ++iy, ++c.y) {
        const int next_x = ac_dc_state.next_x;
        const size_t block_offset = next_x * kDCTBlockSize;
        c.coeffs = m.ac_coeffs + c.y * ac_stride + block_offset;
        c.prev_row_coeffs = c.coeffs - ac_stride;
        c.prev_col_coeffs = c.coeffs - kDCTBlockSize;
        const uint8_t* block_state = m.block_state + c.y * b_stride + next_x;
        c.prev_sgn = &cst.prev_sign[kDCTBlockSize] + block_offset;
        c.prev_abs = &cst.prev_abs_coeff[((c.y & 1u) * (width + 3) + 2) *
                                         kDCTBlockSize] +
                     block_offset;
        for (c.x = next_x; c.x < width; ++c.x) {
          bool is_empty = *(block_state++);
          if (!is_empty) {
            size_t num_nonzeros = DecodeAcBlock(c);
            BRUNSLI_DCHECK(num_nonzeros <= kNumNonZeroTreeSize);
            c.prev_num_nonzeros[c.x] = static_cast<uint8_t>(num_nonzeros);
          } else {
            DecodeEmptyAcBlock(c.prev_sgn, c.prev_abs);
            c.prev_num_nonzeros[c.x] = 0;
          }
          c.coeffs += kDCTBlockSize;
          c.prev_sgn += kDCTBlockSize;
          c.prev_abs += kDCTBlockSize;
          c.prev_row_coeffs += kDCTBlockSize;
          c.prev_col_coeffs += kDCTBlockSize;
        }
        c.prev_row_delta *= -1;
        ac_dc_state.next_x = 0;
      }
      ac_dc_state.next_iy = 0;
    }
    ac_dc_state.next_component = 0;
  }
  // ac_dc_state.next_mcu_y = 0;

  comps.clear();
  comps.shrink_to_fit();

  if (!FinalizeSubdecoders(state)) return false;

  return (in->error_ == 0);
}

static bool CheckCanRead(State* state, size_t required) {
  // TODO(eustas): dcheck len > pos
  size_t available = state->len - state->pos;
  return required <= available;
}

static bool CheckCanReadByte(State* state) {
  // TODO(eustas): dcheck len > pos
  return state->pos != state->len;
}

static uint8_t ReadByte(State* state) {
  // TODO(eustas): dcheck len > pos
  return state->data[state->pos++];
}

static uint8_t PeekByte(State* state, size_t offset) {
  // TODO(eustas): dcheck overflow.
  return state->data[state->pos + offset];
}

static void SkipBytes(State* state, size_t len) {
  // TODO(eustas): dcheck overflow.
  state->pos += len;
}

static size_t GetBytesAvailable(State* state) {
  // TODO(eustas): dcheck len > pos
  return state->len - state->pos;
}

static size_t SkipAvailableBytes(State* state, size_t len) {
  size_t available = GetBytesAvailable(state);
  size_t skip_bytes = std::min(available, len);
  state->pos += skip_bytes;
  return skip_bytes;
}

static BrunsliStatus DecodeBase128(State* state, size_t* val) {
  *val = 0;
  uint64_t b = 0x80;
  size_t i = 0;
  while ((i < 9) && (b & 0x80u)) {
    if (!CheckCanRead(state, i + 1)) return BRUNSLI_NOT_ENOUGH_DATA;
    b = PeekByte(state, i);
    *val |= (b & 0x7Fu) << (i * 7);
    ++i;
  }
  SkipBytes(state, i);
  return ((b & 0x80u) == 0) ? BRUNSLI_OK : BRUNSLI_INVALID_BRN;
}

static bool DecodeDataLength(State* state, size_t* data_len) {
  if (DecodeBase128(state, data_len) != BRUNSLI_OK) return false;
  return CheckCanRead(state, *data_len);
}

static Stage Fail(State* state, BrunsliStatus result) {
  InternalState& s = *state->internal;
  s.result = result;
  // Preserve current stage for continuation / error reporting.
  s.last_stage = state->stage;
  return Stage::ERROR;
}

static BrunsliStatus ReadTag(State* state, SectionState* section) {
  if (!CheckCanReadByte(state)) return BRUNSLI_NOT_ENOUGH_DATA;
  const uint8_t marker = ReadByte(state);

  const size_t tag = marker >> 3u;
  if (tag == 0 || tag > 15) return BRUNSLI_INVALID_BRN;
  section->tag = tag;

  const size_t wiring_type = marker & 0x7u;
  if (wiring_type != kBrunsliWiringTypeVarint &&
      wiring_type != kBrunsliWiringTypeLengthDelimited) {
    return BRUNSLI_INVALID_BRN;
  }
  section->is_section = (wiring_type == kBrunsliWiringTypeLengthDelimited);

  const uint32_t tag_bit = 1u << tag;
  if (section->tags_met & tag_bit) {
    BRUNSLI_LOG_ERROR() << "Duplicate marker " << std::hex
                        << static_cast<int>(marker) << BRUNSLI_ENDL();
    return BRUNSLI_INVALID_BRN;
  }
  section->tags_met |= tag_bit;

  return BRUNSLI_OK;
}

static BrunsliStatus EnterSection(State* state, SectionState* section) {
  size_t section_size;
  BrunsliStatus status = DecodeBase128(state, &section_size);
  if (status != BRUNSLI_OK) return status;
  section->is_active = true;
  section->remaining = section_size;
  section->milestone = state->pos;
  section->projected_end = state->pos + section->remaining;
  return BRUNSLI_OK;
}

static void LeaveSection(SectionState* section) {
  section->is_active = false;
}

static bool IsOutOfSectionBounds(State* state) {
  return state->pos > state->internal->section.projected_end;
}

static size_t RemainingSectionLength(State* state) {
  // TODO(eustas): remove this check?
  if (IsOutOfSectionBounds(state)) return 0;
  return state->internal->section.projected_end - state->pos;
}

static bool IsAtSectionBoundary(State* state) {
  return state->pos == state->internal->section.projected_end;
}

Stage VerifySignature(State* state) {
  InternalState& s = *state->internal;

  if (!CheckCanRead(state, kBrunsliSignatureSize)) {
    return Fail(state, BRUNSLI_NOT_ENOUGH_DATA);
  }
  const bool is_signature_ok =
      (memcmp(state->data + state->pos, kBrunsliSignature,
              kBrunsliSignatureSize) != 0);
  state->pos += kBrunsliSignatureSize;
  s.section.tags_met |= 1u << kBrunsliSignatureTag;
  if (is_signature_ok) return Fail(state, BRUNSLI_INVALID_BRN);
  return Stage::HEADER;
}

// Parses the brunsli header starting at data[*pos] and fills in *jpg.
// Sets *pos to the position after the header.
// Returns BRUNSLI_OK, unless the data is not valid brunsli byte stream
// or is truncated.
Stage DecodeHeader(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;
  HeaderState& hs = s.header;

  while (hs.stage != HeaderState::DONE) {
    switch (hs.stage) {
      case HeaderState::READ_TAG: {
        BrunsliStatus status = ReadTag(state, &s.section);
        if (status != BRUNSLI_OK) return Fail(state, status);
        if (s.section.tag != kBrunsliHeaderTag || !s.section.is_section) {
          return Fail(state, BRUNSLI_INVALID_BRN);
        }
        hs.stage = HeaderState::ENTER_SECTION;
        break;
      }

      case HeaderState::ENTER_SECTION: {
        BrunsliStatus status = EnterSection(state, &s.section);
        if (status != BRUNSLI_OK) return Fail(state, status);
        hs.stage = HeaderState::ITEM_READ_TAG;
        break;
      }

      case HeaderState::ITEM_READ_TAG: {
        if (IsAtSectionBoundary(state)) {
          hs.stage = HeaderState::FINALE;
          break;
        }
        BrunsliStatus status = ReadTag(state, &hs.section);
        if (status != BRUNSLI_OK) return Fail(state, status);
        const uint32_t tag_bit = 1u << hs.section.tag;
        if (hs.section.is_section) {
          if (kKnownHeaderVarintTags & tag_bit) {
            Fail(state, BRUNSLI_INVALID_BRN);
          }
          hs.stage = HeaderState::ITEM_ENTER_SECTION;
          break;
        }
        hs.stage = HeaderState::ITEM_READ_VALUE;
        break;
      }

      case HeaderState::ITEM_ENTER_SECTION: {
        BrunsliStatus status = DecodeBase128(state, &hs.remaining_skip_length);
        if (status != BRUNSLI_OK) return Fail(state, status);
        hs.stage = HeaderState::ITEM_SKIP_CONTENTS;
        break;
      }

      case HeaderState::ITEM_SKIP_CONTENTS: {
        size_t bytes_skipped =
            SkipAvailableBytes(state, hs.remaining_skip_length);
        hs.remaining_skip_length -= bytes_skipped;
        if (hs.remaining_skip_length > 0) {
          return Fail(state, BRUNSLI_NOT_ENOUGH_DATA);
        }
        hs.stage = HeaderState::ITEM_READ_TAG;
        break;
      }

      case HeaderState::ITEM_READ_VALUE: {
        size_t value;
        BrunsliStatus status = DecodeBase128(state, &value);
        if (status != BRUNSLI_OK) return Fail(state, status);
        hs.varint_values[hs.section.tag] = value;
        hs.stage = HeaderState::ITEM_READ_TAG;
        break;
      }

      case HeaderState::FINALE: {
        const bool has_version =
            hs.section.tags_met & (1u << kBrunsliHeaderVersionCompTag);
        if (!has_version) return Fail(state, BRUNSLI_INVALID_BRN);
        const size_t version_and_comp_count =
            hs.varint_values[kBrunsliHeaderVersionCompTag];

        const int version = version_and_comp_count >> 2u;
        jpg->version = version;

        if (version == 1) {  // fallback mode
          // TODO(eustas): do we need this?
          jpg->width = 0;
          jpg->height = 0;
          hs.stage = HeaderState::DONE;
          break;
        }

        if (version != 0) {  // unknown mode
          return Fail(state, BRUNSLI_INVALID_BRN);
        }

        // Otherwise version == 0 i.e. regular brunsli.

        // Do not allow "original_jpg" for regular Brunsli files.
        s.section.tags_met |= 1u << kBrunsliOriginalJpgTag;

        const bool has_width =
            hs.section.tags_met & (1u << kBrunsliHeaderWidthTag);
        if (!has_width) return Fail(state, BRUNSLI_INVALID_BRN);
        const size_t width = hs.varint_values[kBrunsliHeaderWidthTag];
        const bool has_height =
            hs.section.tags_met & (1u << kBrunsliHeaderHeightTag);
        if (!has_height) return Fail(state, BRUNSLI_INVALID_BRN);
        const size_t height = hs.varint_values[kBrunsliHeaderHeightTag];

        if (width == 0 || height == 0) return Fail(state, BRUNSLI_INVALID_BRN);
        if (width > kMaxDimPixels || height > kMaxDimPixels) {
          return Fail(state, BRUNSLI_INVALID_BRN);
        }
        jpg->width = width;
        jpg->height = height;

        const size_t num_components = (version_and_comp_count & 3u) + 1u;
        jpg->components.resize(num_components);

        const bool has_subsampling =
            hs.section.tags_met & (1u << kBrunsliHeaderSubsamplingTag);
        if (!has_subsampling) return Fail(state, BRUNSLI_INVALID_BRN);
        size_t subsampling_code =
            hs.varint_values[kBrunsliHeaderSubsamplingTag];

        for (size_t i = 0; i < jpg->components.size(); ++i) {
          JPEGComponent* c = &jpg->components[i];
          c->v_samp_factor = (subsampling_code & 0xFu) + 1;
          subsampling_code >>= 4u;
          c->h_samp_factor = (subsampling_code & 0xFu) + 1;
          subsampling_code >>= 4u;
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

        hs.stage = HeaderState::DONE;
        break;
      }

      default: return Fail(state, BRUNSLI_DECOMPRESSION_ERROR);
    }
  }

  LeaveSection(&s.section);
  return jpg->version == 1 ? Stage::FALLBACK : Stage::SECTION;
}

static BrunsliStatus DecodeMetaDataSection(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;
  MetadataState& ms = s.metadata;

  if (ms.decompression_stage == MetadataDecompressionStage::DONE) {
    return BRUNSLI_INVALID_BRN;
  }

  if (ms.decompression_stage == MetadataDecompressionStage::INITIAL) {
    if (IsAtSectionBoundary(state)) {
      ms.decompression_stage = MetadataDecompressionStage::DONE;
      return BRUNSLI_OK;
    }
    if (RemainingSectionLength(state) == 1) {
      if (!CheckCanReadByte(state)) {
        return BRUNSLI_NOT_ENOUGH_DATA;
      }
      uint8_t data[1];
      data[0] = ReadByte(state);
      bool ok = ProcessMetaData(data, 1, &ms, jpg) && ms.CanFinish();
      ms.decompression_stage = MetadataDecompressionStage::DONE;
      return ok ? BRUNSLI_OK : BRUNSLI_INVALID_BRN;
    }
    ms.decompression_stage =
        MetadataDecompressionStage::READ_LENGTH;
  }

  if (ms.decompression_stage == MetadataDecompressionStage::READ_LENGTH) {
    BrunsliStatus status = DecodeBase128(state, &ms.metadata_size);
    if (status != BRUNSLI_OK) return status;
    // TODO(eustas): ms.metadata_size should be limited to avoid "zip-bombs".
    if (IsOutOfSectionBounds(state)) return BRUNSLI_INVALID_BRN;
    if (RemainingSectionLength(state) == 0) return BRUNSLI_INVALID_BRN;
    ms.brotli = BrotliDecoderCreateInstance(nullptr, nullptr, nullptr);
    if (ms.brotli == nullptr) return BRUNSLI_DECOMPRESSION_ERROR;
    ms.decompression_stage = MetadataDecompressionStage::DECOMPRESSING;
  }

  if (ms.decompression_stage == MetadataDecompressionStage::DECOMPRESSING) {
    // Free Brotli decoder and return result
    const auto finish_decompression = [&ms] (BrunsliStatus result) {
      BRUNSLI_DCHECK(ms.brotli != nullptr);
      BrotliDecoderDestroyInstance(ms.brotli);
      ms.brotli = nullptr;
      ms.decompression_stage = MetadataDecompressionStage::DONE;
      return result;
    };

    while (true) {
      size_t available_bytes =
          std::min(GetBytesAvailable(state), RemainingSectionLength(state));
      size_t available_in = available_bytes;
      const uint8_t* next_in = state->data + state->pos;
      size_t available_out = 0;
      BrotliDecoderResult result = BrotliDecoderDecompressStream(
          ms.brotli, &available_in, &next_in, &available_out, nullptr, nullptr);
      if (result == BROTLI_DECODER_RESULT_ERROR) {
        return finish_decompression(BRUNSLI_INVALID_BRN);
      }
      size_t chunk_size = 0;
      const uint8_t* chunk_data =
          BrotliDecoderTakeOutput(ms.brotli, &chunk_size);
      ms.decompressed_size += chunk_size;
      if (ms.decompressed_size > ms.metadata_size) {
        return finish_decompression(BRUNSLI_INVALID_BRN);
      }
      size_t consumed_bytes = available_bytes - available_in;
      SkipBytes(state, consumed_bytes);
      bool chunk_ok = ProcessMetaData(chunk_data, chunk_size, &ms, jpg);
      if (!chunk_ok) return finish_decompression(BRUNSLI_INVALID_BRN);
      if (result == BROTLI_DECODER_RESULT_SUCCESS) {
        if (RemainingSectionLength(state) != 0) {
          return finish_decompression(BRUNSLI_INVALID_BRN);
        }
        if (ms.decompressed_size != ms.metadata_size) {
          return finish_decompression(BRUNSLI_INVALID_BRN);
        }
        if (!ms.CanFinish()) return finish_decompression(BRUNSLI_INVALID_BRN);
        return finish_decompression(BRUNSLI_OK);
      }
      if (result == BROTLI_DECODER_RESULT_NEEDS_MORE_OUTPUT) continue;
      BRUNSLI_DCHECK(result == BROTLI_DECODER_RESULT_NEEDS_MORE_INPUT);
      if (RemainingSectionLength(state) == 0) {
        return finish_decompression(BRUNSLI_INVALID_BRN);
      }
      return BRUNSLI_NOT_ENOUGH_DATA;
    }
  }

  // Unreachable.
  BRUNSLI_DCHECK(false);
  return BRUNSLI_DECOMPRESSION_ERROR;
}

static BrunsliStatus DecodeJPEGInternalsSection(State* state, JPEGData* jpg) {
  if (GetBytesAvailable(state) < RemainingSectionLength(state)) {
    return BRUNSLI_NOT_ENOUGH_DATA;
  }
  if (IsAtSectionBoundary(state)) return BRUNSLI_INVALID_BRN;

  // TODO(eustas): merge BitReader into State
  size_t section_len = RemainingSectionLength(state);
  BrunsliBitReader br;
  BrunsliBitReaderInit(&br);
  BrunsliBitReaderResume(&br, state->data + state->pos, section_len);
  if (!DecodeAuxData(&br, jpg)) return BRUNSLI_INVALID_BRN;
  size_t tail_length = BrunsliBitReaderSuspend(&br);
  BrunsliBitReaderFinish(&br);
  if (!BrunsliBitReaderIsHealthy(&br)) return BRUNSLI_INVALID_BRN;
  size_t consumed = section_len - tail_length;
  state->pos += consumed;

  for (size_t i = 0; i < jpg->marker_order.size(); ++i) {
    if (jpg->marker_order[i] != 0xff) {
      continue;
    }
    size_t data_size = 0;
    if (!DecodeDataLength(state, &data_size)) {
      return BRUNSLI_INVALID_BRN;
    }
    jpg->inter_marker_data.emplace_back(
        reinterpret_cast<const char*>(state->data + state->pos), data_size);
    state->pos += data_size;
  }
  return BRUNSLI_OK;
}

static BrunsliStatus DecodeQuantDataSection(State* state, JPEGData* jpg) {
  if (GetBytesAvailable(state) < RemainingSectionLength(state)) {
    return BRUNSLI_NOT_ENOUGH_DATA;
  }
  if (IsAtSectionBoundary(state)) return BRUNSLI_INVALID_BRN;

  // TODO(eustas): merge BitReader into State
  size_t section_len = RemainingSectionLength(state);
  BrunsliBitReader br;
  BrunsliBitReaderInit(&br);
  BrunsliBitReaderResume(&br, state->data + state->pos, section_len);
  if (!DecodeQuantTables(&br, jpg)) return BRUNSLI_INVALID_BRN;
  size_t tail_length = BrunsliBitReaderSuspend(&br);
  BrunsliBitReaderFinish(&br);
  if (!BrunsliBitReaderIsHealthy(&br)) return BRUNSLI_INVALID_BRN;
  if (tail_length != 0) return BRUNSLI_INVALID_BRN;
  state->pos += section_len;
  return BRUNSLI_OK;
}

static BrunsliStatus DecodeHistogramDataSection(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;

  if (GetBytesAvailable(state) < RemainingSectionLength(state)) {
    return BRUNSLI_NOT_ENOUGH_DATA;
  }
  if (IsAtSectionBoundary(state)) return BRUNSLI_INVALID_BRN;

  size_t num_components = jpg->components.size();
  BRUNSLI_DCHECK(num_components != 0);

  std::vector<ComponentMeta>& meta = state->meta;

  // TODO(eustas): merge BitReader into State
  size_t section_len = RemainingSectionLength(state);
  BrunsliBitReader br;
  BrunsliBitReaderInit(&br);
  BrunsliBitReaderResume(&br, state->data + state->pos, section_len);

  size_t num_contexts = num_components;
  for (size_t i = 0; i < num_components; ++i) {
    int scheme = BrunsliBitReaderRead(&br, 3);
    if (scheme >= kNumSchemes) return BRUNSLI_INVALID_BRN;
    meta[i].context_bits = scheme;
    meta[i].context_offset = num_contexts;
    num_contexts += kNumNonzeroContextSkip[scheme];
  }
  s.num_contexts = num_contexts;

  s.num_histograms = DecodeVarLenUint8(&br) + 1;
  if (!BrunsliBitReaderIsHealthy(&br)) return BRUNSLI_INVALID_BRN;

  if (!s.shallow_histograms) {
    s.context_map_.resize(s.num_contexts * kNumAvrgContexts);
    if (!DecodeContextMap(s.num_histograms, s.context_map_.size(),
                          s.context_map_.data(), &br)) {
      return BRUNSLI_INVALID_BRN;
    }
    state->context_map = s.context_map_.data();

    s.entropy_codes_.resize(s.num_histograms);
    for (size_t i = 0; i < s.num_histograms; ++i) {
      if (!s.entropy_codes_[i].ReadFromBitStream(kCoeffAlphabetSize, &br)) {
        return BRUNSLI_INVALID_BRN;
      }
    }
    state->entropy_codes = s.entropy_codes_.data();
  }
  // Finalize bit-reader state anyways.
  size_t tail_length = BrunsliBitReaderSuspend(&br);
  BrunsliBitReaderFinish(&br);
  if (!s.shallow_histograms) {
    if (tail_length != 0) return BRUNSLI_INVALID_BRN;
    if (!BrunsliBitReaderIsHealthy(&br)) return BRUNSLI_INVALID_BRN;
  }

  state->pos += section_len;
  return BRUNSLI_OK;
}

static BrunsliStatus DecodeDCDataSection(State* state) {
  if (GetBytesAvailable(state) < RemainingSectionLength(state)) {
    return BRUNSLI_NOT_ENOUGH_DATA;
  }

  size_t section_len = RemainingSectionLength(state);
  WordSource in(state->data + state->pos, section_len);

  if (!DecodeDC(state, &in)) return BRUNSLI_INVALID_BRN;

  if (in.len_ != in.pos_) return BRUNSLI_INVALID_BRN;
  state->pos += section_len;
  return BRUNSLI_OK;
}

static BrunsliStatus DecodeACDataSection(State* state) {
  if (GetBytesAvailable(state) < RemainingSectionLength(state)) {
    return BRUNSLI_NOT_ENOUGH_DATA;
  }

  size_t section_len = RemainingSectionLength(state);
  WordSource in(state->data + state->pos, section_len);

  if (!DecodeAC(state, &in)) return BRUNSLI_INVALID_BRN;

  if (in.len_ != in.pos_) return BRUNSLI_INVALID_BRN;
  state->pos += section_len;
  return BRUNSLI_OK;
}

static Stage DecodeOriginalJpg(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;
  FallbackState& fs = s.fallback;

  while (fs.stage != FallbackState::DONE) {
    switch (fs.stage) {
      case FallbackState::READ_TAG: {
        BrunsliStatus status = ReadTag(state, &s.section);
        if (status != BRUNSLI_OK) return Fail(state, status);
        if (s.section.tag != kBrunsliOriginalJpgTag || !s.section.is_section) {
          return Fail(state, BRUNSLI_INVALID_BRN);
        }
        fs.stage = FallbackState::ENTER_SECTION;
        break;
      }

      case FallbackState::ENTER_SECTION: {
        BrunsliStatus status = EnterSection(state, &s.section);
        if (status != BRUNSLI_OK) return Fail(state, status);
        jpg->original_jpg_size = s.section.remaining;
        // Edge case - empty payload.
        if (jpg->original_jpg_size == 0) {
          jpg->original_jpg = nullptr;
          fs.stage = FallbackState::DONE;
          break;
        }
        fs.stage = FallbackState::READ_CONTENTS;
        break;
      }

      case FallbackState::READ_CONTENTS: {
        size_t chunk_size = GetBytesAvailable(state);
        if (chunk_size == 0) {
          // TODO(eustas): dcheck s.section.remaining != 0
          return Fail(state, BRUNSLI_NOT_ENOUGH_DATA);
        }
        // Check if it is possible to avoid copy.
        const uint8_t* src = state->data + state->pos;
        if (fs.storage.empty()) {
          if (chunk_size >= jpg->original_jpg_size) {
            jpg->original_jpg = src;
            SkipBytes(state, jpg->original_jpg_size);
            fs.stage = FallbackState::DONE;
            break;
          }
        }
        // Otherwise, copy input.
        size_t remaining = jpg->original_jpg_size - fs.storage.size();
        size_t to_copy = std::min(chunk_size, remaining);
        fs.storage.insert(fs.storage.cend(), src, src + to_copy);
        SkipBytes(state, to_copy);
        if (fs.storage.size() == jpg->original_jpg_size) {
          jpg->original_jpg = fs.storage.data();
          fs.stage = FallbackState::DONE;
          break;
        }
        // TODO(eustas): dcheck GetBytesAvailable(state) == 0
        return Fail(state, BRUNSLI_NOT_ENOUGH_DATA);
      }

      default: return Fail(state, BRUNSLI_DECOMPRESSION_ERROR);
    }
  }

  LeaveSection(&s.section);
  return Stage::DONE;
}

static bool HasSection(State* state, uint32_t tag) {
  return state->internal->section.tags_met & (1u << tag);
}

static Stage ParseSection(State* state) {
  InternalState& s = *state->internal;
  SectionHeaderState& sh = s.section_header;

  Stage result = Stage::ERROR;

  while (sh.stage != SectionHeaderState::DONE) {
    switch (sh.stage) {
      case SectionHeaderState::READ_TAG: {
        BrunsliStatus status = ReadTag(state, &s.section);
        if (status == BRUNSLI_NOT_ENOUGH_DATA) {
          if (HasSection(state, kBrunsliACDataTag)) return Stage::DONE;
        }
        if (status != BRUNSLI_OK) return Fail(state, status);
        if (s.section.is_section) {
          sh.stage = SectionHeaderState::ENTER_SECTION;
          continue;
        }
        const uint32_t tag_bit = 1u << s.section.tag;
        const bool is_known_section_tag = kKnownSectionTags & tag_bit;
        if (is_known_section_tag) return Fail(state, BRUNSLI_INVALID_BRN);
        sh.stage = SectionHeaderState::READ_VALUE;
        continue;
      }

      case SectionHeaderState::READ_VALUE: {
        // No known varint tags on top level.
        size_t dummy;
        BrunsliStatus status = DecodeBase128(state, &dummy);
        if (status != BRUNSLI_OK) return Fail(state, status);
        result = Stage::SECTION;
        sh.stage = SectionHeaderState::DONE;
        continue;
      }

      case SectionHeaderState::ENTER_SECTION: {
        BrunsliStatus status = EnterSection(state, &s.section);
        if (status != BRUNSLI_OK) return Fail(state, status);
        result = Stage::SECTION_BODY;
        sh.stage = SectionHeaderState::DONE;
        continue;
      }

      default: return Fail(state, BRUNSLI_DECOMPRESSION_ERROR);
    }
  }

  sh.stage = SectionHeaderState::READ_TAG;
  BRUNSLI_DCHECK(result != Stage::ERROR);
  return result;
}

static Stage ProcessSection(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;

  const int32_t tag_bit = 1u << s.section.tag;
  const bool is_known_section_tag = kKnownSectionTags & tag_bit;

  const bool skip_section =
      !is_known_section_tag || (state->skip_tags & tag_bit);

  if (skip_section) {
    // Skip section content.
    size_t to_skip =
        std::min(GetBytesAvailable(state), RemainingSectionLength(state));
    state->pos += to_skip;
    if (RemainingSectionLength(state) != 0) {
      BRUNSLI_DCHECK(GetBytesAvailable(state) == 0);
      return Fail(state, BRUNSLI_NOT_ENOUGH_DATA);
    }
    return Stage::SECTION;
  }

  switch (s.section.tag) {
    case kBrunsliMetaDataTag: {
      BrunsliStatus status = DecodeMetaDataSection(state, jpg);
      if (status != BRUNSLI_OK) return Fail(state, status);
      break;
    }

    case kBrunsliJPEGInternalsTag: {
      BrunsliStatus status = DecodeJPEGInternalsSection(state, jpg);
      if (status != BRUNSLI_OK) return Fail(state, status);
      break;
    }

    case kBrunsliQuantDataTag: {
      if (!HasSection(state, kBrunsliJPEGInternalsTag)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      BrunsliStatus status = DecodeQuantDataSection(state, jpg);
      if (status != BRUNSLI_OK) return Fail(state, status);
      break;
    }

    case kBrunsliHistogramDataTag: {
      if (!HasSection(state, kBrunsliJPEGInternalsTag)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      BrunsliStatus status = DecodeHistogramDataSection(state, jpg);
      if (status != BRUNSLI_OK) return Fail(state, status);
      break;
    }

    case kBrunsliDCDataTag: {
      if (!HasSection(state, kBrunsliHistogramDataTag)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      if (!HasSection(state, kBrunsliQuantDataTag)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      internal::dec::WarmupMeta(jpg, state);
      BrunsliStatus status = DecodeDCDataSection(state);
      if (status != BRUNSLI_OK) return Fail(state, status);
      break;
    }

    case kBrunsliACDataTag: {
      if (!HasSection(state, kBrunsliDCDataTag)) {
        return Fail(state, BRUNSLI_INVALID_BRN);
      }
      internal::dec::WarmupMeta(jpg, state);
      BrunsliStatus status = DecodeACDataSection(state);
      if (status != BRUNSLI_OK) return Fail(state, status);
      break;
    }

    default:
      /* Unreachable */
      return Fail(state, BRUNSLI_INVALID_BRN);
  }

  if (!IsAtSectionBoundary(state)) {
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
    // 8205 == max[ceil((65535 / (i * 8)) * i) for i in range(1, 16 + 1)]
    BRUNSLI_DCHECK(c->width_in_blocks <= 8205);
    BRUNSLI_DCHECK(c->height_in_blocks <= 8205);
    uint32_t num_blocks = c->width_in_blocks * c->height_in_blocks;
    if (num_blocks > kBrunsliMaxNumBlocks) {
      return false;
    }
    c->num_blocks = num_blocks;
  }
  return true;
}

void PrepareMeta(const JPEGData* jpg, State* state) {
  InternalState& s = *state->internal;

  size_t num_components = jpg->components.size();
  s.block_state_.resize(num_components);
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
  InternalState& s = *state->internal;
  std::vector<ComponentMeta>& meta = state->meta;
  const size_t num_components = meta.size();

  if (!state->is_storage_allocated) {
    state->is_storage_allocated = true;
    for (size_t i = 0; i < num_components; ++i) {
      size_t num_blocks = meta[i].width_in_blocks * meta[i].height_in_blocks;
      jpg->components[i].coeffs.resize(num_blocks * kDCTBlockSize);
      s.block_state_[i].resize(num_blocks);
      meta[i].block_state = s.block_state_[i].data();
    }
  }

  if (!s.is_meta_warm) {
    s.is_meta_warm = true;
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

BrunsliStatus DoProcessJpeg(State* state, JPEGData* jpg) {
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
        if (state->pos != state->len) {
          state->stage = Fail(state, BRUNSLI_INVALID_BRN);
          break;
        }
        return BRUNSLI_OK;

      case Stage::ERROR:
        return state->internal->result;

      default:
        /* Unreachable */
        state->stage = Fail(state, BRUNSLI_DECOMPRESSION_ERROR);
        break;
    }
  }
}

BrunsliStatus ProcessJpeg(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;

  if (s.section.is_active) {
    s.section.milestone = state->pos;
    s.section.projected_end = s.section.milestone + s.section.remaining;
  }

  if (state->stage == Stage::ERROR) {
    // General error -> no recovery.
    if (s.result != BRUNSLI_NOT_ENOUGH_DATA) {
      return s.result;
    }
    // Continue parsing.
    s.result = BRUNSLI_OK;
    state->stage = s.last_stage;
    s.last_stage = Stage::ERROR;
  }

  s.section.tags_met |= state->tags_met;
  BrunsliStatus result = DoProcessJpeg(state, jpg);

  if (s.section.is_active) {
    // TODO(eustas): dcheck state->pos > s.section.milestone
    size_t processed_len = state->pos - s.section.milestone;
    // TODO(eustas): dcheck processed_len < s.section.remaining
    s.section.remaining -= processed_len;
  }

  return result;
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
  InternalState& s = *state.internal;
  s.shallow_histograms = true;

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
  size_t context_map_size = s.num_contexts * kNumAvrgContexts * sizeof(int32_t);
  size_t histogram_size = s.num_histograms * sizeof(ANSDecodingData);
  size_t decode_peak = context_map_size + histogram_size + component_state_size;
  size_t jpeg_writer_size = (1u << 17u) + (1u << 16u) * sizeof(int32_t);
  return (out_size + jpeg_data_size + std::max(decode_peak, jpeg_writer_size));
}

}  // namespace brunsli
