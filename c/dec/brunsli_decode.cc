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
#include "./histogram_decode.h"
#include "./huffman_table.h"
#include <brunsli/jpeg_data_writer.h>
#include "./state.h"
#include "./state_internal.h"

namespace brunsli {

using ::brunsli::internal::dec::AcDcState;
using ::brunsli::internal::dec::BlockI32;
using ::brunsli::internal::dec::Buffer;
using ::brunsli::internal::dec::ComponentMeta;
using ::brunsli::internal::dec::FallbackState;
using ::brunsli::internal::dec::HeaderState;
using ::brunsli::internal::dec::HistogramDataState;
using ::brunsli::internal::dec::InternalState;
using ::brunsli::internal::dec::JpegInternalsState;
using ::brunsli::internal::dec::MetadataDecompressionStage;
using ::brunsli::internal::dec::MetadataState;
using ::brunsli::internal::dec::PrepareMeta;
using ::brunsli::internal::dec::QuantDataState;
using ::brunsli::internal::dec::SectionHeaderState;
using ::brunsli::internal::dec::SectionState;
using ::brunsli::internal::dec::SerializationStatus;
using ::brunsli::internal::dec::Stage;
using ::brunsli::internal::dec::State;
using ::brunsli::internal::dec::UpdateSubsamplingDerivatives;
using ::brunsli::internal::dec::VarintState;

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

bool DecodeVarint(VarintState* s, BrunsliBitReader* br, size_t max_bits) {
  if (s->stage == VarintState::INIT) {
    s->value = 0;
    s->i = 0;
    s->stage = VarintState::READ_CONTINUATION;
  }

  while (true) {
    switch (s->stage) {
      case VarintState::READ_CONTINUATION: {
        if (s->i >= max_bits) {
          s->stage = VarintState::INIT;
          return true;
        }
        if (s->i + 1 != max_bits) {
          if (!BrunsliBitReaderCanRead(br, 1)) return false;
          if (!BrunsliBitReaderRead(br, 1)) {
            s->stage = VarintState::INIT;
            return true;
          }
        }
        s->stage = VarintState::READ_DATA;
        continue;
      }
      case VarintState::READ_DATA: {
        if (!BrunsliBitReaderCanRead(br, 1)) return false;
        size_t next_bit = BrunsliBitReaderRead(br, 1);
        s->value |= next_bit << s->i;
        ++s->i;
        s->stage = VarintState::READ_CONTINUATION;
        continue;
      }
      default: {
        BRUNSLI_CHECK(false);
        return false;
      }
    }
  }
}

template <size_t kChunkSize>
bool DecodeLimitedVarint(VarintState* s, BrunsliBitReader* br,
                         size_t max_symbols) {
  if (s->stage == VarintState::INIT) {
    s->value = 0;
    s->i = 0;
    s->stage = VarintState::READ_CONTINUATION;
  }
  while (true) {
    switch (s->stage) {
      case VarintState::READ_CONTINUATION: {
        if (s->i < max_symbols) {
          if (!BrunsliBitReaderCanRead(br, 1)) return false;
          if (BrunsliBitReaderRead(br, 1)) {
            s->stage = VarintState::READ_DATA;
            continue;
          }
        }
        s->stage = VarintState::INIT;
        return true;
      }
      case VarintState::READ_DATA: {
        if (!BrunsliBitReaderCanRead(br, kChunkSize)) return false;
        size_t next_bits = BrunsliBitReaderRead(br, kChunkSize);
        s->value |= next_bits << (s->i * kChunkSize);
        ++s->i;
        s->stage = VarintState::READ_CONTINUATION;
        continue;
      }
      default: {
        BRUNSLI_CHECK(false);
        return false;
      }
    }
  }
}

std::vector<uint8_t> GenerateApp0Marker(uint8_t app0_status) {
  std::vector<uint8_t> app0_marker(AppData_0xe0, AppData_0xe0 + 17);
  app0_marker[9] = app0_status & 1u ? 2 : 1;
  app0_status >>= 1u;
  app0_marker[10] = app0_status & 0x3u;
  app0_status >>= 2u;
  uint16_t x_dens = kApp0Densities[app0_status];
  app0_marker[11] = app0_marker[13] = x_dens >> 8u;
  app0_marker[12] = app0_marker[14] = x_dens & 0xFFu;
  return app0_marker;
}

std::vector<uint8_t> GenerateAppMarker(uint8_t marker, uint8_t code) {
  std::vector<uint8_t> s;
  if (marker == 0x80) {
    s = std::vector<uint8_t>(AppData_0xe2, AppData_0xe2 + 3161);
    s[84] = code;
  } else if (marker == 0x81) {
    s = std::vector<uint8_t>(AppData_0xec, AppData_0xec + 18);
    s[15] = code;
  } else {
    BRUNSLI_DCHECK(marker == 0x82);
    s = std::vector<uint8_t>(AppData_0xee, AppData_0xee + 15);
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
          jpg->tail_data = std::vector<uint8_t>();
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
        Append(&jpg->tail_data, data + pos, data + len);
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
        auto* dest = (state->marker == 0xFE) ? &jpg->com_data : &jpg->app_data;
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
        Append(state->multibyte_sink, data + pos, chunk_size);
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

static BrunsliStatus DecodeHuffmanCode(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;
  JpegInternalsState& js = s.internals;
  BrunsliBitReader* br = &js.br;

  while (true) {
    switch (js.stage) {
      case JpegInternalsState::READ_HUFFMAN_LAST: {
        if (!BrunsliBitReaderCanRead(br, 1)) return BRUNSLI_NOT_ENOUGH_DATA;
        js.is_known_last_huffman_code = BrunsliBitReaderRead(br, 1);
        jpg->huffman_code.emplace_back();
        js.stage = JpegInternalsState::READ_HUFFMAN_SIMPLE;
        continue;
      }
      case JpegInternalsState::READ_HUFFMAN_SIMPLE: {
        if (!BrunsliBitReaderCanRead(br, 5 + !js.is_known_last_huffman_code)) {
          return BRUNSLI_NOT_ENOUGH_DATA;
        }
        JPEGHuffmanCode* huff = &jpg->huffman_code.back();

        huff->slot_id = BrunsliBitReaderRead(br, 2);
        js.is_dc_table = (BrunsliBitReaderRead(br, 1) == 0);
        huff->slot_id += js.is_dc_table ? 0 : 0x10;
        huff->is_last =
            js.is_known_last_huffman_code || BrunsliBitReaderRead(br, 1);
        huff->counts[0] = 0;
        int found_match = BrunsliBitReaderRead(br, 1);
        if (found_match) {
          if (js.is_dc_table) {
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
          js.stage = JpegInternalsState::HUFFMAN_UPDATE;
        } else {
          // One less bit is used than requested, but it is guaranteed to be
          // consumed in complex Huffman code case.
          js.p.Init(js.is_dc_table
                        ? std::vector<uint8_t>(kDefaultDCValues,
                                               std::end(kDefaultDCValues))
                        : std::vector<uint8_t>(kDefaultACValues,
                                               std::end(kDefaultACValues)));
          js.stage = JpegInternalsState::READ_HUFFMAN_MAX_LEN;
        }
        continue;
      }
      case JpegInternalsState::READ_HUFFMAN_MAX_LEN: {
        if (!BrunsliBitReaderCanRead(br, 4)) return BRUNSLI_NOT_ENOUGH_DATA;
        js.max_len = BrunsliBitReaderRead(br, 4) + 1;
        js.total_count = 0;
        js.max_count =
            js.is_dc_table ? kJpegDCAlphabetSize : kJpegHuffmanAlphabetSize;
        js.space = (1u << kJpegHuffmanMaxBitLength) -
                   (1u << (kJpegHuffmanMaxBitLength - js.max_len));
        js.i = 1;
        js.stage = JpegInternalsState::READ_HUFFMAN_COUNT;
        continue;
      }
      case JpegInternalsState::READ_HUFFMAN_COUNT: {
        JPEGHuffmanCode* huff = &jpg->huffman_code.back();
        if (js.i <= js.max_len) {
          size_t shift = kJpegHuffmanMaxBitLength - js.i;
          size_t count_limit =
              std::min(js.max_count - js.total_count, js.space >> shift);
          if (count_limit > 0) {
            int nbits =
                Log2FloorNonZero(static_cast<uint32_t>(count_limit)) + 1;
            if (!BrunsliBitReaderCanRead(br, nbits)) {
              return BRUNSLI_NOT_ENOUGH_DATA;
            }
            size_t count = BrunsliBitReaderRead(br, nbits);
            if (count > count_limit) {
              return BRUNSLI_INVALID_BRN;
            }
            huff->counts[js.i] = static_cast<int>(count);
            js.total_count += count;
            js.space -= count * (static_cast<size_t>(1) << shift);
          }
          ++js.i;
          continue;
        }
        ++huff->counts[js.max_len];
        js.i = 0;
        js.stage = JpegInternalsState::READ_HUFFMAN_PERMUTATION;
        continue;
      }
      case JpegInternalsState::READ_HUFFMAN_PERMUTATION: {
        JPEGHuffmanCode* huff = &jpg->huffman_code.back();
        if (js.i < js.total_count) {
          const int nbits = js.p.num_bits();
          if (!DecodeLimitedVarint<2>(&js.varint, br, (nbits + 1) >> 1u)) {
            return BRUNSLI_NOT_ENOUGH_DATA;
          }
          uint8_t value;
          if (!js.p.Remove(js.varint.value, &value)) {
            return BRUNSLI_INVALID_BRN;
          }
          huff->values[js.i] = value;
          ++js.i;
          continue;
        }
        huff->values[js.total_count] = kJpegHuffmanAlphabetSize;
        js.stage = JpegInternalsState::HUFFMAN_UPDATE;
        continue;
      }
      case JpegInternalsState::HUFFMAN_UPDATE: {
        // This stage does not perform reading -> transient.
        if (jpg->huffman_code.back().is_last) {
          js.terminal_huffman_code_count++;
        }
        if (js.is_known_last_huffman_code) {
          js.p.Clear();
          return BRUNSLI_OK;
        }
        if (jpg->huffman_code.size() >= kMaxDHTMarkers) {
          // Too many Huffman codes for a valid bit-stream. Normally, a jpeg
          // file can have any arbitrary number of DHT, DQT, etc. But i prefer
          // we force a reasonable lower bound instead of open door to likely
          // forged BRN input.
          return BRUNSLI_INVALID_BRN;
        }
        js.stage = JpegInternalsState::READ_HUFFMAN_LAST;
        continue;
      }
      default:
        return BRUNSLI_DECOMPRESSION_ERROR;
    }
  }
  return BRUNSLI_OK;
}

BrunsliStatus DecodeScanInfo(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;
  JpegInternalsState& js = s.internals;
  BrunsliBitReader* br = &js.br;

  const auto maybe_add_zero_run = [&js, jpg] () {
    if (js.last_num > 0) {
      JPEGScanInfo::ExtraZeroRunInfo info;
      info.block_idx = js.last_block_idx;
      info.num_extra_zero_runs = js.last_num;
      jpg->scan_info[js.i].extra_zero_runs.push_back(info);
      js.last_num = 0;
    }
  };

  while (true) {
    switch (js.stage) {
      case JpegInternalsState::READ_SCAN_COMMON: {
        JPEGScanInfo* si = &jpg->scan_info[js.i];
        if (!BrunsliBitReaderCanRead(br, 22)) return BRUNSLI_NOT_ENOUGH_DATA;
        si->Ss = BrunsliBitReaderRead(br, 6);
        si->Se = BrunsliBitReaderRead(br, 6);
        si->Ah = BrunsliBitReaderRead(br, 4);
        si->Al = BrunsliBitReaderRead(br, 4);
        si->num_components = BrunsliBitReaderRead(br, 2) + 1;
        js.j = 0;
        js.stage = JpegInternalsState::READ_SCAN_COMPONENT;
        continue;
      }
      case JpegInternalsState::READ_SCAN_COMPONENT: {
        JPEGScanInfo* si = &jpg->scan_info[js.i];
        if (js.j < si->num_components) {
          if (!BrunsliBitReaderCanRead(br, 6)) return BRUNSLI_NOT_ENOUGH_DATA;
          si->components[js.j].comp_idx = BrunsliBitReaderRead(br, 2);
          si->components[js.j].dc_tbl_idx = BrunsliBitReaderRead(br, 2);
          si->components[js.j].ac_tbl_idx = BrunsliBitReaderRead(br, 2);
          js.j++;
        } else {
          js.last_block_idx = -1;
          js.stage = JpegInternalsState::READ_SCAN_RESET_POINT_CONTINUATION;
        }
        continue;
      }
      case JpegInternalsState::READ_SCAN_RESET_POINT_CONTINUATION: {
        if (!BrunsliBitReaderCanRead(br, 1)) return BRUNSLI_NOT_ENOUGH_DATA;
        if (BrunsliBitReaderRead(br, 1)) {
          js.stage = JpegInternalsState::READ_SCAN_RESET_POINT_DATA;
        } else {
          js.last_block_idx = 0;
          js.last_num = 0;
          js.stage = JpegInternalsState::READ_SCAN_ZERO_RUN_CONTINUATION;
        }
        continue;
      }
      case JpegInternalsState::READ_SCAN_RESET_POINT_DATA: {
        JPEGScanInfo* si = &jpg->scan_info[js.i];
        if (!DecodeVarint(&js.varint, br, 28)) return BRUNSLI_NOT_ENOUGH_DATA;
        int block_idx =
            js.last_block_idx + static_cast<int>(js.varint.value) + 1;
        si->reset_points.emplace_back(block_idx);
        js.last_block_idx = block_idx;
        // TODO(eustas): limit to exact number of blocks.
        if (js.last_block_idx > (1 << 30)) {
          // At most 8K x 8K x num_channels blocks are expected. That is,
          // typically, 1.5 * 2^27. 2^30 should be sufficient for any sane
          // image.
          return BRUNSLI_INVALID_BRN;
        }
        js.stage = JpegInternalsState::READ_SCAN_RESET_POINT_CONTINUATION;
        continue;
      }
      case JpegInternalsState::READ_SCAN_ZERO_RUN_CONTINUATION: {
        if (!BrunsliBitReaderCanRead(br, 1)) return BRUNSLI_NOT_ENOUGH_DATA;
        if (BrunsliBitReaderRead(br, 1)) {
          js.stage = JpegInternalsState::READ_SCAN_ZERO_RUN_DATA;
        } else {
          maybe_add_zero_run();
          ++js.i;
          if (js.i < js.num_scans) {
            js.stage = JpegInternalsState::READ_SCAN_COMMON;
            continue;
          }
          return BRUNSLI_OK;
        }
        continue;
      }
      case JpegInternalsState::READ_SCAN_ZERO_RUN_DATA: {
        if (!DecodeVarint(&js.varint, br, 28)) return BRUNSLI_NOT_ENOUGH_DATA;
        int block_idx = js.last_block_idx + static_cast<int>(js.varint.value);
        if (block_idx > js.last_block_idx) maybe_add_zero_run();
        ++js.last_num;
        js.last_block_idx = block_idx;
        // TODO(eustas): limit to exact number of blocks.
        if (js.last_block_idx > (1 << 30)) {
          // At most 8K x 8K x num_channels blocks are expected. That is,
          // typically, 1.5 * 2^27. 2^30 should be sufficient for any sane
          // image.
          return BRUNSLI_INVALID_BRN;
        }
        js.stage = JpegInternalsState::READ_SCAN_ZERO_RUN_CONTINUATION;
        continue;
      }
      default: return BRUNSLI_DECOMPRESSION_ERROR;
    }
  }
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

/** Reads 0..6 words from |in| and returns the value in the range 0..63. */
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

BrunsliStatus DecodeDC(State* state, WordSource* in) {
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

  if (!in->CanRead(5)) return BRUNSLI_NOT_ENOUGH_DATA;
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
      const int ac_stride = static_cast<int>(m.ac_stride);
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
          if (BRUNSLI_PREDICT_FALSE(!in->CanRead(6))) {
            ac_dc_state.next_mcu_y = mcu_y;
            ac_dc_state.next_component = i;
            ac_dc_state.next_iy = iy;
            ac_dc_state.next_x = x;
            s.ans_decoder = ans;
            s.bit_reader = br;
            s.arith_decoder = ac;
            return BRUNSLI_NOT_ENOUGH_DATA;
          }
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
                int nbits = code - kNumDirectCodes;
                Prob* BRUNSLI_RESTRICT p_first_extra_bit =
                    &c->first_extra_bit_prob[nbits];
                int first_extra_bit =
                    ac.ReadBit(p_first_extra_bit->get_proba(), in);
                p_first_extra_bit->Add(first_extra_bit);
                int extra_bits_val = first_extra_bit << nbits;
                if (nbits > 0) {
                  extra_bits_val |= static_cast<int>(br.ReadBits(nbits, in));
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
  if (!FinalizeSubdecoders(state)) return BRUNSLI_INVALID_BRN;

  return BRUNSLI_OK;
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
  const uint8_t* BRUNSLI_RESTRICT context_modes;
  const int* BRUNSLI_RESTRICT mult_col;
  const int* BRUNSLI_RESTRICT mult_row;
  int prev_row_delta;
  Prob* BRUNSLI_RESTRICT sign_prob;
  size_t context_bits;
  const uint8_t* BRUNSLI_RESTRICT context_map;
  const ANSDecodingData* BRUNSLI_RESTRICT entropy_codes;
  Prob* BRUNSLI_RESTRICT first_extra_bit_prob;
};

static size_t BRUNSLI_NOINLINE DecodeAcBlock(const AcBlockCookie& cookie) {
  AcBlockCookie c = cookie;

  BinaryArithmeticDecoder ac = *c.ac;
  WordSource* in = c.in;
  ANSDecoder ans = *c.ans;
  BitSource br = *c.br;

  size_t num_nonzeros = 0;

  const uint8_t nonzero_ctx = NumNonzerosContext(c.prev_num_nonzeros, c.x, c.y);
  size_t last_nz = DecodeNumNonzeros(
      c.num_nonzero_prob + kNumNonZeroTreeSize * nonzero_ctx, &ac, in);
  for (size_t k = last_nz + 1; k < kDCTBlockSize; ++k) {
    c.prev_sgn[k] = 0;
    c.prev_abs[k] = 0;
  }
  for (size_t k = last_nz; k > 0; --k) {
    int is_zero = 0;
    if (k < last_nz) {
      size_t bucket = kNonzeroBuckets[num_nonzeros - 1];
      size_t is_zero_ctx = bucket * kDCTBlockSize + k;
      Prob& p = c.is_zero_prob[is_zero_ctx];
      is_zero = ac.ReadBit(p.get_proba(), in);
      p.Add(is_zero);
    }
    int abs_val = 0;
    int sign = 1;
    const int k_nat = c.order[k];
    if (!is_zero) {
      size_t context_type = c.context_modes[k_nat];
      size_t avg_ctx = 0;
      size_t sign_ctx = kMaxAverageContext;
      if ((context_type & 1) && (c.y > 0)) {
        size_t offset = k_nat & 7;
        ACPredictContextRow(c.prev_row_coeffs + offset, c.coeffs + offset,
                            c.mult_col + offset * 8, &avg_ctx, &sign_ctx);
      } else if ((context_type & 2) && (c.x > 0)) {
        size_t offset = k_nat & ~7;
        ACPredictContextCol(c.prev_col_coeffs + offset, c.coeffs + offset,
                            c.mult_row + offset, &avg_ctx, &sign_ctx);
      } else if (!context_type) {
        avg_ctx = WeightedAverageContext(c.prev_abs + k, c.prev_row_delta);
        sign_ctx =
            c.prev_sgn[k] * 3 + c.prev_sgn[static_cast<int>(k) - kDCTBlockSize];
      }
      sign_ctx = sign_ctx * kDCTBlockSize + k;
      Prob& sign_p = c.sign_prob[sign_ctx];
      sign = ac.ReadBit(sign_p.get_proba(), in);
      sign_p.Add(sign);
      c.prev_sgn[k] = sign + 1;
      sign = 1 - 2 * sign;
      const size_t z_dens_ctx =
          ZeroDensityContext(num_nonzeros, k, c.context_bits);
      size_t histo_ix = z_dens_ctx * kNumAvrgContexts + avg_ctx;
      size_t entropy_ix = c.context_map[histo_ix];
      int code = ans.ReadSymbol(c.entropy_codes[entropy_ix], in);
      if (code < kNumDirectCodes) {
        abs_val = code + 1;
      } else {
        int nbits = code - kNumDirectCodes;
        Prob& p = c.first_extra_bit_prob[k * 10 + nbits];
        int first_extra_bit = ac.ReadBit(p.get_proba(), in);
        p.Add(first_extra_bit);
        int extra_bits_val = first_extra_bit << nbits;
        if (nbits > 0) {
          extra_bits_val |= br.ReadBits(nbits, in);
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

BrunsliStatus DecodeAC(State* state, WordSource* in) {
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

  if (!in->CanRead(5)) return BRUNSLI_NOT_ENOUGH_DATA;
  EnsureSubdecodersInitialized(state, in);

  if (!ac_dc_state.ac_coeffs_order_decoded) {
    while (ac_dc_state.next_component < num_components) {
      if (!in->CanRead(121)) return BRUNSLI_NOT_ENOUGH_DATA;
      if (!DecodeCoeffOrder(comps[ac_dc_state.next_component].order,
                            &s.bit_reader, in)) {
        return BRUNSLI_INVALID_BRN;
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
  c.context_modes =
      kContextAlgorithm + (state->use_legacy_context_model ? 64 : 0);

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
            if (BRUNSLI_PREDICT_FALSE(!in->CanRead(297))) {
              ac_dc_state.next_mcu_y = mcu_y;
              ac_dc_state.next_component = i;
              ac_dc_state.next_iy = iy;
              ac_dc_state.next_x = c.x;
              return BRUNSLI_NOT_ENOUGH_DATA;
            }
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
  ac_dc_state.next_mcu_y = 0;

  comps.clear();
  comps.shrink_to_fit();

  if (!FinalizeSubdecoders(state)) return BRUNSLI_INVALID_BRN;

  return BRUNSLI_OK;
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

        const size_t version = version_and_comp_count >> 2u;
        jpg->version = static_cast<int>(version);

        if (version == 1) {  // fallback mode
          // TODO(eustas): do we need this?
          jpg->width = 0;
          jpg->height = 0;
          hs.stage = HeaderState::DONE;
          break;
        }

        // Wrong mode = fallback + something.
        if ((version & 1u) != 0) {
          return Fail(state, BRUNSLI_INVALID_BRN);
        }
        // Unknown mode - only 3 bits are defined.
        if ((version & ~0x7u) != 0) {
          return Fail(state, BRUNSLI_INVALID_BRN);
        }

        // Otherwise regular brunsli.
        state->use_legacy_context_model = !(version & 2);

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
        jpg->width = static_cast<int>(width);
        jpg->height = static_cast<int>(height);

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
  return (jpg->version == 1) ? Stage::FALLBACK : Stage::SECTION;
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

/**
 * Wraps result, depending on the state of input.
 *
 * If parser needs more data, but section data is depleted,
 * then input is corrupted.
 */
static BrunsliStatus CheckBoundary(State* state, BrunsliStatus result) {
  if (result == BRUNSLI_NOT_ENOUGH_DATA) {
    bool last = (RemainingSectionLength(state) <= GetBytesAvailable(state));
    return last ? BRUNSLI_INVALID_BRN : BRUNSLI_NOT_ENOUGH_DATA;
  } else {
    return result;
  }
}

static void PrepareBitReader(BrunsliBitReader* br, State* state) {
  size_t chunk_len =
      std::min(GetBytesAvailable(state), RemainingSectionLength(state));
  BrunsliBitReaderResume(br, state->data + state->pos, chunk_len);
  BRUNSLI_DCHECK(BrunsliBitReaderIsHealthy(br));
}

/**
 * Marks data used by bit-reader as consumed.
 */
static BrunsliStatus SuspendBitReader(BrunsliBitReader* br, State* state,
                                      BrunsliStatus result) {
  size_t chunk_len =
      std::min(GetBytesAvailable(state), RemainingSectionLength(state));
  size_t unused_bytes = BrunsliBitReaderSuspend(br);
  size_t consumed_bytes = chunk_len - unused_bytes;
  SkipBytes(state, consumed_bytes);
  result = CheckBoundary(state, result);
  // Once BitReader becomes unhealthy, further decoding should be impossible.
  BRUNSLI_DCHECK(
      BrunsliBitReaderIsHealthy(br) ||
      ((result != BRUNSLI_OK) && (result != BRUNSLI_NOT_ENOUGH_DATA)));
  return result;
}

static BrunsliStatus DecodeJPEGInternalsSection(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;
  JpegInternalsState& js = s.internals;
  BrunsliBitReader* br = &js.br;

  if (js.stage == JpegInternalsState::INIT) {
    BrunsliBitReaderInit(br);
    js.stage = JpegInternalsState::READ_MARKERS;
  }
  PrepareBitReader(br, state);

  const auto suspend_bit_reader = [&](BrunsliStatus result) -> BrunsliStatus {
    return SuspendBitReader(br, state, result);
  };

  if (js.stage == JpegInternalsState::READ_MARKERS) {
    while (true) {
      if (!BrunsliBitReaderCanRead(br, 6)) {
        return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
      }
      uint8_t marker = 0xc0 + BrunsliBitReaderRead(br, 6);
      jpg->marker_order.push_back(marker);
      if (marker == 0xc4) ++js.dht_count;
      if (marker == 0xdd) js.have_dri = true;
      if (marker == 0xda) ++js.num_scans;
      if (marker == 0xd9) break;
    }
    js.stage = JpegInternalsState::READ_DRI;
  }

  if (js.stage == JpegInternalsState::READ_DRI) {
    if (js.have_dri) {
      if (!BrunsliBitReaderCanRead(br, 16)) {
        return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
      }
      jpg->restart_interval = BrunsliBitReaderRead(br, 16);
    }
    js.stage = JpegInternalsState::READ_HUFFMAN_LAST;
  }

  if (js.stage & JpegInternalsState::DECODE_HUFFMAN_MASK) {
    BrunsliStatus status = DecodeHuffmanCode(state, jpg);
    if (status != BRUNSLI_OK) return suspend_bit_reader(status);
    js.stage = JpegInternalsState::PREPARE_READ_SCANS;
  }

  if (js.stage == JpegInternalsState::PREPARE_READ_SCANS) {
    if (js.dht_count != js.terminal_huffman_code_count) {
      BRUNSLI_LOG_ERROR() << "Invalid number of DHT markers" << BRUNSLI_ENDL();
      return suspend_bit_reader(BRUNSLI_INVALID_BRN);
    }
    if (js.num_scans > 0) {
      jpg->scan_info.resize(js.num_scans);
      js.i = 0;
      js.stage = JpegInternalsState::READ_SCAN_COMMON;
    } else {
      js.stage = JpegInternalsState::READ_NUM_QUANT;
    }
  }

  if (js.stage & JpegInternalsState::DECODE_SCAN_MASK) {
    BrunsliStatus status = DecodeScanInfo(state, jpg);
    if (status != BRUNSLI_OK) return suspend_bit_reader(status);
    js.stage = JpegInternalsState::READ_NUM_QUANT;
  }

  if (js.stage == JpegInternalsState::READ_NUM_QUANT) {
    if (!BrunsliBitReaderCanRead(br, 2)) {
      return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
    }
    int num_quant_tables = BrunsliBitReaderRead(br, 2) + 1;
    jpg->quant.resize(num_quant_tables);
    js.i = 0;
    js.stage = JpegInternalsState::READ_QUANT;
  }

  while (js.stage == JpegInternalsState::READ_QUANT) {
    if (js.i >= jpg->quant.size()) {
      js.stage = JpegInternalsState::READ_COMP_ID_SCHEME;
      break;
    }
    // 6 or 7 bits are used, but we know that at least one more bit is
    // guaranteed to be used by varint out of the loop.
    if (!BrunsliBitReaderCanRead(br, 7)) {
      return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
    }
    JPEGQuantTable* q = &jpg->quant[js.i];
    q->index = BrunsliBitReaderRead(br, 2);
    q->is_last = (js.i == jpg->quant.size() - 1) || BrunsliBitReaderRead(br, 1);
    q->precision = BrunsliBitReaderRead(br, 4);
    if (q->precision > 1) {
      BRUNSLI_LOG_ERROR() << "Invalid quantization table precision: "
                          << q->precision << BRUNSLI_ENDL();
      return suspend_bit_reader(BRUNSLI_INVALID_BRN);
    }
    // note that q->values[] are initialized to invalid 0 values.
    ++js.i;
  }

  if (js.stage == JpegInternalsState::READ_COMP_ID_SCHEME) {
    if (!BrunsliBitReaderCanRead(br, 2)) {
      return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
    }
    int comp_ids = BrunsliBitReaderRead(br, 2);
    static const size_t kMinRequiredComponents[4] = {
        3 /* Ids123*/, 1 /* IdsGray */, 3 /* IdsRGB */, 0 /* IdsCustom */
    };
    if (jpg->components.size() < kMinRequiredComponents[comp_ids]) {
      BRUNSLI_LOG_ERROR() << "Insufficient number of components for ColorId #"
                          << comp_ids << BRUNSLI_ENDL();
      return suspend_bit_reader(BRUNSLI_INVALID_BRN);
    }
    js.stage = JpegInternalsState::READ_NUM_PADDING_BITS;
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
      js.i = 0;
      js.stage = JpegInternalsState::READ_COMP_ID;
    }
  }

  if (js.stage == JpegInternalsState::READ_COMP_ID) {
    while (js.i < jpg->components.size()) {
      if (!BrunsliBitReaderCanRead(br, 8)) {
        return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
      }
      jpg->components[js.i].id = BrunsliBitReaderRead(br, 8);
      ++js.i;
    }
    js.stage = JpegInternalsState::READ_NUM_PADDING_BITS;
  }

  if (js.stage == JpegInternalsState::READ_NUM_PADDING_BITS) {
    // TODO(eustas): sanitize: should not be bigger than
    //               7 x (num_scans + num_blocks / dri)
    // security: limit is 32b for n_size
    if (!DecodeLimitedVarint<8>(&js.varint, br, 4)) {
      return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
    }
    js.num_padding_bits = js.varint.value;
    jpg->has_zero_padding_bit = (js.num_padding_bits > 0);
    if (js.num_padding_bits > PaddingBitsLimit(*jpg)) {
      BRUNSLI_LOG_ERROR() << "Suspicious number of padding bits "
                          << js.num_padding_bits << BRUNSLI_ENDL();
      return suspend_bit_reader(BRUNSLI_INVALID_BRN);
    }
    js.i = 0;
    js.stage = JpegInternalsState::READ_PADDING_BITS;
  }

  if (js.stage == JpegInternalsState::READ_PADDING_BITS) {
    while (js.i < js.num_padding_bits) {
      if (!BrunsliBitReaderCanRead(br, 1)) {
        return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
      }
      jpg->padding_bits.emplace_back(BrunsliBitReaderRead(br, 1));
      ++js.i;
    }
    suspend_bit_reader(BRUNSLI_OK);
    BrunsliBitReaderFinish(br);
    if (!BrunsliBitReaderIsHealthy(br)) return BRUNSLI_INVALID_BRN;
    js.i = 0;
    js.stage = JpegInternalsState::ITERATE_MARKERS;
  } else {
    // no-op
    suspend_bit_reader(BRUNSLI_OK);
  }

  while (true) {
    switch (js.stage) {
      case JpegInternalsState::ITERATE_MARKERS: {
        if (js.i >= jpg->marker_order.size()) {
          js.stage = JpegInternalsState::DONE;
        } else if (jpg->marker_order[js.i] == 0xFF) {
          js.stage = JpegInternalsState::READ_INTERMARKER_LENGTH;
        } else {
          ++js.i;
        }
        continue;
      }

      case JpegInternalsState::READ_INTERMARKER_LENGTH: {
        BrunsliStatus status = DecodeBase128(state, &js.intermarker_length);
        if (status != BRUNSLI_OK) return CheckBoundary(state, status);
        if (js.intermarker_length > RemainingSectionLength(state)) {
          return BRUNSLI_INVALID_BRN;
        }
        jpg->inter_marker_data.emplace_back();
        js.stage = JpegInternalsState::READ_INTERMARKER_DATA;
        continue;
      }

      case JpegInternalsState::READ_INTERMARKER_DATA: {
        auto& dest = jpg->inter_marker_data.back();
        size_t piece_limit = js.intermarker_length - dest.size();
        size_t piece_size = std::min(piece_limit, GetBytesAvailable(state));
        Append(&dest, state->data + state->pos, piece_size);
        SkipBytes(state, piece_size);
        if (dest.size() < js.intermarker_length) {
          BRUNSLI_DCHECK(GetBytesAvailable(state) == 0);
          BRUNSLI_DCHECK(RemainingSectionLength(state) > 0);
          return BRUNSLI_NOT_ENOUGH_DATA;
        }
        ++js.i;
        js.stage = JpegInternalsState::ITERATE_MARKERS;
        continue;
      }

      default: { /* no-op */ }
    }
    break;  // no matching stage has been found; exit the loop.
  }

  if (!IsAtSectionBoundary(state)) return BRUNSLI_INVALID_BRN;

  return BRUNSLI_OK;
}

static BrunsliStatus DecodeQuantDataSection(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;
  QuantDataState& qs = s.quant;
  BrunsliBitReader* br = &qs.br;

  if (qs.stage == QuantDataState::INIT) {
    BrunsliBitReaderInit(br);
    qs.stage = QuantDataState::READ_NUM_QUANT;
  }
  PrepareBitReader(br, state);

  const auto suspend_bit_reader = [&](BrunsliStatus result) -> BrunsliStatus {
    return SuspendBitReader(br, state, result);
  };

  if (qs.stage == QuantDataState::READ_NUM_QUANT) {
    if (!BrunsliBitReaderCanRead(br, 2)) {
      return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
    }
    size_t num_quant_tables = BrunsliBitReaderRead(br, 2) + 1;
    if (jpg->quant.size() != num_quant_tables) {
      return suspend_bit_reader(BRUNSLI_INVALID_BRN);
    }
    qs.predictor.resize(kDCTBlockSize);
    qs.i = 0;
    qs.stage = QuantDataState::READ_STOCK;
  }

  while (true) {
    switch (qs.stage) {
      case QuantDataState::READ_STOCK: {
        if (qs.i >= jpg->quant.size()) {
          std::vector<uint8_t>().swap(qs.predictor);
          qs.i = 0;
          qs.stage = QuantDataState::READ_QUANT_IDX;
          continue;
        }
        // Depending on еру 1-st bit, it is guaranteed that we will need to read
        // at least 3 or 6 more bits.
        if (!BrunsliBitReaderCanRead(br, 4)) {
          return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
        }
        qs.data_precision = 0;
        bool is_short = !BrunsliBitReaderRead(br, 1);
        if (is_short) {
          const size_t short_code = BrunsliBitReaderRead(br, 3);
          int32_t* table = jpg->quant[qs.i].values.data();
          size_t selector = (qs.i > 0) ? 1 : 0;
          for (size_t k = 0; k < kDCTBlockSize; ++k) {
            table[k] = kStockQuantizationTables[selector][short_code][k];
          }
          qs.stage = QuantDataState::UPDATE;
        } else {
          qs.stage = QuantDataState::READ_Q_FACTOR;
        }
        continue;
      }

      case QuantDataState::READ_Q_FACTOR: {
        if (!BrunsliBitReaderCanRead(br, 6)) {
          return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
        }
        const uint32_t q_factor = BrunsliBitReaderRead(br, 6);
        FillQuantMatrix(qs.i > 0, q_factor, qs.predictor.data());
        qs.j = 0;
        qs.delta = 0;
        qs.stage = QuantDataState::READ_DIFF_IS_ZERO;
        continue;
      }

      case QuantDataState::READ_DIFF_IS_ZERO: {
        if (qs.j >= kDCTBlockSize) {
          qs.stage = QuantDataState::UPDATE;
          continue;
        }
        if (!BrunsliBitReaderCanRead(br, 1)) {
          return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
        }
        if (BrunsliBitReaderRead(br, 1)) {
          qs.stage = QuantDataState::READ_DIFF_SIGN;
        } else {
          qs.stage = QuantDataState::APPLY_DIFF;
        }
        continue;
      }

      case QuantDataState::READ_DIFF_SIGN: {
        if (!BrunsliBitReaderCanRead(br, 1)) {
          return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
        }
        qs.sign = BrunsliBitReaderRead(br, 1) ? -1 : 1;
        qs.stage = QuantDataState::READ_DIFF;
        continue;
      }

      case QuantDataState::READ_DIFF: {
        if (!DecodeVarint(&qs.vs, br, 16)) {
          return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
        }
        int diff = static_cast<int>(qs.vs.value) + 1;
        qs.delta += qs.sign * diff;
        qs.stage = QuantDataState::APPLY_DIFF;
        continue;
      }

      case QuantDataState::APPLY_DIFF: {
        const int k = kJPEGNaturalOrder[qs.j];
        const int quant_value = qs.predictor[k] + qs.delta;
        jpg->quant[qs.i].values[k] = quant_value;
        if (quant_value <= 0) {
          return suspend_bit_reader(BRUNSLI_INVALID_BRN);
        }
        if (quant_value >= 256) {
          qs.data_precision = 1;
        }
        if (quant_value >= 65536) {
          return suspend_bit_reader(BRUNSLI_INVALID_BRN);
        }
        ++qs.j;
        qs.stage = QuantDataState::READ_DIFF_IS_ZERO;
        continue;
      }

      case QuantDataState::UPDATE: {
        if (jpg->quant[qs.i].precision != qs.data_precision) {
          return suspend_bit_reader(BRUNSLI_INVALID_BRN);
        }
        ++qs.i;
        qs.stage = QuantDataState::READ_STOCK;
        continue;
      }

      default: { /* no-op */ }
    }
    break;  // no matching stage has been found; exit the loop.
  }

  while (qs.stage == QuantDataState::READ_QUANT_IDX) {
    if (qs.i >= jpg->components.size()) {
      qs.stage = QuantDataState::FINISH;
      continue;
    }
    JPEGComponent* c = &jpg->components[qs.i];
    if (!BrunsliBitReaderCanRead(br, 2)) {
       return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
    }
    c->quant_idx = BrunsliBitReaderRead(br, 2);
    if (c->quant_idx >= jpg->quant.size()) {
      return suspend_bit_reader(BRUNSLI_INVALID_BRN);
    }
    ++qs.i;
  }

  BRUNSLI_DCHECK(qs.stage == QuantDataState::FINISH);
  suspend_bit_reader(BRUNSLI_OK);
  BrunsliBitReaderFinish(br);
  if (!BrunsliBitReaderIsHealthy(br)) return BRUNSLI_INVALID_BRN;
  if (!IsAtSectionBoundary(state)) return BRUNSLI_INVALID_BRN;
  return BRUNSLI_OK;
}

static BrunsliStatus DecodeHistogramDataSection(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;
  HistogramDataState& hs = s.histogram;
  BrunsliBitReader* br = &hs.br;

  if (hs.stage == HistogramDataState::INIT) {
    BrunsliBitReaderInit(br);
    BRUNSLI_DCHECK(!jpg->components.empty());
    s.num_contexts = jpg->components.size();
    hs.stage = HistogramDataState::READ_SCHEME;
    /* Optimization: hint arena about maximal used alphabet size: 272.
     * 648 = 272 + 376, where 376 is a "universal" overhead for max-15-bit
     * root-8-bit 2-level Huffman tables. */
    hs.arena.reserve(648);
  }
  PrepareBitReader(br, state);
  if (RemainingSectionLength(state) <= GetBytesAvailable(state)) {
    // If end of section is reachable, then we could parse the remainings in
    // non-streaming mode.
    BrunsliBitReaderSetOptimistic(br);
  }
  const auto suspend_bit_reader = [&](BrunsliStatus result) -> BrunsliStatus {
    return SuspendBitReader(br, state, result);
  };

  if (hs.stage == HistogramDataState::READ_SCHEME) {
    const size_t num_components = jpg->components.size();
    BRUNSLI_DCHECK(num_components <= 4);
    if (!BrunsliBitReaderCanRead(br, 3 * num_components)) {
      return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
    }
    for (size_t i = 0; i < num_components; ++i) {
      size_t scheme = BrunsliBitReaderRead(br, 3);
      if (scheme >= kNumSchemes) return suspend_bit_reader(BRUNSLI_INVALID_BRN);
      ComponentMeta& m = state->meta[i];
      m.context_bits = scheme;
      m.context_offset = s.num_contexts;
      s.num_contexts += kNumNonzeroContextSkip[scheme];
    }
    if (!BrunsliBitReaderIsHealthy(br)) {
      return suspend_bit_reader(BRUNSLI_INVALID_BRN);
    }
    hs.stage = HistogramDataState::READ_NUM_HISTOGRAMS;
  }

  if (hs.stage == HistogramDataState::READ_NUM_HISTOGRAMS) {
    if (!BrunsliBitReaderCanRead(br, 11)) {
      return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
    }
    s.num_histograms = DecodeVarLenUint8(br) + 1;
    if (!BrunsliBitReaderIsHealthy(br)) {
      return suspend_bit_reader(BRUNSLI_INVALID_BRN);
    }
    if (s.shallow_histograms) {
      hs.stage = HistogramDataState::SKIP_CONTENT;
    } else {
      s.context_map_.resize(s.num_contexts * kNumAvrgContexts);
      state->context_map = s.context_map_.data();
      s.entropy_codes_.resize(s.num_histograms);
      state->entropy_codes = s.entropy_codes_.data();
      if (s.num_histograms > 1) {
        hs.stage = HistogramDataState::READ_CONTEXT_MAP_CODE;
      } else {
        hs.i = 0;
        hs.counts.resize(kCoeffAlphabetSize);
        hs.stage = HistogramDataState::READ_HISTOGRAMS;
      }
    }
  }

  if (hs.stage == HistogramDataState::SKIP_CONTENT) {
    suspend_bit_reader(BRUNSLI_OK);
    if (!BrunsliBitReaderIsHealthy(br)) {
      return suspend_bit_reader(BRUNSLI_INVALID_BRN);
    }
    SkipAvailableBytes(state, RemainingSectionLength(state));
    if (!IsAtSectionBoundary(state)) return BRUNSLI_NOT_ENOUGH_DATA;
    hs.stage = HistogramDataState::DONE;
  }

  if (hs.stage == HistogramDataState::READ_CONTEXT_MAP_CODE) {
    if (!BrunsliBitReaderCanRead(br, 207 + s.num_histograms * 8)) {
      return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
    }
    hs.max_run_length_prefix = 0;
    bool use_rle_for_zeros = !!BrunsliBitReaderRead(br, 1);
    if (use_rle_for_zeros) {
      hs.max_run_length_prefix = BrunsliBitReaderRead(br, 4) + 1;
    }
    size_t alphabet_size = s.num_histograms + hs.max_run_length_prefix;
    hs.entropy.reset(new HuffmanDecodingData);
    if (!hs.entropy->ReadFromBitStream(alphabet_size, br, &hs.arena)) {
      return suspend_bit_reader(BRUNSLI_INVALID_BRN);
    }
    hs.i = 0;
    hs.stage = HistogramDataState::READ_CONTEXT_MAP;
  }

  if (hs.stage == HistogramDataState::READ_CONTEXT_MAP) {
    BrunsliStatus status = DecodeContextMap(
        *hs.entropy, hs.max_run_length_prefix, &hs.i, &s.context_map_, br);
    if (status != BRUNSLI_OK) return suspend_bit_reader(status);
    hs.i = 0;
    hs.counts.resize(kCoeffAlphabetSize);
    hs.stage = HistogramDataState::READ_HISTOGRAMS;
  }

  if (hs.stage == HistogramDataState::READ_HISTOGRAMS) {
    while (hs.i < s.num_histograms) {
      if (!BrunsliBitReaderCanRead(br, 9 + kCoeffAlphabetSize * 11)) {
        return suspend_bit_reader(BRUNSLI_NOT_ENOUGH_DATA);
      }
      if (!ReadHistogram(BRUNSLI_ANS_LOG_TAB_SIZE, &hs.counts, br)) {
        return suspend_bit_reader(BRUNSLI_INVALID_BRN);
      }
      if (!s.entropy_codes_[hs.i].Init(hs.counts)) {
        return suspend_bit_reader(BRUNSLI_INVALID_BRN);
      }
      ++hs.i;
    }
    hs.entropy.reset();
    std::vector<uint32_t>().swap(hs.counts);
    suspend_bit_reader(BRUNSLI_OK);
    BrunsliBitReaderFinish(br);
    if (!BrunsliBitReaderIsHealthy(br)) return BRUNSLI_INVALID_BRN;
    if (!IsAtSectionBoundary(state)) return BRUNSLI_INVALID_BRN;
    hs.stage = HistogramDataState::DONE;
  }

  hs.arena.reset();
  BRUNSLI_DCHECK(hs.stage == HistogramDataState::DONE);
  return BRUNSLI_OK;
}

static BrunsliStatus DecodeDCDataSection(State* state) {
  size_t available = GetBytesAvailable(state) & ~1;
  size_t limit = RemainingSectionLength(state);
  BRUNSLI_DCHECK((limit & 1) == 0);
  size_t chunk_len = std::min(available, limit);
  // If end of section is reachable, then we could parse the remainings in
  // non-streaming mode.
  bool is_last_chunk = (chunk_len == limit);
  WordSource in(state->data + state->pos, chunk_len, is_last_chunk);

  BrunsliStatus status = DecodeDC(state, &in);

  BRUNSLI_DCHECK((in.pos_ & 1) == 0);
  if (in.error_) return BRUNSLI_INVALID_BRN;
  BRUNSLI_DCHECK(in.pos_ <= chunk_len);
  SkipBytes(state, in.pos_);
  if (is_last_chunk) {
    BRUNSLI_DCHECK(status != BRUNSLI_NOT_ENOUGH_DATA);
    if (!IsAtSectionBoundary(state)) return BRUNSLI_INVALID_BRN;
  }
  return status;
}

static BrunsliStatus DecodeACDataSection(State* state) {
  size_t available = GetBytesAvailable(state) & ~1;
  size_t limit = RemainingSectionLength(state);
  BRUNSLI_DCHECK((limit & 1) == 0);
  size_t chunk_len = std::min(available, limit);
  // If end of section is reachable, then we could parse the remainings in
  // non-streaming mode.
  bool is_last_chunk = (chunk_len == limit);
  WordSource in(state->data + state->pos, chunk_len, is_last_chunk);

  BrunsliStatus status = DecodeAC(state, &in);

  BRUNSLI_DCHECK((in.pos_ & 1) == 0);
  if (in.error_) return BRUNSLI_INVALID_BRN;
  BRUNSLI_DCHECK(in.pos_ <= chunk_len);
  SkipBytes(state, in.pos_);
  if (is_last_chunk) {
    BRUNSLI_DCHECK(status != BRUNSLI_NOT_ENOUGH_DATA);
    if (!IsAtSectionBoundary(state)) return BRUNSLI_INVALID_BRN;
  }
  return status;
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
      // This section reads input word by word.
      if ((RemainingSectionLength(state) & 1) != 0) {
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
      // This section reads input word by word.
      if ((RemainingSectionLength(state) & 1) != 0) {
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

  // Nothing is expected after the AC data.
  if (s.section.tag == kBrunsliACDataTag) {
    return Stage::DONE;
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

/** Adds new input to buffer. */
void ChargeBuffer(State* state) {
  InternalState& s = *state->internal;
  Buffer& b = s.buffer;

  b.borrowed_len = 0;
  b.external_data = state->data;
  b.external_pos = state->pos;
  b.external_len = state->len;
}

constexpr size_t kBufferMaxReadAhead = 600;

/** Sets input source either to buffered, or to external data. */
void LoadInput(State* state) {
  InternalState& s = *state->internal;
  Buffer& b = s.buffer;

  // No data buffered. Just pass external data as is.
  if (b.data_len == 0) {
    state->data = b.external_data;
    state->pos = b.external_pos;
    state->len = b.external_len;
    return;
  }

  BRUNSLI_DCHECK(b.data_len <= kBufferMaxReadAhead);

  // Otherwise use buffered data.
  size_t available = b.external_len - b.external_pos;
  // Always try to borrow as much as parser could require. This way, when
  // buffer is unable to provide enough input, we could switch to unbuffered
  // input.
  b.borrowed_len = std::min(kBufferMaxReadAhead, available);
  memcpy(b.data.data() + b.data_len, b.external_data + b.external_pos,
         b.borrowed_len);
  state->data = b.data.data();
  state->pos = 0;
  state->len = b.data_len + b.borrowed_len;
}

/**
 * Cancel borrowed bytes, if any.
 *
 * Returns false, if it is impossible to continue parsing.
 */
bool UnloadInput(State* state, BrunsliStatus result) {
  InternalState& s = *state->internal;
  Buffer& b = s.buffer;

  // Non-buffered input; put tail to buffer.
  if (state->data == b.external_data) {
    b.external_pos = state->pos;
    BRUNSLI_DCHECK(b.external_pos <= b.external_len);
    if (result != BRUNSLI_NOT_ENOUGH_DATA) return true;
    BRUNSLI_DCHECK(b.data_len == 0);
    size_t available = b.external_len - b.external_pos;
    BRUNSLI_DCHECK(available < kBufferMaxReadAhead);
    if (b.data.empty()) b.data.resize(2 * kBufferMaxReadAhead);
    b.data_len = available;
    memcpy(b.data.data(), b.external_data + b.external_pos, b.data_len);
    b.external_pos += available;
    return false;
  }

  // Buffer depleted; switch to non-buffered input.
  if (state->pos >= b.data_len) {
    size_t used_borrowed_bytes = state->pos - b.data_len;
    b.data_len = 0;
    b.external_pos += used_borrowed_bytes;
    return true;
  }

  // Buffer not depleted; either problem discovered was already buffered data,
  // or extra input was too-short.
  b.data_len -= state->pos;
  if (result == BRUNSLI_NOT_ENOUGH_DATA) {
    // We couldn't have taken more bytes.
    BRUNSLI_DCHECK(b.external_pos + b.borrowed_len == b.external_len);
    // Remaining piece is not too large.
    BRUNSLI_DCHECK(b.data_len + b.borrowed_len < kBufferMaxReadAhead);
    b.data_len += b.borrowed_len;
    b.external_pos += b.borrowed_len;
  }
  BRUNSLI_DCHECK(!b.data.empty());
  if (state->pos > 0 && b.data_len > 0) {
    memmove(b.data.data(), b.data.data() + state->pos, b.data_len);
  }
  BRUNSLI_DCHECK(b.data_len <= kBufferMaxReadAhead);

  return (result != BRUNSLI_NOT_ENOUGH_DATA);
}

/** Sets back user-provided input. */
void UnchargeBuffer(State* state) {
  InternalState& s = *state->internal;
  Buffer& b = s.buffer;

  state->data = b.external_data;
  state->pos = b.external_pos;
  state->len = b.external_len;
}

BrunsliStatus ProcessJpeg(State* state, JPEGData* jpg) {
  InternalState& s = *state->internal;

  if (state->pos > state->len) return BRUNSLI_INVALID_PARAM;
  ChargeBuffer(state);

  BrunsliStatus result = BRUNSLI_NOT_ENOUGH_DATA;
  while (result == BRUNSLI_NOT_ENOUGH_DATA) {
    if (state->stage == Stage::ERROR) {
      // General error -> no recovery.
      if (s.result != BRUNSLI_NOT_ENOUGH_DATA) return s.result;
      // Continue parsing.
      s.result = BRUNSLI_OK;
      state->stage = s.last_stage;
      s.last_stage = Stage::ERROR;
    }

    LoadInput(state);
    if (s.section.is_active) {
      s.section.milestone = state->pos;
      s.section.projected_end = s.section.milestone + s.section.remaining;
    }

    s.section.tags_met |= state->tags_met;
    result = DoProcessJpeg(state, jpg);

    if (s.section.is_active) {
      // TODO(eustas): dcheck state->pos > s.section.milestone
      size_t processed_len = state->pos - s.section.milestone;
      // TODO(eustas): dcheck processed_len < s.section.remaining
      s.section.remaining -= processed_len;
    }

    if (!UnloadInput(state, result)) break;
  }
  UnchargeBuffer(state);
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

BrunsliDecoder::BrunsliDecoder() {
  jpg_.reset(new JPEGData);
  state_.reset(new State);
}

BrunsliDecoder::~BrunsliDecoder() {}

BrunsliDecoder::Status BrunsliDecoder::Decode(size_t* available_in,
                                              const uint8_t** next_in,
                                              size_t* available_out,
                                              uint8_t** next_out) {
  JPEGData* jpg = jpg_.get();
  BRUNSLI_DCHECK(jpg);
  State* state = state_.get();
  BRUNSLI_DCHECK(state);

  state->data = *next_in;
  state->pos = 0;
  state->len = *available_in;
  BrunsliStatus parse_status = internal::dec::ProcessJpeg(state, jpg);
  size_t consumed_bytes = state->pos;
  *available_in -= consumed_bytes;
  *next_in += consumed_bytes;

  if ((parse_status != BRUNSLI_OK) &&
      (parse_status != BRUNSLI_NOT_ENOUGH_DATA)) {
    return BrunsliDecoder::ERROR;
  }

  // All the input given input should be consumed.
  BRUNSLI_DCHECK(*available_in == 0);

  SerializationStatus serialization_status =
      SerializeJpeg(state, *jpg, available_out, next_out);
  if (serialization_status == SerializationStatus::ERROR) {
    return BrunsliDecoder::ERROR;
  }

  switch (serialization_status) {
    case SerializationStatus::DONE:
      // Should be impossible to finish serialization without finishing parsing.
      BRUNSLI_DCHECK(parse_status == BRUNSLI_OK);
      return BrunsliDecoder::DONE;

    case SerializationStatus::NEEDS_MORE_INPUT:
      // If serializer says that data is incomplete, parser should say the same.
      BRUNSLI_DCHECK(parse_status == BRUNSLI_NOT_ENOUGH_DATA);
      return BrunsliDecoder::NEEDS_MORE_INPUT;

    case SerializationStatus::NEEDS_MORE_OUTPUT:
      // TODO(eustas): make sure that serializer could produce more bytes
      //               without providing more bytes to parser
      BRUNSLI_DCHECK(*available_out == 0);
      return BrunsliDecoder::NEEDS_MORE_OUTPUT;

    case SerializationStatus::ERROR:
      return BrunsliDecoder::ERROR;

    default:
      /* Unreachable */
      BRUNSLI_DCHECK(false);
      return BrunsliDecoder::ERROR;
  }
}

}  // namespace brunsli
