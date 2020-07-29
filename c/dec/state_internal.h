// Copyright (c) Google LLC 2020
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_STATE_INTERNAL_H_
#define BRUNSLI_DEC_STATE_INTERNAL_H_

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "../common/context.h"
#include "../common/lehmer_code.h"
#include <brunsli/status.h>
#include <brunsli/types.h>
#include "./ans_decode.h"
#include "./arith_decode.h"
#include "./bit_reader.h"
#include "./brunsli_input.h"
#include "./huffman_decode.h"
#include <brunsli/jpeg_data_writer.h>
#include "./serialization_state.h"
#include "./state.h"

struct BrotliDecoderStateStruct;

namespace brunsli {

struct HuffmanDecodingData;

namespace internal {
namespace dec {

struct AcDcState {
  int next_mcu_y = 0;
  size_t next_component = 0;
  int next_iy = 0;
  int next_x = 0;
  bool ac_coeffs_order_decoded = false;

  std::vector<ComponentState> ac;
  std::vector<ComponentStateDC> dc;
};

// Aid for section / subsection parsing.
struct SectionState {
  // Current value tag and type.
  size_t tag = 0;
  // True, if section is entered.
  bool is_active = false;
  // True, if "message" is actually "section", not a primitive value.
  bool is_section = false;

  // Encountered tags tracker.
  uint32_t tags_met = 0;

  // Remaining section length. Actual only when outside of workflow.
  size_t remaining = 0;

  // Position in current input, for which |remaining| was actual.
  size_t milestone = 0;

  // Projected section end, given enough input is provided.
  // |projected_end| == |milestone| + |provided|
  size_t projected_end = 0;
};

// Fields used for "Header" section parsing.
struct HeaderState {
  enum Stage {
    // Check section tag.
    READ_TAG,
    // Read section length.
    ENTER_SECTION,
    // Read value marker.
    ITEM_READ_TAG,
    // Read subsection length.
    ITEM_ENTER_SECTION,
    // Skip subsection payload.
    ITEM_SKIP_CONTENTS,
    // Read value.
    ITEM_READ_VALUE,
    // Verify values and apply to decoder state
    FINALE,
    // Finish section decoding.
    DONE
  };

  size_t stage = READ_TAG;

  // Subsection properties.
  SectionState section;
  // Length of subsection remaining to skip.
  size_t remaining_skip_length = 0;

  // Collected data (values).
  std::array<size_t, 16> varint_values;
};

// Fields used for "Fallback" section parsing.
struct FallbackState {
  enum Stage {
    // Check section tag.
    READ_TAG,
    // Read section length.
    ENTER_SECTION,
    // Copy "original JPEG" contents to internal storage (if necessary).
    READ_CONTENTS,
    // Finish section decoding.
    DONE
  };

  size_t stage = READ_TAG;

  // Storage for original JPEG contents.
  std::vector<uint8_t> storage;
};

// Fields used for section header parsing.
struct SectionHeaderState {
  enum Stage {
    // Check section tag.
    READ_TAG,
    // Read (dummy) value.
    READ_VALUE,
    // Read section length.
    ENTER_SECTION,
    // Finish section header decoding.
    DONE
  };

  size_t stage = READ_TAG;
};

enum class MetadataDecompressionStage {
  // Initial state in which it is decided which one of 3 processing variants to
  // use.
  INITIAL,
  // Read the length of uncompressed payload.
  READ_LENGTH,
  // Continuing as stream-decompressing/-parsing of Brotli-compressed metadata.
  DECOMPRESSING,
  // Parsing is finished, no further processing expected.
  DONE,
};

struct MetadataState {
  enum Stage {
    // Parse sequence type.
    READ_MARKER,
    // Dump the remaining of metadata to tail sequence.
    READ_TAIL,
    // Parse second byte of 2-byte sequence.
    READ_CODE,
    // Parse multi-byte sequence length.
    READ_LENGTH_HI,
    READ_LENGTH_LO,
    // Parse multi-byte sequence.
    READ_MULTIBYTE,
  };

  size_t short_marker_count = 0;
  uint8_t marker;
  uint8_t length_hi;
  size_t remaining_multibyte_length;
  std::vector<uint8_t>* multibyte_sink;
  size_t stage = READ_MARKER;

  BrotliDecoderStateStruct* brotli = nullptr;
  size_t metadata_size;
  size_t decompressed_size = 0;
  BrunsliStatus result = BRUNSLI_DECOMPRESSION_ERROR;
  MetadataDecompressionStage decompression_stage =
      MetadataDecompressionStage::INITIAL;

  ~MetadataState();

  bool CanFinish() { return (stage == READ_MARKER) || (stage == READ_TAIL); }
};

/**
 * Fits both DecodeVarint and DecodeLimitedVarint workflows.
 *
 * TODO(eustas): we could turn those methods back to stateless,
 * when "mark / rewind" utilities are added to BrunsliBitReader, and outer
 * parsing workflow supports input buffering.
 */
struct VarintState {
  enum Stage {
    INIT,
    READ_CONTINUATION,
    READ_DATA
  };

  Stage stage = INIT;
  size_t value;
  size_t i;
};

struct JpegInternalsState {
  enum Stage {
    INIT = 0,
    READ_MARKERS,
    READ_DRI,

    DECODE_HUFFMAN_MASK = 0x10,
    READ_HUFFMAN_LAST,
    READ_HUFFMAN_SIMPLE,
    READ_HUFFMAN_MAX_LEN,
    READ_HUFFMAN_COUNT,
    READ_HUFFMAN_PERMUTATION,
    HUFFMAN_UPDATE,

    PREPARE_READ_SCANS = 0x20,

    DECODE_SCAN_MASK = 0x40,
    READ_SCAN_COMMON,
    READ_SCAN_COMPONENT,
    READ_SCAN_RESET_POINT_CONTINUATION,
    READ_SCAN_RESET_POINT_DATA,
    READ_SCAN_ZERO_RUN_CONTINUATION,
    READ_SCAN_ZERO_RUN_DATA,

    READ_NUM_QUANT = 0x80,
    READ_QUANT,
    READ_COMP_ID_SCHEME,
    READ_COMP_ID,
    READ_NUM_PADDING_BITS,
    READ_PADDING_BITS,

    ITERATE_MARKERS,
    READ_INTERMARKER_LENGTH,
    READ_INTERMARKER_DATA,

    DONE
  };

  Stage stage = INIT;

  bool have_dri = false;
  size_t num_scans = 0;
  size_t dht_count = 0;

  BrunsliBitReader br;
  size_t is_known_last_huffman_code;
  size_t terminal_huffman_code_count = 0;
  bool is_dc_table;
  size_t total_count;
  size_t space;
  size_t max_len;
  size_t max_count;
  size_t i;
  PermutationCoder p;
  VarintState varint;

  size_t j;
  int last_block_idx;
  int last_num;

  size_t num_padding_bits;
  size_t intermarker_length;
};

struct QuantDataState {
  enum Stage {
    INIT,

    READ_NUM_QUANT,

    READ_STOCK,
    READ_Q_FACTOR,
    READ_DIFF_IS_ZERO,
    READ_DIFF_SIGN,
    READ_DIFF,
    APPLY_DIFF,
    UPDATE,

    READ_QUANT_IDX,

    FINISH
  };

  Stage stage = INIT;

  BrunsliBitReader br;
  size_t i;
  size_t j;
  uint8_t data_precision;
  VarintState vs;
  int delta;
  int sign;
  std::vector<uint8_t> predictor;
};

struct HistogramDataState {
  enum Stage {
    INIT,

    READ_SCHEME,
    READ_NUM_HISTOGRAMS,
    READ_CONTEXT_MAP_CODE,
    READ_CONTEXT_MAP,
    READ_HISTOGRAMS,

    SKIP_CONTENT,

    DONE
  };

  Stage stage = INIT;

  BrunsliBitReader br;
  size_t max_run_length_prefix;
  std::unique_ptr<HuffmanDecodingData> entropy;
  size_t i;
  std::vector<uint32_t> counts;
  Arena<HuffmanCode> arena;
};

struct Buffer {
  size_t data_len = 0;
  size_t borrowed_len;
  std::vector<uint8_t> data;

  const uint8_t* external_data;
  size_t external_pos;
  size_t external_len;
};

struct InternalState {
  /* Parsing */

  AcDcState ac_dc;
  SectionState section;

  // Sections.
  HeaderState header;
  FallbackState fallback;
  SectionHeaderState section_header;
  MetadataState metadata;
  JpegInternalsState internals;
  QuantDataState quant;
  HistogramDataState histogram;

  // "JPEGDecodingState" storage.
  std::vector<uint8_t> context_map_;
  std::vector<ANSDecodingData> entropy_codes_;
  std::vector<std::vector<uint8_t>> block_state_;

  bool is_meta_warm = false;

  // For "estimate peak memory".
  bool shallow_histograms = false;
  size_t num_contexts = 0;
  size_t num_histograms = 0;

  // Sub-decoders.
  bool subdecoders_initialized = false;
  ANSDecoder ans_decoder;
  BitSource bit_reader;
  BinaryArithmeticDecoder arith_decoder;

  BrunsliStatus result = BRUNSLI_OK;

  Stage last_stage = Stage::ERROR;

  Buffer buffer;

  /* Serialization */

  SerializationState serialization;
};

}  // namespace dec
}  // namespace internal
}  // namespace brunsli

#endif  // BRUNSLI_DEC_STATE_INTERNAL_H_
