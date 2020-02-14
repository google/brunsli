// Copyright (c) Google LLC 2020
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_STATE_INTERNAL_H_
#define BRUNSLI_DEC_STATE_INTERNAL_H_

#include <array>
#include <memory>
#include <vector>

#include "../common/context.h"
#include <brunsli/types.h>
#include <brunsli/status.h>
#include "./ans_decode.h"
#include "./arith_decode.h"
#include "./brunsli_input.h"
#include "./state.h"

namespace brunsli {
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

struct SectionState {
  // Current value tag and type.
  size_t tag = 0;
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

struct HeaderState {
  enum Stage {
    READ_TAG,
    ENTER_SECTION,
    ITEM_READ_TAG,
    ITEM_ENTER_SECTION,
    ITEM_SKIP_CONTENTS,
    ITEM_READ_VALUE,
    FINALE,
    DONE
  };

  // Workflow.
  size_t stage = READ_TAG;

  // Value parsing.
  SectionState section;
  size_t remaining_skip_length = 0;

  // Collected data.
  std::array<size_t, 16> varint_values;
};

struct InternalState {
  AcDcState ac_dc;
  SectionState section;
  std::unique_ptr<HeaderState> header;

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
};

}  // namespace dec
}  // namespace internal
}  // namespace brunsli

#endif  // BRUNSLI_DEC_STATE_INTERNAL_H_
