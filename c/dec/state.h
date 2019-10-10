// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Functions for reading a brunsli byte stream into a JPEGData object and
// converting a brunsli byte stream to a jpeg byte stream.

#ifndef BRUNSLI_DEC_STATE_H_
#define BRUNSLI_DEC_STATE_H_

#include <array>
#include <vector>

#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>
#include "./ans_decode.h"

namespace brunsli {
namespace internal {
namespace dec {

typedef std::array<int32_t, kDCTBlockSize> BlockI32;

struct ComponentMeta {
  size_t context_offset;
  int32_t h_samp;
  int32_t v_samp;
  int32_t context_bits;
  int32_t ac_stride;
  int32_t b_stride;
  int32_t width_in_blocks;
  int32_t height_in_blocks;
  coeff_t* ac_coeffs;
  // TODO: investigate bit fields.
  uint8_t* block_state;
  BlockI32 quant;
};

enum struct Stage {
  SIGNATURE = 0,
  HEADER,
  FALLBACK,
  SECTION,
  SECTION_BODY,
  DONE,
  ERROR
};

struct State {
  Stage stage = Stage::SIGNATURE;
  BrunsliStatus result = BRUNSLI_OK;
  int32_t tags_met = 0;
  int32_t skip_tags = 0;

  // "JPEGDecodingState" storage.
  std::vector<uint8_t> context_map_;
  std::vector<ANSDecodingData> entropy_codes_;
  std::vector<std::vector<uint8_t>> block_state_;

  // "JPEGDecodingState" view.
  const uint8_t* context_map;
  const ANSDecodingData* entropy_codes;

  bool is_meta_warm = false;
  bool is_storage_allocated = false;
  std::vector<ComponentMeta> meta;

  // That is not exactly state, but very convenient for passing the input.
  const uint8_t* data = nullptr;
  size_t len = 0;
  size_t pos = 0;

  // Section
  size_t tag = 0;
  size_t section_end = 0;

  // For "estimate peak memory".
  bool shallow_histograms = false;
  size_t num_contexts = 0;
  size_t num_histograms = 0;
};

// Use in "headerless" mode, after jpg is filled, but before decoding.
bool UpdateSubsamplingDerivatives(JPEGData* jpg);

// Use in "headerless" mode, after UpdateSubsamplingDerivatives.
void PrepareMeta(const JPEGData* jpg, State* state);

// Use in "groups" mode after views are derived, but before decoding.
// Make sure to set State::is_storage_allocated appropriately.
void WarmupMeta(JPEGData* jpg, State* state);

// Core decoding loop.
BrunsliStatus ProcessJpeg(State* state, JPEGData* jpg);

}  // namespace dec
}  // namespace internal
}  // namespace brunsli

#endif  // BRUNSLI_DEC_BRUNSLI_DECODE_H_
