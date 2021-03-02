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
#include <memory>
#include <vector>

// TODO(eustas): cut - used only for "coeff_t*" and "JPEGData*"
#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>

namespace brunsli {

struct ANSDecodingData;

namespace internal {
namespace dec {

typedef std::array<int32_t, kDCTBlockSize> BlockI32;

struct ComponentMeta {
  size_t context_offset;
  int32_t h_samp;
  int32_t v_samp;
  size_t context_bits;
  int32_t ac_stride;
  int32_t b_stride;
  int32_t width_in_blocks;
  int32_t height_in_blocks;
  coeff_t* ac_coeffs;
  // TODO(eustas): investigate bit fields.
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

enum struct SerializationStatus {
  NEEDS_MORE_INPUT,
  NEEDS_MORE_OUTPUT,
  ERROR,
  DONE
};

struct InternalState;

class State {
 public:
  State();
  State(State&&);
  ~State();

  // Public workflow knobs.
  Stage stage = Stage::SIGNATURE;
  // NB: this |tags_met| is not updated by decoder.
  uint32_t tags_met = 0;
  uint32_t skip_tags = 0;

  // Public input knobs.
  const uint8_t* data = nullptr;
  size_t len = 0;
  size_t pos = 0;

  // "JPEGDecodingState" view.
  const uint8_t* context_map;
  const ANSDecodingData* entropy_codes;
  bool use_legacy_context_model = false;

  bool is_storage_allocated = false;
  std::vector<ComponentMeta> meta;

  // Private state parts.
  std::unique_ptr<InternalState> internal;
};

// Use in "headerless" mode, after jpg is filled, but before decoding.
bool UpdateSubsamplingDerivatives(JPEGData* jpg);

// Use in "headerless" mode, after UpdateSubsamplingDerivatives.
void PrepareMeta(const JPEGData* jpg, State* state);

// Use in "groups" mode after views are derived, but before decoding.
// Make sure to set State::is_storage_allocated appropriately.
void WarmupMeta(JPEGData* jpg, State* state);

bool HasSection(const State* state, uint32_t tag);

// Core decoding loop.
BrunsliStatus ProcessJpeg(State* state, JPEGData* jpg);

// Core serialization loop.
SerializationStatus SerializeJpeg(State* state, const JPEGData& jpg,
                                  size_t* available_out, uint8_t** next_out);

}  // namespace dec
}  // namespace internal
}  // namespace brunsli

#endif  // BRUNSLI_DEC_BRUNSLI_DECODE_H_
