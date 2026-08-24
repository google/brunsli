// Copyright (c) Google LLC 2020
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./state.h"

#include <brotli/decode.h>
#include <brunsli/brunsli_decode.h>
#include <brunsli/types.h>

#include "./state_internal.h"

namespace brunsli {
namespace internal {
namespace dec {

State::State() : internal(new InternalState()) {}

State::State(State&&) = default;

State::~State() {}

void State::ApplyOptions(const BrunsliDecodeOptions& options) {
  internal->metadata.metadata_size_limit =
      options.max_uncompressed_metadata_size;
}

MetadataState::~MetadataState() {
  if (brotli != nullptr) {
    BrotliDecoderDestroyInstance(brotli);
    brotli = nullptr;
  }
}

bool HasSection(const State* state, uint32_t tag) {
  return state->internal->section.tags_met & (1u << tag);
}

}  // namespace dec
}  // namespace internal
}  // namespace brunsli
