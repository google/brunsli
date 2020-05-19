// Copyright (c) Google LLC 2020
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <brotli/decode.h>
#include "./state.h"
#include "./state_internal.h"

namespace brunsli {
namespace internal {
namespace dec {

State::State() : internal(new InternalState()) {}

State::State(State&&) = default;

State::~State() {}

MetadataState::~MetadataState() {
  if (brotli != nullptr) {
    BrotliDecoderDestroyInstance(brotli);
    brotli = nullptr;
  }
}

}  // namespace dec
}  // namespace internal
}  // namespace brunsli
