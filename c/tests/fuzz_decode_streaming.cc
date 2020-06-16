// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <brunsli/jpeg_data.h>
#include "../common/platform.h"
#include <brunsli/status.h>
#include <brunsli/types.h>
#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>
#include "../dec/state.h"

size_t DiscardOutputFunction(void* data, const uint8_t* buf, size_t count) {
  BRUNSLI_UNUSED(data);
  BRUNSLI_UNUSED(buf);
  return count;
}

// Entry point for LibFuzzer.
extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  brunsli::JPEGOutput out(DiscardOutputFunction, nullptr);
  brunsli::JPEGData jpg;
  brunsli::internal::dec::State state;
  brunsli::BrunsliStatus status;

  size_t start = 0;
  for (size_t end = 0; end <= size; ++end) {
    state.data = data + start;
    state.pos = 0;
    state.len = end - start;
    status = brunsli::internal::dec::ProcessJpeg(&state, &jpg);
    brunsli::BrunsliStatus expected_status =
        end < size ? brunsli::BRUNSLI_NOT_ENOUGH_DATA : brunsli::BRUNSLI_OK;
    if (status != expected_status) return 0;
    start += state.pos;
  }

  // TODO(eustas): check that streaming and non-streaming decoders produce the
  //               same output
  brunsli::WriteJpeg(jpg, out);
  return 0;
}
