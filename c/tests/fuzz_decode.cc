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

size_t DiscardOutputFunction(void* data, const uint8_t* buf, size_t count) {
  BRUNSLI_UNUSED(data);
  BRUNSLI_UNUSED(buf);
  return count;
}

// Entry point for LibFuzzer.
extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  brunsli::JPEGOutput out(DiscardOutputFunction, nullptr);
  brunsli::JPEGData jpg;
  brunsli::BrunsliStatus status;
  status = brunsli::BrunsliDecodeJpeg(data, size, &jpg);
  if (status == brunsli::BRUNSLI_OK) {
    brunsli::WriteJpeg(jpg, out);
  }
  return 0;
}
