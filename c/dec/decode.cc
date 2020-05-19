// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <brunsli/decode.h>

#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>
#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>

/* C API for brunsli encoder */

extern "C" {

struct OutputStruct {
  size_t (*fun)(void* out_data, const uint8_t* buf, size_t size);
  void* data;
};

int DecodeBrunsli(size_t in_size, const uint8_t* in, void* out_data,
                  DecodeBrunsliSink out_fun) {
  OutputStruct out = {out_fun, out_data};
  brunsli::JPEGData jpg;
  brunsli::BrunsliStatus status = brunsli::BrunsliDecodeJpeg(in, in_size, &jpg);
  if (status != brunsli::BRUNSLI_OK) {
    return 0;
  }
  brunsli::JPEGOutput writer(
      [](void* data, const uint8_t* buf, size_t count) {
        OutputStruct* sink = (OutputStruct*)data;
        return sink->fun(sink->data, buf, count);
      },
      &out);
  if (!brunsli::WriteJpeg(jpg, writer)) {
    return 0;
  }
  return 1;  // ok
}

} /* extern "C" */
