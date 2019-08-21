// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <brunsli/decode.h>

#include "../common/jpeg_data.h"
#include "../common/status.h"
#include "../common/types.h"
#include "./brunsli_decode.h"
#include "./jpeg_data_writer.h"

/* C API for brunsli encoder */

extern "C" {

struct OutputStruct {
  size_t (*fun)(void* outdata, const unsigned char* buf, size_t size);
  void* data;
};

int DecodeBrunsli(size_t insize, const unsigned char* in, void* outdata,
    size_t (*outfun)(void* outdata, const unsigned char* buf, size_t size)) {
  OutputStruct out = {outfun, outdata};
  brunsli::JPEGData jpg;
  brunsli::BrunsliStatus status = brunsli::BrunsliDecodeJpeg(
      reinterpret_cast<const uint8_t*>(in), insize,
      brunsli::BRUNSLI_READ_ALL, &jpg, nullptr);
  if (status != brunsli::BRUNSLI_OK) {
    return 0;
  }
  brunsli::JPEGOutput writer([](void* data, const uint8_t* buf, size_t count) {
    OutputStruct* out = (OutputStruct*)data;
    int result = out->fun(out->data, buf, count) == count ? count : -1;
    return result;
  }, &out);
  if (!brunsli::WriteJpeg(jpg, writer)) {
    return 0;
  }
  return 1;  // ok
}

}  /* extern "C" */
