// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <brunsli/encode.h>

#include <brunsli/jpeg_data.h>
#include <brunsli/types.h>
#include <brunsli/brunsli_encode.h>
#include <brunsli/jpeg_data_reader.h>

/* C API for brunsli encoder */

extern "C" {

int EncodeBrunsli(size_t insize, const unsigned char* in, void* outdata,
    size_t (*outfun)(void* outdata, const unsigned char* buf, size_t size)) {
  std::vector<uint8_t> output;
  brunsli::JPEGData jpg;
  if (!brunsli::ReadJpeg(in, insize, brunsli::JPEG_READ_ALL, &jpg)) {
    return 0;
  }
  size_t output_size = brunsli::GetMaximumBrunsliEncodedSize(jpg);
  output.resize(output_size);
  if (!brunsli::BrunsliEncodeJpeg(jpg, output.data(), &output_size)) {
    return 0;
  }
  output.resize(output_size);
  if (!outfun(outdata, reinterpret_cast<const unsigned char*>(output.data()),
      output.size())) {
    return 0;
  }
  return 1;  /* ok */
}

}  /* extern "C" */
