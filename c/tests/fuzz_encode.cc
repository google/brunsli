// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <memory>
#include <string>
#include <vector>

#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>
#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>
#include <brunsli/brunsli_encode.h>
#include <brunsli/jpeg_data_reader.h>
#include "./test_utils.h"

// Entry point for LibFuzzer.
extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
  // Encode.
  std::unique_ptr<brunsli::JPEGData> enc_jpg(new brunsli::JPEGData);
  if (!brunsli::ReadJpeg(data, size, brunsli::JPEG_READ_ALL, enc_jpg.get())) {
    return 0;
  }
  size_t enc_output_size = brunsli::GetMaximumBrunsliEncodedSize(*enc_jpg);
  std::vector<uint8_t> enc_output(enc_output_size);
  bool enc_ok = brunsli::BrunsliEncodeJpeg(*enc_jpg, enc_output.data(),
                                           &enc_output_size);
  enc_jpg.reset();
  if (!enc_ok) {
    // It is OK, when regular encoder fails.
    // BrunsliEncodeJpegBypass could be used to wrap "broken" JPEGs.
    return 0;
  }
  enc_output.resize(enc_output_size);

  // Decode.
  brunsli::JPEGData dec_jpg;
  brunsli::BrunsliStatus dec_status;
  dec_status =
      brunsli::BrunsliDecodeJpeg(enc_output.data(), enc_output_size, &dec_jpg);
  if (dec_status != brunsli::BRUNSLI_OK) {
    __builtin_trap();
  }
  std::string dec_output;
  brunsli::JPEGOutput dec_out(brunsli::StringOutputFunction, &dec_output);
  bool dec_ok = brunsli::WriteJpeg(dec_jpg, dec_out);
  if (!dec_ok) {
    __builtin_trap();
  }

  // Compare.
  if (dec_output.size() != size) {
    __builtin_trap();
  }
  const uint8_t* dec_data = reinterpret_cast<const uint8_t*>(dec_output.data());
  for (size_t i = 0; i < size; ++i) {
    if (data[i] != dec_data[i]) {
      __builtin_trap();
    }
  }

  return 0;
}
