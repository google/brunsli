// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <string>

#include "gtest/gtest.h"
#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>
#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>
#include "./test_utils.h"

namespace brunsli {

TEST(DecodeTest, TestFallback) {
  uint8_t signature[] = {/* marker */ 0x0A, /* length */ 0x04,
      0x42, 0xD2, 0xD5, 0x4E};
  uint8_t header[] = {/* marker */ 0x12, /* length */ 0x02,
      /* submarker */ 0x18, /* version */ 0x04};
  std::string payload = "KOTLOMOMKOLOLSLONA";
  uint8_t payload_size = static_cast<uint8_t>(payload.size());
  uint8_t fallback[] = {/* marker */ 0x4A, /* length */ payload_size};
  std::string encoded =
      std::string(reinterpret_cast<char*>(signature), sizeof(signature)) +
      std::string(reinterpret_cast<char*>(header), sizeof(header)) +
      std::string(reinterpret_cast<char*>(fallback), sizeof(fallback)) +
      payload;

  JPEGData jpg;
  ASSERT_EQ(BRUNSLI_OK,
            BrunsliDecodeJpeg(reinterpret_cast<const uint8_t*>(encoded.data()),
                              encoded.size(), &jpg));

  std::string decoded;
  JPEGOutput out(StringOutputFunction, &decoded);
  ASSERT_TRUE(WriteJpeg(jpg, out));

  ASSERT_EQ(payload, decoded);
}

}  // namespace brunsli
