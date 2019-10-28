// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "../dec/bit_reader.h"

#include "gtest/gtest.h"
#include <brunsli/types.h>

namespace brunsli {

TEST(BitReader, ReadOneByte) {
  const uint8_t data[2] = {1, 2};
  BrunsliBitReader br;
  BrunsliBitReaderInit(&br, data, 1);

  uint32_t firstByte = BrunsliBitReaderRead(&br, 8);
  ASSERT_EQ(1u, firstByte);
  ASSERT_TRUE(BrunsliBitReaderIsHealthy(&br));

  // It is legal to "peek" after the end of data.
  uint32_t secondByte = BrunsliBitReaderGet(&br, 8);
  ASSERT_EQ(0u, secondByte);
  ASSERT_TRUE(BrunsliBitReaderIsHealthy(&br));

  // It is illegal to read past input, but bit reader gives us some guarantees.
  secondByte = BrunsliBitReaderRead(&br, 8);
  ASSERT_EQ(0u, secondByte);
  ASSERT_FALSE(BrunsliBitReaderIsHealthy(&br));
}

}  // namespace brunsli
