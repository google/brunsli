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
  // Fill with garbage.
  for (size_t i = 0; i < sizeof(br.tail_); ++i) {
    br.tail_[i] = 42;
  }
  BrunsliBitReaderInit(&br, data, 1);

  uint32_t firstByte = BrunsliBitReaderReadBits(&br, 8);
  ASSERT_EQ(1u, firstByte);

  // It is illegal to read past input, but bit reader gives us some guarantees.
  uint32_t secondByte = BrunsliBitReaderReadBits(&br, 8);
  ASSERT_EQ(0u, secondByte);
}

}  // namespace brunsli
