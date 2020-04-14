// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <string>

#include "gtest/gtest.h"
#include "../common/constants.h"
#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>
#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>
#include "../dec/state.h"
#include "./test_utils.h"

namespace brunsli {

using ::brunsli::internal::dec::State;

TEST(StreamDecodeTest, DoneDone) {
  std::vector<uint8_t> src = GetSmallBrunsliFile();

  JPEGData jpg;
  State state;
  state.data = src.data();
  state.len = src.size();

  uint8_t foo[] = {42};

  // Decoding is finished.
  ASSERT_EQ(BRUNSLI_OK, internal::dec::ProcessJpeg(&state, &jpg));

  // It is OK to "continue" decoding, result is still "OK"...

  ASSERT_EQ(BRUNSLI_OK, internal::dec::ProcessJpeg(&state, &jpg));

  // ... unless more data is added.
  state.data = foo;
  state.len = 1;
  ASSERT_EQ(BRUNSLI_INVALID_BRN, internal::dec::ProcessJpeg(&state, &jpg));
}

TEST(StreamDecodeTest, ErrorError) {
  std::vector<uint8_t> src = GetSmallBrunsliFile();

  JPEGData jpg;
  State state;
  state.data = src.data();
  state.len = src.size();

  // Once decoder detects corrupted input...
  uint8_t original = src[0];
  src[0] = 42;
  ASSERT_EQ(BRUNSLI_INVALID_BRN, internal::dec::ProcessJpeg(&state, &jpg));

  // ... passing fixed input does not switch decoder to a good state.
  src[0] = original;
  ASSERT_EQ(BRUNSLI_INVALID_BRN, internal::dec::ProcessJpeg(&state, &jpg));
}

TEST(StreamDecodeTest, BytewiseInput) {
  std::vector<uint8_t> src = GetSmallBrunsliFile();

  JPEGData jpg;
  State state;
  size_t start = 0;
  // TODO(eustas): should be src.size() when streaming support is complete.
  size_t last_byte = kSmallBrunsliSignatuteSize + kSmallBrunsliHeaderSize;
  for (size_t end = 0; end < last_byte; ++end) {
    state.data = src.data() + start;
    state.pos = 0;
    state.len = end - start;
    ASSERT_EQ(BRUNSLI_NOT_ENOUGH_DATA,
              internal::dec::ProcessJpeg(&state, &jpg));
    start += state.pos;
  }
}

TEST(StreamDecodeTest, BytewiseFallbackInput) {
  std::vector<uint8_t> src = GetFallbackBrunsliFile();

  JPEGData jpg;
  State state;
  size_t start = 0;
  for (size_t end = 0; end < src.size(); ++end) {
    state.data = src.data() + start;
    state.pos = 0;
    state.len = end - start;
    ASSERT_EQ(BRUNSLI_NOT_ENOUGH_DATA,
              internal::dec::ProcessJpeg(&state, &jpg));
    start += state.pos;
  }
  state.data = src.data() + start;
  state.pos = 0;
  state.len = src.size() - start;
  ASSERT_EQ(BRUNSLI_OK, internal::dec::ProcessJpeg(&state, &jpg));
}
}  // namespace brunsli
