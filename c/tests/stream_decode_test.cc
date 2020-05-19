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

  // TODO(eustas): should be more fine-grained when streaming support is
  //               complete.
  std::vector<size_t> deltas = {
      6,  // signature

      1,  // header tag
      1,  // header size
      1, 1, 1, 1, 1, 1, 1, 1,  // header fields

      1,  // internal tag
      2,  // internal size
      40,  // internal contents

      1,  // metadata tag
      1,  // metadata size
      1,  // uncompressed size
      1, 1, 1, 1, 1,  // metadata compressed contents

      1,  // quant tag
      2,  // quant size
      2,  // quant contents

      1,  // histo tag
      3,  // histo size
      56,  // histo contents

      1,  // unknown tag
      1,  // unknown size
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  // unknown contents
      1, 1, 1,  // unknown contents

      1,  // DC tag
      3,  // DC size
      14,  // DC contents

      1,  // AC tag
      3,  // AC size
      0x142  // AC contents
  };
  size_t delta_idx = 0;

  JPEGData jpg;
  State state;
  size_t start = 0;
  for (size_t end = 0; end <= src.size(); ++end) {
    state.data = src.data() + start;
    state.pos = 0;
    state.len = end - start;
    ASSERT_EQ(end < src.size() ? BRUNSLI_NOT_ENOUGH_DATA : BRUNSLI_OK,
              internal::dec::ProcessJpeg(&state, &jpg));
    size_t delta = state.pos;
    if (delta != 0) {
      ASSERT_LT(delta_idx, deltas.size());
      ASSERT_EQ(deltas[delta_idx], delta);
      delta_idx++;
      start += delta;
    }
  }
  ASSERT_EQ(deltas.size(), delta_idx);
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
