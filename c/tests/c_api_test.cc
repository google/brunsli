// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <brunsli/decode.h>
#include <brunsli/encode.h>

#include "gtest/gtest.h"

namespace {
  size_t OutputToString(void* data, const uint8_t* buf, size_t count) {
    std::string* output = reinterpret_cast<std::string*>(data);
    output->append(reinterpret_cast<const char*>(buf), count);
    return count;
  }
}  // namespace

TEST(CApiTest, Roundtrip) {
  // An 8x8 pixel JPEG image
  const std::string jpeg(
    "\377\330\377\340\0\20JFIF\0\1\1\1\0`\0`\0\0\377\333\0C\0\b\6\6\a"
    "\6\5\b\a\a\a\t\t\b\n\f\24\r\f\v\v\f\31\22\23\17\24\35\32\37\36\35"
    "\32\34\34 $.' \",#\34\34(7),01444\37'9=82<.342\377\333\0C\1\t\t\t"
    "\f\v\f\30\r\r\0302!\34!22222222222222222222222222222222222222222"
    "222222222\377\300\0\21\b\0\b\0\b\3\1\"\0\2\21\1\3\21\1\377\304\0"
    "\25\0\1\1\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\6\377\304\0\32\20\1\1\1\1"
    "\0\3\0\0\0\0\0\0\0\0\0\0\1\2\0\3\0211A\377\304\0\25\1\1\1\0\0\0\0"
    "\0\0\0\0\0\0\0\0\0\0\6\a\377\304\0\25\21\1\1\0\0\0\0\0\0\0\0\0\0"
    "\0\0\0\0\1\0\377\332\0\f\3\1\0\2\21\3\21\0?\0\216\273\276\267W\322"
    "\232\272U\245\362\253\355_\256fc5\270\0\277\377\331",
    311);

  // Test encoding to brunsli string
  std::string brunsli;
  ASSERT_EQ(1, EncodeBrunsli(jpeg.size(),
      reinterpret_cast<const uint8_t*>(jpeg.data()), &brunsli, OutputToString));

  // Brunsli should have made it smaller
  ASSERT_LT(brunsli.size(), jpeg.size());

  // Test decoding to jpeg string
  std::string jpeg2;
  ASSERT_EQ(1, DecodeBrunsli(brunsli.size(),
      reinterpret_cast<const uint8_t*>(brunsli.data()),
      &jpeg2, OutputToString));

  // Roundtrip should be equal to initial string
  ASSERT_EQ(jpeg, jpeg2);
}
