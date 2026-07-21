// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <string>

#include "gtest/gtest.h"
#include <brunsli/jpeg_data.h>
#include <brunsli/jpeg_data_reader.h>
#include <brunsli/jpeg_data_writer.h>
#include "./test_utils.h"

namespace brunsli {

// An 8x8 pixel baseline JPEG.
static const char kJpeg[] =
    "\377\330\377\340\0\20JFIF\0\1\1\1\0`\0`\0\0\377\333\0C\0\b\6\6\a"
    "\6\5\b\a\a\a\t\t\b\n\f\24\r\f\v\v\f\31\22\23\17\24\35\32\37\36\35"
    "\32\34\34 $.' \",#\34\34(7),01444\37'9=82<.342\377\333\0C\1\t\t\t"
    "\f\v\f\30\r\r\0302!\34!22222222222222222222222222222222222222222"
    "222222222\377\300\0\21\b\0\b\0\b\3\1\"\0\2\21\1\3\21\1\377\304\0"
    "\25\0\1\1\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\6\377\304\0\32\20\1\1\1\1"
    "\0\3\0\0\0\0\0\0\0\0\0\0\1\2\0\3\0211A\377\304\0\25\1\1\1\0\0\0\0"
    "\0\0\0\0\0\0\0\0\0\0\6\a\377\304\0\25\21\1\1\0\0\0\0\0\0\0\0\0\0"
    "\0\0\0\0\1\0\377\332\0\f\3\1\0\2\21\3\21\0?\0\216\273\276\267W\322"
    "\232\272U\245\362\253\355_\256fc5\270\0\277\377\331";

// A brunsli stream can decode a DC coefficient that no valid JPEG would carry.
// The first block's DC difference equals coeffs[0] (last_dc_coeff starts at 0),
// so a coefficient of INT16_MIN produces a difference whose magnitude is not
// representable. WriteJpeg must reject it rather than compute a 32-bit Huffman
// magnitude category and perform an out-of-range bit shift.
TEST(DctBlockTest, NonRepresentableDcCoeff) {
  std::string jpeg(kJpeg, sizeof(kJpeg) - 1);
  JPEGData jpg;
  ASSERT_TRUE(ReadJpeg(reinterpret_cast<const uint8_t*>(jpeg.data()),
                       jpeg.size(), JPEG_READ_ALL, &jpg));
  ASSERT_FALSE(jpg.components.empty());
  ASSERT_FALSE(jpg.components[0].coeffs.empty());

  jpg.components[0].coeffs[0] = static_cast<coeff_t>(-32768);

  std::string out;
  JPEGOutput sink(StringOutputFunction, &out);
  EXPECT_FALSE(WriteJpeg(jpg, sink));
}

}  // namespace brunsli
