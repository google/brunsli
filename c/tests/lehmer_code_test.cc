// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "../common/lehmer_code.h"

#include "gtest/gtest.h"
#include "../common/types.h"

namespace brunsli {

// Ironically, we use Lehmer random number generator for testing Lehmer code...
uint32_t LehmerRng(uint32_t* state) {
  uint64_t next = *state;
  next = (next * 16807u) % 0x7fffffffu;
  *state = static_cast<uint32_t>(next);
  return *state;
}

TEST(CommonTest, TestPermutation) {
  const size_t kSize = 1000;
  const size_t kNumIterations = 10;
  for (int nbits = 0; nbits < 8; ++nbits) {
    const int max_value = 1 << nbits;
    uint32_t seed = 624;
    std::vector<int> values(kSize);
    for (size_t iter = 0; iter < kNumIterations; ++iter) {
      for (size_t n = 0; n < values.size(); ++n) {
        values[n] = (LehmerRng(&seed) % (256 + 10)) - 10;
      }
      std::vector<int> codes;
      std::vector<int> code_num_bits;
      std::vector<unsigned char> good_values;

      // Compression.
      {
        PermutationCoder p(nbits);
        ASSERT_EQ(p.num_bits(), nbits);

        std::vector<int> value_already_used(max_value, false);
        std::vector<int> bits(kSize);
        int last_num_bits = 10000;
        for (size_t n = 0; n < kSize; ++n) {
          const int v = values[n];
          int cur_code, cur_num_bits;
          const bool ok = p.RemoveValue(v, &cur_code, &cur_num_bits);
          if (!ok) {
            ASSERT_TRUE(v >= max_value || v < 0 || value_already_used[v]);
          } else {
            ASSERT_TRUE(cur_num_bits <= last_num_bits)
                << "number of code bits should decrease.";
            ASSERT_GT(1 << cur_num_bits, cur_code);
            last_num_bits = cur_num_bits;
            codes.push_back(cur_code);
            good_values.push_back(v);
            code_num_bits.push_back(cur_num_bits);
            value_already_used[v] = true;
          }
        }
      }
      // Decompression.
      {
        PermutationCoder p(nbits);

        ASSERT_EQ(p.Remove(-1), -1);
        ASSERT_EQ(p.Remove(10000), -1);
        ASSERT_EQ(p.Remove(max_value), -1);

        for (size_t n = 0; n < codes.size(); ++n) {
          ASSERT_EQ(p.num_bits(), code_num_bits[n]);
          const int v = p.Remove(codes[n]);
          ASSERT_EQ(v, good_values[n]);
        }
      }
    }
  }
}

}  // namespace brunsli
