// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <array>
#include <vector>

#include "gtest/gtest.h"
#include <brunsli/types.h>
#include "../dec/huffman_table.h"

namespace brunsli {

void CheckTableSize(std::array<uint16_t, 16> counts, size_t num_symbols,
                    size_t expected_size) {
  std::vector<HuffmanCode> table(expected_size);
  std::vector<uint8_t> code_lengths;

  for (size_t i = 0; i < counts.size(); ++i) {
    for (size_t j = 0; j < counts[i]; ++j) {
      code_lengths.emplace_back(i);
    }
  }

  ASSERT_EQ(num_symbols, code_lengths.size());

  size_t table_size = BuildHuffmanTable(table.data(), 8, code_lengths.data(),
                                        code_lengths.size(), counts.data());
  EXPECT_EQ(table.size(), table_size);
}

TEST(BuildHuffmanTable, Max2048) {
  std::array<uint16_t, 16> counts = {0, 0, 0, 0, 0,  0,  0,  0,
                                     0, 1, 1, 9, 17, 33, 65, 1922};
  CheckTableSize(counts, 2048, 2424);
}

TEST(BuildHuffmanTable, Max2560) {
  std::array<uint16_t, 16> counts = {0, 0, 0, 0, 0,  0,  0,  0,
                                     0, 1, 1, 9, 17, 33, 65, 2434};
  CheckTableSize(counts, 2560, 2936);
}

TEST(BuildHuffmanTable, Max272) {
  std::array<uint16_t, 16> counts = {0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 1, 1, 9, 1, 1, 1, 258};
  CheckTableSize(counts, 272, 648);
}

}  // namespace brunsli
