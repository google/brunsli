// Copyright (c) Google LLC 2020
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "../common/quant_matrix.h"

#include <vector>
#include <utility>

#include "gtest/gtest.h"
#include <brunsli/types.h>

namespace brunsli {

TEST(QuantMatrixTest, TestFindQ) {
  for (size_t c = 0; c < 2; ++c) {
    for (uint32_t q = 0; q < kQFactorLimit; ++q) {
      uint8_t src[kDCTBlockSize];
      FillQuantMatrix(!!c, q, src);
      int src_extended[kDCTBlockSize];
      for (size_t k = 0; k < kDCTBlockSize; ++k) src_extended[k] = src[k];
      uint8_t dst[kDCTBlockSize];
      uint32_t best_q = FindBestMatrix(src_extended, !!c, dst);
      EXPECT_EQ(q, best_q);
    }
  }
}

}  // namespace brunsli
