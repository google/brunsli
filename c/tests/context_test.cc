// Copyright (c) Google LLC 2020
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "../common/context.h"

#include "gtest/gtest.h"
#include <brunsli/types.h>

namespace brunsli {

// This is a copy of code that is known to work as expected.
int VanillaACPredictContext(int64_t p_orig) {
  uint64_t p;
  static const int64_t kMaxPred = (1 << (kMaxAverageContext + 1)) - 1;
  if (p_orig >= 0) {
    p = static_cast<uint64_t>(p_orig) << 1;
    if (p > kMaxPred) {
      return kMaxAverageContext;
    }
    return Log2FloorNonZero(++p);
  }
  p = static_cast<uint64_t>(-p_orig) << 1;
  if (p > kMaxPred) {
    return -static_cast<int>(kMaxAverageContext);
  } else {
    return -Log2FloorNonZero(++p);
  }
}

TEST(ContextTest, GoldenACPredictContext) {
  ASSERT_EQ(VanillaACPredictContext(-1024), -kMaxAverageContext);
  ASSERT_EQ(VanillaACPredictContext(1024), kMaxAverageContext);
  for (int i = -4096; i <= 4096; ++i) {
    int vanilla_ctx = VanillaACPredictContext(i);
    size_t vanilla_avg_ctx = std::abs(vanilla_ctx);
    size_t vanilla_sign_ctx = kMaxAverageContext + vanilla_ctx;
    size_t avg_ctx;
    size_t sign_ctx;
    ACPredictContext(i, &avg_ctx, &sign_ctx);
    EXPECT_EQ(vanilla_avg_ctx, avg_ctx);
    EXPECT_EQ(vanilla_sign_ctx, sign_ctx);
  }
}

}  // namespace brunsli
