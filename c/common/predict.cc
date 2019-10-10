// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./predict.h"

#include <brunsli/jpeg_data.h>

namespace brunsli {

namespace {

int AdaptiveMedian(int w, int n, int nw) {
  const int mx = (w > n) ? w : n;
  const int mn = w + n - mx;
  if (nw > mx) {
    return mn;
  } else if (nw < mn) {
    return mx;
  } else {
    return n + w - nw;
  }
}

}  // namespace

int PredictWithAdaptiveMedian(const coeff_t* coeffs, int x, int y, int stride) {
  const int offset1 = -kDCTBlockSize;
  const int offset2 = -stride;
  const int offset3 = offset2 + offset1;
  if (y != 0) {
    if (x != 0) {
      return AdaptiveMedian(coeffs[offset1], coeffs[offset2], coeffs[offset3]);
    } else {
      return coeffs[offset2];
    }
  } else {
    return x ? coeffs[offset1] : 0;
  }
}

}  // namespace brunsli
