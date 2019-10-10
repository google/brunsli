// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library to predict DCT coefficents based on previous blocks.

#ifndef BRUNSLI_COMMON_PREDICT_H_
#define BRUNSLI_COMMON_PREDICT_H_

#include <brunsli/jpeg_data.h>

namespace brunsli {

// Returns a prediction for the DC coefficient in *coeffs, with block
// cooridinates (x,y). The coefficients in the coeffs array are laid out
// block-by-block, i.e. the DC coefficient in the block above is
// *(coeffs - 64 * w).
int PredictWithAdaptiveMedian(const coeff_t* coeffs, int x, int y, int stride);

}  // namespace brunsli

#endif  // BRUNSLI_COMMON_PREDICT_H_
