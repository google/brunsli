// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Data structures that represent the contents of a jpeg file.

#ifndef BRUNSLI_COMMON_QUANT_MATRIX_H_
#define BRUNSLI_COMMON_QUANT_MATRIX_H_

#include "./constants.h"
#include <brunsli/types.h>

namespace brunsli {

static const size_t kQFactorBits = 6;
static const size_t kQFactorLimit = 1u << kQFactorBits;

void FillQuantMatrix(bool is_chroma, uint32_t q, uint8_t dst[kDCTBlockSize]);
uint32_t FindBestMatrix(const int* src, bool is_chroma,
                        uint8_t dst[kDCTBlockSize]);

}  // namespace brunsli

#endif  // BRUNSLI_COMMON_QUANT_MATRIX_H_
