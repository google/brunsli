// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./quant_matrix.h"

#include "./constants.h"
#include "./platform.h"
#include <brunsli/types.h>

namespace brunsli {

// TODO(eustas): consider high-precision (16-bit) tables in Brunsli v3.
void FillQuantMatrix(bool is_chroma, uint32_t q,
                     uint8_t dst[kDCTBlockSize]) {
  BRUNSLI_DCHECK(q >= 0 && q < kQFactorLimit);
  const uint8_t* const in = kDefaultQuantMatrix[is_chroma];
  for (int i = 0; i < kDCTBlockSize; ++i) {
    const uint32_t v = (in[i] * q + 32) >> 6;
    // clamp to prevent illegal quantizer values
    dst[i] = (v < 1) ? 1 : (v > 255) ? 255u : v;
  }
}

// TODO(eustas): consider high-precision (16-bit) tables in Brunsli v3.
uint32_t FindBestMatrix(const int* src, bool is_chroma,
                        uint8_t dst[kDCTBlockSize]) {
  uint32_t best_q = 0;
  const size_t kMaxDiffCost = 33;
  const size_t kWorstLen = (kDCTBlockSize + 1) * (kMaxDiffCost + 1);
  size_t best_len = kWorstLen;
  for (uint32_t q = 0; q < kQFactorLimit; ++q) {
    FillQuantMatrix(is_chroma, q, dst);
    // Copycat encoder behavior.
    int last_diff = 0;  // difference predictor
    size_t len = 0;
    for (int k = 0; k < kDCTBlockSize; ++k) {
      const int j = kJPEGNaturalOrder[k];
      const int new_diff = src[j] - dst[j];
      int diff = new_diff - last_diff;
      last_diff = new_diff;
      if (diff != 0) {
        len += 1;
        if (diff < 0) diff = -diff;
        diff -= 1;
        if (diff == 0) {
          len++;
        } else if (diff > 65535) {
          len = kWorstLen;
          break;
        } else {
          uint32_t diff_len = Log2FloorNonZero(diff) + 1;
          if (diff_len == 16) diff_len--;
          len += 2 * diff_len + 1;
        }
      }
    }
    if (len < best_len) {
      best_len = len;
      best_q = q;
    }
  }
  FillQuantMatrix(is_chroma, best_q, dst);
  return best_q;
}

}  // namespace brunsli
