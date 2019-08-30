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

static const int kMaxQFactor = 64;

// TODO: consider high-precision (16-bit) tables in Brunsli v3.
void FillQuantMatrix(bool is_chroma, uint32_t q,
                     uint8_t dst[kDCTBlockSize]) {
  BRUNSLI_DCHECK(q >= 0 && q < kMaxQFactor);
  const uint8_t* const in = kDefaultQuantMatrix[is_chroma];
  for (int i = 0; i < kDCTBlockSize; ++i) {
    const uint32_t v = (in[i] * q + 32) >> 6;
    // clamp to prevent illegal quantizer values
    dst[i] = (v < 1) ? 1 : (v > 255) ? 255u : v;
  }
}

// TODO: consider high-precision (16-bit) tables in Brunsli v3.
uint32_t FindBestMatrix(const int* src, bool is_chroma,
                        uint8_t dst[kDCTBlockSize]) {
  uint32_t best_q = 0;
  float best_err = 274877906944.0f;  // == 64 * 65536 * 65536
  for (uint32_t q = 0; q < kMaxQFactor; ++q) {
    FillQuantMatrix(is_chroma, q, dst);
    float err = 0.0f;
    for (size_t k = 0; k < kDCTBlockSize; ++k) {
      // TODO: is it possible to replace L2 with Linf?
      float delta = src[k] - dst[k];
      err += delta * delta;
      if (err >= best_err) break;
    }
    if (err < best_err) {
      best_err = err;
      best_q = q;
    }
  }
  FillQuantMatrix(is_chroma, best_q, dst);
  return best_q;
}

}  // namespace brunsli
