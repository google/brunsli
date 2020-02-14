// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_ARITH_DECODE_H_
#define BRUNSLI_DEC_ARITH_DECODE_H_

#include "../common/distributions.h"
#include <brunsli/types.h>
#include "./brunsli_input.h"

namespace brunsli {

// A class used for entropy decoding a sequence of binary values.
// skal@ wrote the original version, szabadka@ ported it for brunsli.
class BinaryArithmeticDecoder {
 public:
  BinaryArithmeticDecoder() {}

  void Init(WordSource* in) {
    low_ = 0;
    high_ = ~0u;
    value_ = in->GetNextWord();
    value_ = (value_ << 16u) | in->GetNextWord();
  }

  // Returns the next bit decoded from the bit stream, based on the given 8-bit
  // precision probability, i.e. P(bit = 0) = prob / 256. This probability must
  // be the same as the one used by the encoder.
  int ReadBit(int prob, WordSource* in) {
    const uint32_t diff = high_ - low_;
    const uint32_t split = low_ + (((uint64_t)diff * prob) >> 8u);
    int bit;
    if (value_ > split) {
      low_ = split + 1;
      bit = 1;
    } else {
      high_ = split;
      bit = 0;
    }
    if (((low_ ^ high_) >> 16u) == 0) {
      value_ = (value_ << 16u) | in->GetNextWord();
      low_ <<= 16u;
      high_ <<= 16u;
      high_ |= 0xFFFFu;
    }
    return bit;
  }

 private:
  uint32_t low_;
  uint32_t high_;
  uint32_t value_;
};

}  // namespace brunsli

#endif  // BRUNSLI_DEC_ARITH_DECODE_H_
