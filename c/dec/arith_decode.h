// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_ARITH_DECODE_H_
#define BRUNSLI_DEC_ARITH_DECODE_H_

#include "../common/distributions.h"
#include "../common/types.h"
#include "./brunsli_input.h"

namespace brunsli {

// A class used for entropy decoding a sequence of binary values.
// skal@ wrote the original version, szabadka@ ported it for brunsli.
class BinaryArithmeticDecoder {
 public:
  BinaryArithmeticDecoder() : low_(0), high_(~0), value_(0) {}

  void Init(BrunsliInput* in) {
    value_ = in->GetNextWord();
    value_ = (value_ << 16) | in->GetNextWord();
  }

  // Returns the next bit decoded from the bit stream, based on the given 8-bit
  // precision probability, i.e. P(bit = 0) = prob / 256. This probability must
  // be the same as the one used by the encoder.
  int ReadBit(int prob, BrunsliInput* in) {
    const uint32_t diff = high_ - low_;
    const uint32_t split = low_ + (((uint64_t)diff * prob) >> 8);
    int bit;
    if (value_ > split) {
      low_ = split + 1;
      bit = 1;
    } else {
      high_ = split;
      bit = 0;
    }
    if (((low_ ^ high_) >> 16) == 0) {
      value_ = (value_ << 16) | in->GetNextWord();
      low_ <<= 16;
      high_ <<= 16;
      high_ |= 0xffff;
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
