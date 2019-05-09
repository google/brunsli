// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library to encode the ANS population counts to the bit-stream and encode
// symbols based on the respective distributions.

#ifndef BRUNSLI_ENC_ANS_ENCODE_H_
#define BRUNSLI_ENC_ANS_ENCODE_H_

#include "../common/ans_params.h"
#include "../common/types.h"

namespace brunsli {

// #define USE_MULT_BY_RECIPROCAL

// precision must be equal to:  #bits(state_) + #bits(freq)
#define RECIPROCAL_PRECISION 42

// Data structure representing one element of the encoding table built
// from a distribution.
struct ANSEncSymbolInfo {
  uint16_t freq_;
  uint16_t start_;
#ifdef USE_MULT_BY_RECIPROCAL
  uint64_t ifreq_;
#endif
};

struct ANSTable {
  ANSEncSymbolInfo info_[ANS_MAX_SYMBOLS];
};

class ANSCoder {
 public:
  ANSCoder() : state_(ANS_SIGNATURE << 16) {}

  uint32_t PutSymbol(const ANSEncSymbolInfo t, uint8_t* nbits) {
    uint32_t bits = 0;
    *nbits = 0;
    if ((state_ >> (32 - ANS_LOG_TAB_SIZE)) >= t.freq_) {
      bits = state_ & 0xffff;
      state_ >>= 16;
      *nbits = 16;
    }
#ifdef USE_MULT_BY_RECIPROCAL
    // We use mult-by-reciprocal trick, but that requires 64b calc.
    const uint32_t v = (state_ * t.ifreq_) >> RECIPROCAL_PRECISION;
    const uint32_t offset = state_ - v * t.freq_ + t.start_;
    state_ = (v << ANS_LOG_TAB_SIZE) + offset;
#else
    state_ = ((state_ / t.freq_) << ANS_LOG_TAB_SIZE)
           + (state_ % t.freq_) + t.start_;
#endif
    return bits;
  }

  uint32_t GetState() const { return state_; }

 private:
  uint32_t state_;
};

void BuildAndStoreANSEncodingData(const int* histogram,
                                  ANSTable* table,
                                  size_t* storage_ix, uint8_t* storage);

}  // namespace brunsli

#endif  // BRUNSLI_ENC_ANS_ENCODE_H_
