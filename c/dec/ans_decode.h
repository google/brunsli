// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library to decode the ANS population counts from the bit-stream and build a
// decoding table from them.

#ifndef BRUNSLI_DEC_ANS_DECODE_H_
#define BRUNSLI_DEC_ANS_DECODE_H_

#include "../common/ans_params.h"
#include "../common/types.h"
#include "./bit_reader.h"
#include "./brunsli_input.h"

namespace brunsli {

typedef struct {
  uint16_t offset_;
  uint16_t freq_;
  uint8_t symbol_;
} ANSSymbolInfo;

struct ANSDecodingData {
  ANSDecodingData() {}

  bool ReadFromBitStream(int alphabet_size, BrunsliBitReader* br);

  ANSSymbolInfo map_[ANS_TAB_SIZE];
};

class ANSDecoder {
 public:
  ANSDecoder() : state_(0) {}

  void Init(BrunsliInput* in) {
    state_ = in->GetNextWord();
    state_ = (state_ << 16) | in->GetNextWord();
  }

  int ReadSymbol(const ANSDecodingData& code, BrunsliInput* in) {
    const uint32_t res = state_ & (ANS_TAB_SIZE - 1);
    const ANSSymbolInfo& s = code.map_[res];
    state_ = s.freq_ * (state_ >> ANS_LOG_TAB_SIZE) + s.offset_;
    if (state_ < (1u << 16)) {
      state_ = (state_ << 16) | in->GetNextWord();
    }
    return s.symbol_;
  }
  bool CheckCRC() const { return state_ == (ANS_SIGNATURE << 16); }

 private:
  uint32_t state_;
};

}  // namespace brunsli

#endif  // BRUNSLI_DEC_ANS_DECODE_H_
