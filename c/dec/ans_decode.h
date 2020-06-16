// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library to decode the ANS population counts from the bit-stream and build a
// decoding table from them.

#ifndef BRUNSLI_DEC_ANS_DECODE_H_
#define BRUNSLI_DEC_ANS_DECODE_H_

#include <vector>

#include "../common/ans_params.h"
#include <brunsli/types.h>
#include "./brunsli_input.h"

namespace brunsli {

typedef struct {
  uint16_t offset_;
  uint16_t freq_;
  uint8_t symbol_;
} ANSSymbolInfo;

struct ANSDecodingData {
  ANSDecodingData() {}

  bool Init(const std::vector<uint32_t>& counts);

  ANSSymbolInfo map_[BRUNSLI_ANS_TAB_SIZE];
};

class ANSDecoder {
 public:
  ANSDecoder() {}

  void Init(WordSource* in) {
    state_ = in->GetNextWord();
    state_ = (state_ << 16u) | in->GetNextWord();
  }

  int ReadSymbol(const ANSDecodingData& code, WordSource* in) {
    const uint32_t res = state_ & (BRUNSLI_ANS_TAB_SIZE - 1);
    const ANSSymbolInfo& s = code.map_[res];
    state_ = s.freq_ * (state_ >> BRUNSLI_ANS_LOG_TAB_SIZE) + s.offset_;
    if (state_ < (1u << 16u)) {
      state_ = (state_ << 16u) | in->GetNextWord();
    }
    return s.symbol_;
  }
  bool CheckCRC() const { return state_ == (BRUNSLI_ANS_SIGNATURE << 16u); }

 private:
  uint32_t state_;
};

}  // namespace brunsli

#endif  // BRUNSLI_DEC_ANS_DECODE_H_
