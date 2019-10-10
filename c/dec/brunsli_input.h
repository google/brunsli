// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_BRUNSLI_INPUT_H_
#define BRUNSLI_DEC_BRUNSLI_INPUT_H_

#include <brunsli/types.h>

namespace brunsli {

static const int kBitMask[] = {
    0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767,
};

struct BrunsliInput {
  BrunsliInput(const uint8_t* data, size_t len)
      : data_(reinterpret_cast<const uint16_t*>(data)),
        len_(len >> 1),
        pos_(0),
        val_(0),
        bit_pos_(0),
        error_(len & 1) {}

  void InitBitReader() { val_ = GetNextWord(); }

  uint16_t GetNextWord() {
    uint16_t val = 0;
    if (pos_ < len_) {
      val = data_[pos_];
    } else {
      error_ = 1;
    }
    ++pos_;
    return val;
  }

  int ReadBits(int nbits) {
    if (bit_pos_ + nbits > 16) {
      uint32_t new_bits = GetNextWord();
      val_ |= new_bits << 16;
    }
    int result = (val_ >> bit_pos_) & kBitMask[nbits];
    bit_pos_ += nbits;
    if (bit_pos_ > 16) {
      bit_pos_ -= 16;
      val_ >>= 16;
    }
    return result;
  }

  const uint16_t* data_;
  size_t len_;
  size_t pos_;
  uint32_t val_;
  int bit_pos_;
  int error_;
};

}  // namespace brunsli

#endif  // BRUNSLI_DEC_BRUNSLI_INPUT_H_
