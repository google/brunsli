// Copyright (c) Google LLC 2020
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./bit_reader.h"

#include "../common/platform.h"

namespace brunsli {

void BrunsliBitReaderInit(BrunsliBitReader* br) {
  br->num_bits_ = 0;
  br->bits_ = 0;
  br->num_debt_bytes_ = 0;
  br->is_healthy = true;
}

void BrunsliBitReaderResume(BrunsliBitReader* br, const uint8_t* buffer,
                              size_t length) {
  br->next_ = buffer;
  br->end_ = buffer + length;
}

/*
   Tries to return "debt" if any, and normalize the state.

   Normal state means that less than 8 bits are held in bit buffer.
   Peeking (BrunsliBitReaderGet) more bits than actually using
   (BrunsliBitReaderDrop) could put bit reader into denormal state.
*/
static BRUNSLI_INLINE void BrunsliBitReaderUnload(BrunsliBitReader* br) {
  // Cancel the overdraft.
  while ((br->num_debt_bytes_ > 0) && (br->num_bits_ >= 8)) {
    br->num_debt_bytes_--;
    br->num_bits_ -= 8;
  }
  // Return unused bytes.
  while (br->num_bits_ >= 8) {
    br->next_--;
    br->num_bits_ -= 8;
  }
  br->bits_ &= BrunsliBitReaderBitMask(br->num_bits_);
}

size_t BrunsliBitReaderSuspend(BrunsliBitReader* br) {
  BrunsliBitReaderUnload(br);
  size_t unused_bytes = br->end_ - br->next_;
  br->next_ = nullptr;
  br->end_ = nullptr;
  return unused_bytes;
}

void BrunsliBitReaderFinish(BrunsliBitReader* br) {
  uint32_t n_bits = br->num_bits_;
  // Likely, did not invoke Suspend before.
  if (n_bits >= 8) {
    br->is_healthy = false;
    return;
  }
  if (n_bits > 0) {
    uint32_t padding_bits = BrunsliBitReaderRead(br, n_bits);
    if (padding_bits != 0) br->is_healthy = false;
  }
}

size_t BrunsliBitReaderIsHealthy(BrunsliBitReader* br) {
  BrunsliBitReaderUnload(br);
  return (br->num_debt_bytes_ == 0) && (br->is_healthy);
}

}  // namespace brunsli
