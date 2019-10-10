// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

/* Bit reading helpers */

#include "./bit_reader.h"

#include <cstring> /* for memset, memcpy */

#include "../common/platform.h"

namespace brunsli {

void BrunsliBitReaderInit(BrunsliBitReader* const br, const uint8_t* buffer,
                          size_t length) {
  br->src_ = buffer;
  br->available_ = length;

  memset(br->tail_, 0, sizeof(br->tail_));
  size_t tail_size =
      length > BRUNSLI_BIT_READER_SLACK ? BRUNSLI_BIT_READER_SLACK : length;
  memcpy(&br->tail_[BRUNSLI_BIT_READER_SLACK - tail_size],
         &buffer[length - tail_size], tail_size);
  // Switch to tail right away, if input is too short.
  int tail_offset = BRUNSLI_BIT_READER_SLACK - br->available_;
  if (tail_offset >= 0) {
    br->src_ = &br->tail_[tail_offset];
  }

  br->val_ = 0;
  for (size_t i = 0; i < sizeof(br->val_); ++i) {
#if (BRUNSLI_64_BITS_LITTLE_ENDIAN)
    br->val_ |= ((uint64_t)(*br->src_)) << (8 * i);
#else
    br->val_ |= ((uint32_t)(*br->src_)) << (8 * i);
#endif
    ++br->src_;
    --br->available_;
  }
  br->bit_pos_ = 0;
}

}  // namespace brunsli
