// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

/* Bit reading helpers */

#ifndef BRUNSLI_DEC_BIT_READER_H_
#define BRUNSLI_DEC_BIT_READER_H_

#include "../common/platform.h"
#include <brunsli/types.h>

namespace brunsli {

#define BRUNSLI_BIT_READER_SLACK (32 + 8)

/* Masking with this expression turns to a single "Unsigned Bit Field Extract"
   UBFX instruction on ARM. */
static BRUNSLI_INLINE uint32_t BrunsliBitMask(uint32_t n) {
  return ~((0xFFFFFFFFu) << n);
}

// TODO: remove this definition when it becomes unused.
#if (BRUNSLI_64_BITS && BRUNSLI_LITTLE_ENDIAN)
#define BRUNSLI_64_BITS_LITTLE_ENDIAN 1
#else
#define BRUNSLI_64_BITS_LITTLE_ENDIAN 0
#endif

typedef struct {
#if (BRUNSLI_64_BITS_LITTLE_ENDIAN)
  uint64_t val_; /* pre-fetched bits */
#else
  uint32_t val_; /* pre-fetched bits */
#endif
  const uint8_t* src_; /* next input */
  int available_;      /* number of bytes available to read */
  uint32_t bit_pos_;   /* current bit-reading position in val_ */

  /* Tail is used as data source when input is almost depleted. At most
   * BRUNSLI_BIT_READER_SLACK are user input, the remaining are zeros to make
   * a guarantee that BRUNSLI_BIT_READER_SLACK bytes would be safely read,
   * plus extra space for bulk (uint64_t) reads. */
  uint8_t tail_[BRUNSLI_BIT_READER_SLACK * 2];
} BrunsliBitReader;

/* Initializes the bit-reader fields. Returns 0 in case of failure. */
void BrunsliBitReaderInit(BrunsliBitReader* br, const uint8_t* buffer,
                          size_t length);

/*
 * Reload up to 32 bits byte-by-byte.
 * This function works on both little- and big-endian.
 */
static BRUNSLI_INLINE void BrunsliShiftBytes32(BrunsliBitReader* br) {
  while (br->bit_pos_ >= 8) {
    br->val_ >>= 8u;
    br->val_ |= ((uint32_t)*br->src_) << 24u;
    ++br->src_;
    --br->available_;
    br->bit_pos_ -= 8;
  }
}

/* Prepare for further input reading. */
static BRUNSLI_INLINE int BrunsliBitReaderReadMoreInput(
    BrunsliBitReader* const br) {
  if (br->available_ <= -8) {
    return 0;
  }
  int tail_offset = BRUNSLI_BIT_READER_SLACK - br->available_;
  if (tail_offset >= 0) {
    br->src_ = &br->tail_[tail_offset];
  }
  return 1;
}

/* Guarantees that there are at least n_bits in the buffer.
   n_bits should be in the range [1..24] */
static BRUNSLI_INLINE void BrunsliBitReaderFillWindow(
    BrunsliBitReader* const br, int n_bits) {
#if (BRUNSLI_64_BITS_LITTLE_ENDIAN)
  // We guarantee 32 available bits, no matter how much is requested.
  (void)n_bits;
  if (br->bit_pos_ >= 32) {
    br->val_ >>= 32u;
    br->bit_pos_ ^= 32u; /* here same as -= 32 because of the if condition */
    br->val_ |= ((uint64_t)(*((const uint32_t*)br->src_))) << 32u;
    br->available_ -= 4;
    br->src_ += 4;
  }
#elif (BRUNSLI_LITTLE_ENDIAN)
  if (br->bit_pos_ >= 16) {
    br->val_ >>= 16;
    br->bit_pos_ ^= 16; /* here same as -= 16 because of the if condition */
    br->val_ |= ((uint32_t)(*((const uint16_t*)br->src_))) << 16;
    br->available_ -= 2;
    br->src_ += 2;
  }
  if (!BRUNSLI_IS_CONSTANT(n_bits) || (n_bits > 16)) {
    if (br->bit_pos_ >= 8) {
      br->val_ >>= 8;
      br->bit_pos_ ^= 8; /* here same as -= 8 because of the if condition */
      br->val_ |= ((uint32_t)(*br->src_)) << 24;
      --br->available_;
      ++br->src_;
    }
  }
#else
  BrunsliShiftBytes32(br);
#endif
}

/* Reads the specified number of bits from Read Buffer. */
static BRUNSLI_INLINE uint32_t
BrunsliBitReaderReadBits(BrunsliBitReader* const br, uint32_t n_bits) {
  uint32_t val;
  BrunsliBitReaderFillWindow(br, n_bits);
  val = (uint32_t)(br->val_ >> br->bit_pos_) & BrunsliBitMask(n_bits);
  BRUNSLI_LOG_DEBUG() << "[BrunsliReadBits]  " << br->bit_pos_ << " "
                      << std::setw(2) << n_bits << "  val: " << std::setw(6)
                      << std::hex << val << BRUNSLI_ENDL();
  br->bit_pos_ += (uint32_t)n_bits;
  return val;
}

/* Returns the number of bytes available skipping bits. */
static BRUNSLI_INLINE int BrunsliBitReaderJumpToByteBoundary(
    BrunsliBitReader* br) {
  uint32_t nbits = br->bit_pos_ & 7u;
  if (nbits > 0) {
    BrunsliBitReaderReadBits(br, 8 - nbits);
  }
  return br->available_ + sizeof(br->val_) - (br->bit_pos_ >> 3u);
}

}  // namespace brunsli

#endif  // BRUNSLI_DEC_BIT_READER_H_
