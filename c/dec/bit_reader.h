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

struct BrunsliBitReader {
  const uint8_t* next_;
  const uint8_t* end_;
  uint32_t num_bits_;
  uint32_t bits_;
  /*
     Number of "virtual" zero bytes located after the end of buffer being
     put into the bit buffer.

     The concept of debt is useful for two purposes:
      - postpone handling of "unexpected end of input"
      - allow efficient peeking of input near the end

     It is OK to have debt, while "balance" is non-negative.
     BrunsliBitReaderUnload returns debt, if possible.
   */
  uint32_t num_debt_bytes_;
  bool is_healthy_;
  bool is_optimistic_;
};

/**
 * Prepare instance.
 *
 * The instance lifecycle looks like:
 *  - Init
 *  - Resume
 *  - (read bits)
 *  - Suspend
 *  - (optionally, go back to Resume stage)
 *  - Finish
 */
void BrunsliBitReaderInit(BrunsliBitReader* br);
/**
 * Supply instance with new chunk of input.
 */
void BrunsliBitReaderResume(BrunsliBitReader* br, const uint8_t* buffer,
                          size_t length);
/**
 * Returns the number of unused bytes.
 */
size_t BrunsliBitReaderSuspend(BrunsliBitReader* br);
/**
 * Drops unused bits of the last used byte.
 *
 * Marks instance as unhealthy, if unused bits are not all 0.
 */
void BrunsliBitReaderFinish(BrunsliBitReader* br);

bool BrunsliBitReaderIsHealthy(BrunsliBitReader* br);

static BRUNSLI_INLINE uint32_t BrunsliBitReaderBitMask(uint32_t n) {
  return ~((0xFFFFFFFFu) << n);
}

/* Internal. */
static BRUNSLI_INLINE void BrunsliBitReaderOweByte(BrunsliBitReader* br) {
  br->num_bits_ += 8;
  br->num_debt_bytes_++;
}

/* Internal. */
static BRUNSLI_INLINE void BrunsliBitReaderMaybeFetchByte(BrunsliBitReader* br,
                                                          uint32_t n_bits) {
  if (br->num_bits_ < n_bits) {
    if (BRUNSLI_PREDICT_FALSE(br->next_ >= br->end_)) {
      BrunsliBitReaderOweByte(br);
    } else {
      br->bits_ |= static_cast<uint32_t>(*br->next_) << br->num_bits_;
      br->num_bits_ += 8;
      br->next_++;
    }
  }
}

/**
 * Turns instance to optimistic mode.
 *
 * The only difference with regular operating mode is that CanRead always
 * returns "true". Use this mode only when instance is supplied with complete
 * input.
 */
void BrunsliBitReaderSetOptimistic(BrunsliBitReader* br);

/**
 * Returns true if there is enough input.
 *
 * Guaranteed to work correctly only in normalized state.
 */
bool BrunsliBitReaderCanRead(BrunsliBitReader* br, size_t n_bits);

static BRUNSLI_INLINE uint32_t BrunsliBitReaderGet(BrunsliBitReader* br,
                                                   uint32_t n_bits) {
  BRUNSLI_DCHECK(n_bits <= 24);
  BrunsliBitReaderMaybeFetchByte(br, n_bits);
  if (n_bits > 8) {
    BrunsliBitReaderMaybeFetchByte(br, n_bits);
    if (n_bits > 16) BrunsliBitReaderMaybeFetchByte(br, n_bits);
  }
  return br->bits_ & BrunsliBitReaderBitMask(n_bits);
}

static BRUNSLI_INLINE void BrunsliBitReaderDrop(BrunsliBitReader* br,
                                                uint32_t n_bits) {
  BRUNSLI_DCHECK(n_bits <= br->num_bits_);
  br->bits_ >>= n_bits;
  br->num_bits_ -= n_bits;
}

BRUNSLI_INLINE uint32_t BrunsliBitReaderRead(BrunsliBitReader* br,
                                                    uint32_t n_bits) {
  uint32_t result = BrunsliBitReaderGet(br, n_bits);
  BrunsliBitReaderDrop(br, n_bits);
  return result;
}

}  // namespace brunsli

#endif  // BRUNSLI_DEC_BIT_READER_H_
