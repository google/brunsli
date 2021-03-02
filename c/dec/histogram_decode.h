// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library to decode a histogram from the bit-stream.

#ifndef BRUNSLI_DEC_HISTOGRAM_DECODE_H_
#define BRUNSLI_DEC_HISTOGRAM_DECODE_H_

#include <vector>

#include <brunsli/types.h>

namespace brunsli {

struct BrunsliBitReader;

/**
 * Decodes a histogram from the bit-stream.
 *
 * The sum of all population |counts| is 1 << |precision_bits|.
 * Fills in |counts| with the decoded population count values.
 *
 * @returns false on decoding error
 */
bool ReadHistogram(uint32_t precision_bits, std::vector<uint32_t>* counts,
                   BrunsliBitReader* br);

}  // namespace brunsli

#endif  // BRUNSLI_DEC_HISTOGRAM_DECODE_H_
