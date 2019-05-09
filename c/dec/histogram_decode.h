// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library to decode a histogram from the bit-stream.

#ifndef BRUNSLI_DEC_HISTOGRAM_DECODE_H_
#define BRUNSLI_DEC_HISTOGRAM_DECODE_H_

#include "./bit_reader.h"

namespace brunsli {

// Decodes a histogram from the bit-stream where the sum of all population
// counts is 1 << precision_bits.
// Fills in counts[0 .. length) with the decoded population count values.
// Returns false on decoding error.
bool ReadHistogram(int precision_bits, int length, int* counts,
                   BrunsliBitReader* br);

}  // namespace brunsli

#endif  // BRUNSLI_DEC_HISTOGRAM_DECODE_H_
