// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library to store a histogram to the bit-stream.

#ifndef BRUNSLI_ENC_HISTOGRAM_ENCODE_H_
#define BRUNSLI_ENC_HISTOGRAM_ENCODE_H_

#include "../common/ans_params.h"
#include <brunsli/types.h>
#include "./write_bits.h"

namespace brunsli {

static const int kMaxNumSymbolsForSmallCode = 4;

// Normalizes the population counts in counts[0 .. length) so that the sum of
// all counts will be 1 << precision_bits.
// Sets *num_symbols to the number of symbols in the range [0 .. length) with
// non-zero population counts.
// Fills in symbols[0 .. kMaxNumSymbolsForSmallCode) with the first few symbols
// with non-zero population counts.
// Each count will all be rounded to multiples of
// 1 << GetPopulationCountPrecision(count), except possibly for one. The index
// of that count will be stored in *omit_pos.
void NormalizeCounts(int* counts, int* omit_pos, const int length,
                     const int precision_bits, int* num_symbols, int* symbols);

// Stores a histogram in counts[0 .. BRUNSLI_ANS_MAX_SYMBOLS) to the bit-stream
// where the sum of all population counts is BRUNSLI_ANS_TAB_SIZE and the number
// of symbols with non-zero counts is num_symbols. symbols[0 ..
// kMaxNumSymbolsForSmallCode) contains the first few symbols with non-zero
// population counts. Each count must be rounded to a multiple of 1 <<
// GetPopulationCountPrecision(count), except possibly counts[omit_pos].
void EncodeCounts(const int* counts, const int omit_pos, const int num_symbols,
                  const int* symbols, Storage* storage);

// Returns an estimate of the number of bits required to encode the given
// histogram (header bits plus data bits).
double PopulationCost(const int* data, int total_count);

}  // namespace brunsli

#endif  // BRUNSLI_ENC_HISTOGRAM_ENCODE_H_
