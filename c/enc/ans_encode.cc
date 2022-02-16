// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./ans_encode.h"

#include <vector>

#include <brunsli/types.h>
#include "./histogram_encode.h"
#include "./write_bits.h"

namespace brunsli {

namespace {

void ANSBuildInfoTable(const int* counts, int alphabet_size,
                       ANSEncSymbolInfo info[BRUNSLI_ANS_MAX_SYMBOLS]) {
  int total = 0;
  for (int s = 0; s < alphabet_size; ++s) {
    const uint32_t freq = counts[s];
    info[s].freq_ = counts[s];
    info[s].start_ = total;
    total += freq;
#ifdef USE_MULT_BY_RECIPROCAL
    if (freq != 0) {
      info[s].ifreq_ =
          ((1ull << RECIPROCAL_PRECISION) + info[s].freq_ - 1) / info[s].freq_;
    } else {
      // Shouldn't matter for valid streams (symbol shouldn't occur).
      // Initialization allows to avoid undefined behavior on corrupted streams.
      info[s].ifreq_ = 1;
    }
#endif
  }
}

}  // namespace

void BuildAndStoreANSEncodingData(const int* histogram, ANSTable* table,
                                  Storage* storage) {
  int num_symbols;
  int symbols[kMaxNumSymbolsForSmallCode] = {0};
  std::vector<int> counts(histogram, histogram + BRUNSLI_ANS_MAX_SYMBOLS);
  int omit_pos = 0;  // Initialization is not necessary, but makes MSAN happy.
  NormalizeCounts(&counts[0], &omit_pos, BRUNSLI_ANS_MAX_SYMBOLS,
                  BRUNSLI_ANS_LOG_TAB_SIZE, &num_symbols, symbols);
  ANSBuildInfoTable(&counts[0], BRUNSLI_ANS_MAX_SYMBOLS, table->info_);
  EncodeCounts(&counts[0], omit_pos, num_symbols, symbols, storage);
}

}  // namespace brunsli
