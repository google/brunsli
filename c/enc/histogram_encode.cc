// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./histogram_encode.h"

#include <algorithm>

#include "../common/ans_params.h"
#include "../common/histogram.h"
#include "../common/platform.h"
#include "../common/types.h"
#include "./fast_log.h"
#include "./write_bits.h"

namespace brunsli {

// Static Huffman code for encoding histogram length.
static const uint8_t kHistogramLengthBitLengths[ANS_MAX_SYMBOLS - 2] = {
  8, 8, 6, 6, 6, 5, 4, 3, 3, 3, 3, 3, 3, 4, 5, 7,
};
static const uint16_t kHistogramLengthSymbols[ANS_MAX_SYMBOLS - 2] = {
  127, 255, 15, 47, 31, 7, 3, 0, 4, 2, 6, 1, 5, 11, 23, 63,
};

// Static Huffman code for encoding logcounts.
static const  uint8_t kLogCountBitLengths[ANS_LOG_TAB_SIZE + 1] = {
  5, 4, 4, 4, 3, 3, 2, 3, 3, 6, 6,
};
static const uint16_t kLogCountSymbols[ANS_LOG_TAB_SIZE + 1] = {
  15, 3, 11, 7, 2, 6, 0, 1, 5, 31, 63,
};

// Returns the difference between largest count that can be represented and is
// smaller than "count" and smallest representable count larger than "count".
static int SmallestIncrement(int count) {
  BRUNSLI_DCHECK(count > 0);
  int bits = Log2FloorNonZero(count);
  int drop_bits = bits - GetPopulationCountPrecision(bits);
  return (1 << drop_bits);
}

template<bool minimize_error_of_sum> bool RebalanceHistogram(
    const float* targets, int max_symbol, int table_size, int* omit_pos,
    int* counts) {
  BRUNSLI_DCHECK(table_size >= 2);
  int sum = 0;
  float sum_nonrounded = 0.0;
  int remainder_pos = -1;
  int remainder_log = -1;
  // Invariant for minimize_error_of_sum == true:
  // abs(sum - sum_nonrounded)
  //   <= SmallestIncrement(max(targets[])) + max_symbol
  for (int n = 0; n < max_symbol; ++n) {
    if (targets[n] > 0) {
      sum_nonrounded += targets[n];
      counts[n] = static_cast<uint32_t>(targets[n] + .5);  // round
      if (counts[n] == 0) counts[n] = 1;
      if (counts[n] == table_size) counts[n] = table_size - 1;
      // Round the count to the closest nonzero multiple of SmallestIncrement
      // (when minimize_error_of_sum is false) or one of two closest so as to
      // keep the sum as close as possible to sum_nonrounded.
      int inc = SmallestIncrement(counts[n]);
      counts[n] -= counts[n] & (inc - 1);
      const float target =
          minimize_error_of_sum ? (sum_nonrounded - sum) : targets[n];
      if (counts[n] == 0 || (target > counts[n] + inc / 2 &&
                             counts[n] + inc < table_size)) {
        counts[n] += inc;
      }
      sum += counts[n];
      const int count_log = Log2FloorNonZero(counts[n]);
      if (count_log > remainder_log) {
        remainder_pos = n;
        remainder_log = count_log;
      }
    }
  }
  BRUNSLI_DCHECK(remainder_pos != -1);
  counts[remainder_pos] -= sum - table_size;
  *omit_pos = remainder_pos;
  return counts[remainder_pos] > 0;
}


void NormalizeCounts(int* counts, int* omit_pos, const int length,
                     const int precision_bits, int* num_symbols, int* symbols) {
  BRUNSLI_DCHECK(precision_bits > 0);
  const int table_size = 1 << precision_bits;  // target sum / table size
  uint64_t total = 0;
  int max_symbol = 0;
  int symbol_count = 0;
  for (int n = 0; n < length; ++n) {
    total += counts[n];
    if (counts[n] > 0) {
      if (symbol_count < kMaxNumSymbolsForSmallCode) {
        symbols[symbol_count] = n;
      }
      ++symbol_count;
      max_symbol = n + 1;
    }
  }
  *num_symbols = symbol_count;
  if (symbol_count == 0) {
    return;
  }
  if (symbol_count == 1) {
    counts[symbols[0]] = table_size;
    return;
  }
  BRUNSLI_DCHECK(symbol_count <= table_size);

  const float norm = 1.f * table_size / total;
  float targets[ANS_MAX_SYMBOLS];
  for (int n = 0; n < max_symbol; ++n) {
    targets[n] = norm * counts[n];
  }
  if (!RebalanceHistogram<false>(targets, max_symbol, table_size, omit_pos,
                                 counts)) {
    // Use an alternative rebalancing mechanism if the one above failed
    // to create a histogram that is positive wherever the original one was.
    // TODO: it is better to report failure than crash here.
    // TODO: fuzz RebalanceHistogram<true>(...).
    BRUNSLI_CHECK(RebalanceHistogram<true>(targets, max_symbol, table_size,
                                           omit_pos, counts));
  }
}

void EncodeCounts(const int* counts,
                  const int omit_pos,
                  const int num_symbols,
                  const int* symbols,
                  size_t* storage_ix,
                  uint8_t* storage) {
  int max_bits = 5;  // = 1 + Log2Floor(ANS_MAX_SYMBOLS - 1);
  if (num_symbols <= 2) {
    // Small tree marker to encode 1-2 symbols.
    WriteBits(1, 1, storage_ix, storage);
    if (num_symbols == 0) {
      WriteBits(max_bits + 1, 0, storage_ix, storage);
    } else {
      WriteBits(1, num_symbols - 1, storage_ix, storage);
      for (int i = 0; i < num_symbols; ++i) {
        WriteBits(max_bits, symbols[i], storage_ix, storage);
      }
    }
    if (num_symbols == 2) {
      WriteBits(ANS_LOG_TAB_SIZE, counts[symbols[0]], storage_ix, storage);
    }
  } else {
    // Mark non-small tree.
    WriteBits(1, 0, storage_ix, storage);

    int length = 0;
    int logcounts[ANS_MAX_SYMBOLS] = { 0 };
    int omit_log = 0;
    for (int i = 0; i < ANS_MAX_SYMBOLS; ++i) {
      BRUNSLI_DCHECK(counts[i] <= ANS_TAB_SIZE);
      BRUNSLI_DCHECK(counts[i] >= 0);
      if (i == omit_pos) {
        length = i + 1;
      } else if (counts[i] > 0) {
        logcounts[i] = Log2FloorNonZero(counts[i]) + 1;
        length = i + 1;
        if (i < omit_pos) {
          omit_log = std::max(omit_log, logcounts[i] + 1);
        } else {
          omit_log = std::max(omit_log, logcounts[i]);
        }
      }
    }
    logcounts[omit_pos] = omit_log;
    // Since num_symbols >= 3, we know that length >= 3, therefore we encode
    // length - 3 with a static Huffman code.
    WriteBits(kHistogramLengthBitLengths[length - 3],
              kHistogramLengthSymbols[length - 3],
              storage_ix, storage);

    // The logcount values are encoded with a static Huffman code.
    for (int i = 0; i < length; ++i) {
      WriteBits(kLogCountBitLengths[logcounts[i]],
                kLogCountSymbols[logcounts[i]],
                storage_ix, storage);
    }
    for (int i = 0; i < length; ++i) {
      if (logcounts[i] > 1 && i != omit_pos) {
        int bitcount = GetPopulationCountPrecision(logcounts[i] - 1);
        int drop_bits = logcounts[i] - 1 - bitcount;
        BRUNSLI_CHECK((counts[i] & ((1 << drop_bits) - 1)) == 0);
        WriteBits(bitcount, (counts[i] >> drop_bits) - (1 << bitcount),
                  storage_ix, storage);
      }
    }
  }
}

double PopulationCost(const int* data, int total_count) {
  if (total_count == 0) {
    return 7;
  }

  double entropy_bits = total_count * ANS_LOG_TAB_SIZE;
  int histogram_bits = 0;
  int count = 0;
  int length = 0;
  if (total_count > ANS_TAB_SIZE) {
    uint64_t total = total_count;
    for (int i = 0; i < ANS_MAX_SYMBOLS; ++i) {
      if (data[i] > 0) {
        ++count;
        length = i;
      }
    }
    if (count == 1) {
      return 7;
    }
    ++length;
    const uint64_t max0 = (total * length) >> ANS_LOG_TAB_SIZE;
    const uint64_t max1 = (max0 * length) >> ANS_LOG_TAB_SIZE;
    const uint32_t min_base = (total + max0 + max1) >> ANS_LOG_TAB_SIZE;
    total += min_base * count;
    const int64_t kFixBits = 32;
    const int64_t kFixOne = 1LL << kFixBits;
    const int64_t kDescaleBits = kFixBits - ANS_LOG_TAB_SIZE;
    const int64_t kDescaleOne = 1LL << kDescaleBits;
    const int64_t kDescaleMask = kDescaleOne - 1;
    const uint32_t mult = kFixOne / total;
    const uint32_t error = kFixOne % total;
    uint32_t cumul = error;
    if (error < kDescaleOne) {
      cumul += (kDescaleOne - error) >> 1;
    }
    if (data[0] > 0) {
      uint64_t c = (uint64_t)(data[0] + min_base) * mult + cumul;
      double log2count = FastLog2(c >> kDescaleBits);
      entropy_bits -= data[0] * log2count;
      cumul = c & kDescaleMask;
    }
    for (int i = 1; i < length; ++i) {
      if (data[i] > 0) {
        uint64_t c = (uint64_t)(data[i] + min_base) * mult + cumul;
        double log2count = FastLog2(c >> kDescaleBits);
        int log2floor = static_cast<int>(log2count);
        entropy_bits -= data[i] * log2count;
        histogram_bits += log2floor;
        histogram_bits += kLogCountBitLengths[log2floor + 1];
        cumul = c & kDescaleMask;
      } else {
        histogram_bits += kLogCountBitLengths[0];
      }
    }
  } else {
    double log2norm = ANS_LOG_TAB_SIZE - FastLog2(total_count);
    if (data[0] > 0) {
      double log2count = FastLog2(data[0]) + log2norm;
      entropy_bits -= data[0] * log2count;
      length = 0;
      ++count;
    }
    for (int i = 1; i < ANS_MAX_SYMBOLS; ++i) {
      if (data[i] > 0) {
        double log2count = FastLog2(data[i]) + log2norm;
        int log2floor = static_cast<int>(log2count);
        entropy_bits -= data[i] * log2count;
        if (log2floor >= ANS_LOG_TAB_SIZE) {
          log2floor = ANS_LOG_TAB_SIZE - 1;
        }
        histogram_bits += GetPopulationCountPrecision(log2floor);
        histogram_bits += kLogCountBitLengths[log2floor + 1];
        length = i;
        ++count;
      } else {
        histogram_bits += kLogCountBitLengths[0];
      }
    }
    ++length;
  }

  if (count == 1) {
    return 7;
  }

  if (count == 2) {
    return static_cast<int>(entropy_bits) + 1 + 12 + ANS_LOG_TAB_SIZE;
  }

  histogram_bits += kHistogramLengthBitLengths[length - 3];

  return histogram_bits + static_cast<int>(entropy_bits) + 1;
}

}  // namespace brunsli
