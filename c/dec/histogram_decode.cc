// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./histogram_decode.h"

#include "../common/ans_params.h"
#include "../common/histogram.h"
#include "../common/platform.h"
#include <brunsli/types.h>
#include "./bit_reader.h"

namespace brunsli {

namespace {

// Binary trees; positive - offset to child nodes; other - minus leaf value.
// 3-8 bits encoding value 3..18.
const int8_t kLengthTree[] = {1, 2, 3, 4, 5, 6, 7, -10, -11, -12, -13,
  -14, -15, 2, 3, -9, -16, 2, 3, -8, -17, 2, 3, -5, -6, -7, 1, -18, 1, -3, -4};
// 2..6 bits encoding value 0..10.
const int8_t kLogCountTree[] = {1, 2, 3, -6, 3, 4, 5, -4, -5, -7, -8, 2,
  3, -1, -2, -3, 1, 0, 1, -9, -10};

size_t ReadShortHuffmanCode(BrunsliBitReader* br, const int8_t* tree) {
  size_t pos = 0;
  int8_t delta = 1;
  while (delta > 0) {
    pos += delta + BrunsliBitReaderRead(br, 1);
    delta = tree[pos];
  }
  return static_cast<size_t>(-delta);
}

}  // namespace

bool ReadHistogram(uint32_t precision_bits, std::vector<uint32_t>* counts,
                   BrunsliBitReader* br) {
  BRUNSLI_DCHECK(!counts->empty());
  uint32_t space = 1u << precision_bits;
  const size_t length = counts->size();
  std::fill(counts->begin(), counts->end(), 0);
  uint32_t* histogram = counts->data();
  int simple_code = BrunsliBitReaderRead(br, 1);
  if (simple_code == 1) {
    size_t max_bits_counter = length - 1;
    uint32_t max_bits = 0;
    int symbols[2] = {0};
    const size_t num_symbols = BrunsliBitReaderRead(br, 1) + 1u;
    while (max_bits_counter) {
      max_bits_counter >>= 1;
      ++max_bits;
    }
    for (size_t i = 0; i < num_symbols; ++i) {
      symbols[i] = BrunsliBitReaderRead(br, max_bits) % length;
    }
    if (num_symbols == 1) {
      histogram[symbols[0]] = space;
    } else {
      if (symbols[0] == symbols[1]) {  // corrupt data
        return false;
      }
      uint32_t value = BrunsliBitReaderRead(br, precision_bits);
      histogram[symbols[0]] = value;
      histogram[symbols[1]] = space - value;
    }
  } else {
    size_t real_length = ReadShortHuffmanCode(br, kLengthTree);
    uint32_t total_count = 0;
    uint32_t log_counts[BRUNSLI_ANS_MAX_SYMBOLS];
    size_t omit_pos = 0;
    BRUNSLI_DCHECK(real_length > 2);
    for (size_t i = 0; i < real_length; ++i) {
      log_counts[i] =
          static_cast<uint32_t>(ReadShortHuffmanCode(br, kLogCountTree));
      if (log_counts[i] > log_counts[omit_pos]) omit_pos = i;
    }
    BRUNSLI_DCHECK(omit_pos >= 0);
    for (size_t i = 0; i < real_length; ++i) {
      uint32_t code = log_counts[i];
      if (i == omit_pos) {
        continue;
      } else if (code == 0) {
        continue;
      } else if (code == 1) {
        histogram[i] = 1;
      } else {
        uint32_t bit_count = GetPopulationCountPrecision(code - 1);
        histogram[i] = (1u << (code - 1)) + (BrunsliBitReaderRead(br, bit_count)
                                             << (code - 1 - bit_count));
      }
      total_count += histogram[i];
    }
    if (total_count >= space) {
      // The histogram we've read sums to more than total_count (including at
      // least 1 for the omitted value).
      return false;
    }
    histogram[omit_pos] = space - total_count;
  }
  return BrunsliBitReaderIsHealthy(br);
}

}  // namespace brunsli
