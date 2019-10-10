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
#include "./huffman_decode.h"
#include "./huffman_table.h"

namespace brunsli {

namespace {

int ReadHistogramLength(BrunsliBitReader* br) {
  // TODO: direct prefix code decoding would be faster / more readable.
  static const uint8_t kHistogramLengthBitLengths[16] = {
      8, 8, 6, 6, 6, 5, 4, 3, 3, 3, 3, 3, 3, 4, 5, 7,
  };
  HuffmanCode table[256];
  uint16_t counts[16] = {0};
  for (size_t i = 0; i < 16; ++i) {
    ++counts[kHistogramLengthBitLengths[i]];
  }
  BuildHuffmanTable(table, 8, kHistogramLengthBitLengths, 16, &counts[0]);
  const HuffmanCode* p = table;
  BrunsliBitReaderFillWindow(br, 8);
  p += (br->val_ >> br->bit_pos_) & 255;
  br->bit_pos_ += p->bits;
  return p->value + 3;
}

}  // namespace

bool ReadHistogram(int precision_bits, int length, int* counts,
                   BrunsliBitReader* br) {
  if (!BrunsliBitReaderReadMoreInput(br)) {
    BRUNSLI_LOG_DEBUG() << "[ReadHistogram] Unexpected end of input."
                        << BRUNSLI_ENDL();
    return false;
  }
  int simple_code = BrunsliBitReaderReadBits(br, 1);
  if (simple_code == 1) {
    int i;
    int max_bits_counter = length - 1;
    int max_bits = 0;
    int symbols[2] = {0};
    const int num_symbols = BrunsliBitReaderReadBits(br, 1) + 1;
    while (max_bits_counter) {
      max_bits_counter >>= 1;
      ++max_bits;
    }
    memset(counts, 0, length * sizeof(counts[0]));
    for (i = 0; i < num_symbols; ++i) {
      symbols[i] = BrunsliBitReaderReadBits(br, max_bits) % length;
    }
    if (num_symbols == 1) {
      counts[symbols[0]] = 1u << precision_bits;
    } else {
      if (symbols[0] == symbols[1]) {  // corrupt data
        return false;
      }
      counts[symbols[0]] = BrunsliBitReaderReadBits(br, precision_bits);
      counts[symbols[1]] = (1u << precision_bits) - counts[symbols[0]];
    }
  } else {
    int real_length = ReadHistogramLength(br);
    memset(counts, 0, length * sizeof(counts[0]));
    int total_count = 0;
    static const HuffmanCode huff[64] = {
        {2, 6}, {3, 7}, {3, 4}, {4, 1}, {2, 6}, {3, 8}, {3, 5}, {4, 3},
        {2, 6}, {3, 7}, {3, 4}, {4, 2}, {2, 6}, {3, 8}, {3, 5}, {5, 0},
        {2, 6}, {3, 7}, {3, 4}, {4, 1}, {2, 6}, {3, 8}, {3, 5}, {4, 3},
        {2, 6}, {3, 7}, {3, 4}, {4, 2}, {2, 6}, {3, 8}, {3, 5}, {6, 9},
        {2, 6}, {3, 7}, {3, 4}, {4, 1}, {2, 6}, {3, 8}, {3, 5}, {4, 3},
        {2, 6}, {3, 7}, {3, 4}, {4, 2}, {2, 6}, {3, 8}, {3, 5}, {5, 0},
        {2, 6}, {3, 7}, {3, 4}, {4, 1}, {2, 6}, {3, 8}, {3, 5}, {4, 3},
        {2, 6}, {3, 7}, {3, 4}, {4, 2}, {2, 6}, {3, 8}, {3, 5}, {6, 10},
    };
    int log_counts[ANS_MAX_SYMBOLS];
    int omit_log = -1;
    int omit_pos = -1;
    BRUNSLI_DCHECK(real_length > 2);
    for (int i = 0; i < real_length; ++i) {
      const HuffmanCode* p = huff;
      BrunsliBitReaderFillWindow(br, 6);
      p += (br->val_ >> br->bit_pos_) & 63;
      br->bit_pos_ += p->bits;
      log_counts[i] = p->value;
      if (log_counts[i] > omit_log) {
        omit_log = log_counts[i];
        omit_pos = i;
      }
    }
    BRUNSLI_DCHECK(omit_pos >= 0);
    for (int i = 0; i < real_length; ++i) {
      int code = log_counts[i];
      if (i == omit_pos) {
        continue;
      } else if (code == 0) {
        continue;
      } else if (code == 1) {
        counts[i] = 1;
      } else {
        int bit_count = GetPopulationCountPrecision(code - 1);
        counts[i] =
            (1u << (code - 1)) +
            (BrunsliBitReaderReadBits(br, bit_count) << (code - 1 - bit_count));
      }
      total_count += counts[i];
    }
    if (total_count >= (1u << precision_bits)) {
      // The histogram we've read sums to more than total_count (including at
      // least 1 for the omitted value).
      return false;
    }
    counts[omit_pos] = (1u << precision_bits) - total_count;
  }
  return true;
}

}  // namespace brunsli
