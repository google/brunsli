// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./jpeg_huffman_decode.h"

#include "../common/jpeg_data.h"

namespace brunsli {

// Returns the table width of the next 2nd level table, count is the histogram
// of bit lengths for the remaining symbols, len is the code length of the next
// processed symbol.
static inline int NextTableBitSize(const int* count, int len) {
  int left = 1 << (len - kJpegHuffmanRootTableBits);
  while (len < kJpegHuffmanMaxBitLength) {
    left -= count[len];
    if (left <= 0) break;
    ++len;
    left <<= 1;
  }
  return len - kJpegHuffmanRootTableBits;
}

void BuildJpegHuffmanTable(const int* count_in, const int* symbols,
                           HuffmanTableEntry* lut) {
  HuffmanTableEntry code;    // current table entry
  HuffmanTableEntry* table;  // next available space in table
  int len;         // current code length
  int idx;         // symbol index
  int key;         // prefix code
  int reps;        // number of replicate key values in current table
  int low;         // low bits for current root entry
  int table_bits;  // key length of current table
  int table_size;  // size of current table

  // Make a local copy of the input bit length histogram.
  int count[kJpegHuffmanMaxBitLength + 1] = { 0 };
  int total_count = 0;
  for (len = 1; len <= kJpegHuffmanMaxBitLength; ++len) {
    count[len] = count_in[len];
    total_count += count[len];
  }

  table = lut;
  table_bits = kJpegHuffmanRootTableBits;
  table_size = 1 << table_bits;

  // Special case code with only one value.
  if (total_count == 1) {
    code.bits = 0;
    code.value = symbols[0];
    for (key = 0; key < table_size; ++key) {
      table[key] = code;
    }
    return;
  }

  // Fill in root table.
  key = 0;
  idx = 0;
  for (len = 1; len <= kJpegHuffmanRootTableBits; ++len) {
    for (; count[len] > 0; --count[len]) {
      code.bits = len;
      code.value = symbols[idx++];
      reps = 1 << (kJpegHuffmanRootTableBits - len);
      while (reps--) {
        table[key++] = code;
      }
    }
  }

  // Fill in 2nd level tables and add pointers to root table.
  table += table_size;
  table_size = 0;
  low = 0;
  for (len = kJpegHuffmanRootTableBits + 1;
       len <= kJpegHuffmanMaxBitLength; ++len) {
    for (; count[len] > 0; --count[len]) {
      // Start a new sub-table if the previous one is full.
      if (low >= table_size) {
        table += table_size;
        table_bits = NextTableBitSize(count, len);
        table_size = 1 << table_bits;
        low = 0;
        lut[key].bits = table_bits + kJpegHuffmanRootTableBits;
        lut[key].value = (table - lut) - key;
        ++key;
      }
      code.bits = len - kJpegHuffmanRootTableBits;
      code.value = symbols[idx++];
      reps = 1 << (table_bits - code.bits);
      while (reps--) {
        table[low++] = code;
      }
    }
  }
}

}  // namespace brunsli
