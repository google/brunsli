// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Utility function for building a Huffman lookup table for the jpeg decoder.

#ifndef BRUNSLI_ENC_JPEG_HUFFMAN_DECODE_H_
#define BRUNSLI_ENC_JPEG_HUFFMAN_DECODE_H_

#include <brunsli/types.h>

namespace brunsli {

static const int kJpegHuffmanRootTableBits = 8;
// Maximum huffman lookup table size.
// Requirements: alphabet of 257 symbols (256 + 1 special symbol for the all 1s
// code) and max bit length 16, the root table has 8 bits.
// zlib/examples/enough.c works with an assumption that Huffman code is
// "complete". Input JPEGs might have this assumption broken, hence the
// following sum is used as estimate:
//  + number of 1-st level cells
//  + number of symbols
//  + assymptotic amount of repeated 2-nd level cells
// The third number is 1 + 3 + ... + 255 i.e. it is assumed that sub-table of
// each "size" might be almost completely be filled with repetitions.
// Total sum is slightly less than 1024,...
static const int kJpegHuffmanLutSize = 1024;

struct HuffmanTableEntry {
  // Initialize the value to an invalid symbol so that we can recognize it
  // when reading the bit stream using a Huffman code with space > 0.
  HuffmanTableEntry() : bits(0), value(0xffff) {}

  uint8_t bits;     // number of bits used for this symbol
  uint16_t value;   // symbol value or table offset
};

// Builds jpeg-style Huffman lookup table from the given symbols.
// The symbols are in order of increasing bit lengths. The number of symbols
// with bit length n is given in counts[n] for each n >= 1.
void BuildJpegHuffmanTable(const int* counts, const int* symbols,
                           HuffmanTableEntry* lut);

}  // namespace brunsli

#endif  // BRUNSLI_ENC_JPEG_HUFFMAN_DECODE_H_
