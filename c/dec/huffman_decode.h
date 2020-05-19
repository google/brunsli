// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library to decode the Huffman code lengths from the bit-stream and build a
// decoding table from them.

#ifndef BRUNSLI_DEC_HUFFMAN_DECODE_H_
#define BRUNSLI_DEC_HUFFMAN_DECODE_H_

#include <vector>

#include "./huffman_table.h"

namespace brunsli {

static const size_t kHuffmanTableMask = 0xFFu;
static const size_t kHuffmanTableBits = 8u;
static const int kMaxHuffmanTableSize = 2048;

struct BrunsliBitReader;

struct HuffmanDecodingData {
  HuffmanDecodingData() : table_(kMaxHuffmanTableSize) {}

  // Decodes the Huffman code lengths from the bit-stream and fills in the
  // pre-allocated table with the corresponding 2-level Huffman decoding table.
  // Returns false if the Huffman code lengths can not de decoded.
  bool ReadFromBitStream(int alphabet_size, BrunsliBitReader* br);

  std::vector<HuffmanCode> table_;
};

}  // namespace brunsli

#endif  // BRUNSLI_DEC_HUFFMAN_DECODE_H_
