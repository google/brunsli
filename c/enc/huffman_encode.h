// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_ENC_HUFFMAN_ENCODE_H_
#define BRUNSLI_ENC_HUFFMAN_ENCODE_H_

#include <brunsli/types.h>
#include "./write_bits.h"

namespace brunsli {

// Builds a Huffman tree for the given histogram, and encodes it into storage
// in a format that can be read by HuffmanDecodingData::ReadFromBitstream.
void BuildAndStoreHuffmanTree(const uint32_t* histogram, const size_t length,
                              uint8_t* depth, uint16_t* bits, Storage* storage);

}  // namespace brunsli

#endif  // BRUNSLI_ENC_HUFFMAN_ENCODE_H_
