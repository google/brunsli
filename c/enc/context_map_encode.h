// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_ENC_CONTEXT_MAP_ENCODE_H_
#define BRUNSLI_ENC_CONTEXT_MAP_ENCODE_H_

#include <vector>

#include <brunsli/types.h>
#include "./write_bits.h"

namespace brunsli {

// Builds a Huffman tree for the given histogram, and encodes it into storage
// in a format that can be read by HuffmanDecodingData::ReadFromBitstream.
void BuildAndStoreHuffmanTree(const uint32_t* histogram, const size_t length,
                              uint8_t* depth, uint16_t* bits,
                              Storage* storage);

// Encodes the given context map to the bit stream. The number of different
// histogram ids is given by num_clusters.
void EncodeContextMap(const std::vector<uint32_t>& context_map,
                      size_t num_clusters, Storage* storage);

}  // namespace brunsli

#endif  // BRUNSLI_ENC_CONTEXT_MAP_ENCODE_H_
