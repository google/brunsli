// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_HUFFMAN_TABLE_H_
#define BRUNSLI_DEC_HUFFMAN_TABLE_H_

#include <brunsli/types.h>

namespace brunsli {

struct HuffmanCode {
  uint8_t bits;    /* number of bits used for this symbol */
  uint16_t value;  /* symbol value or table offset */
};

/* Builds Huffman lookup table assuming code lengths are in symbol order. */
/* Returns 0 in case of error (invalid tree or memory error), otherwise
   populated size of table. */
uint32_t BuildHuffmanTable(HuffmanCode* root_table, size_t root_bits,
                           const uint8_t* code_lengths,
                           size_t code_lengths_size, uint16_t* count);

}  // namespace brunsli

#endif  // BRUNSLI_DEC_HUFFMAN_TABLE_H_
