// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_HUFFMAN_TABLE_H_
#define BRUNSLI_DEC_HUFFMAN_TABLE_H_

#include "../common/types.h"

namespace brunsli {

typedef struct {
  uint8_t bits;     /* number of bits used for this symbol */
  uint16_t value;   /* symbol value or table offset */
} HuffmanCode;

/* Builds Huffman lookup table assuming code lengths are in symbol order. */
/* Returns false in case of error (invalid tree or memory error). */
int BuildHuffmanTable(HuffmanCode* root_table,
                      int root_bits,
                      const uint8_t* const code_lengths,
                      int code_lengths_size,
                      uint16_t* count_arg);

}  // namespace brunsli

#endif  // BRUNSLI_DEC_HUFFMAN_TABLE_H_
