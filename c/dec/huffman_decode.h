// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library to decode the Huffman code lengths from the bit-stream and build a
// decoding table from them.

#ifndef BRUNSLI_DEC_HUFFMAN_DECODE_H_
#define BRUNSLI_DEC_HUFFMAN_DECODE_H_

#include <memory>
#include <vector>

#include "./huffman_table.h"

namespace brunsli {

struct BrunsliBitReader;

template<typename T>
struct Arena {
  size_t capacity = 0;
  // TODO(eustas): use "char" storage to ensure new[] does not initialize...
  std::unique_ptr<T[]> storage;

  void reserve(size_t limit) {
    if (capacity < limit) {
      capacity = limit;
      storage.reset(new T[capacity]);
    }
  }

  T* data() {
    return storage.get();
  }

  void reset() {
    capacity = 0;
    storage.reset();
  }
};

struct HuffmanDecodingData {
  // Decodes the Huffman code lengths from the bit-stream and fills in the
  // pre-allocated table with the corresponding 2-level Huffman decoding table.
  // |arena| is used as an intermediate output for BuildHuffmanTable.
  // Returns false if the Huffman code lengths can not de decoded.
  bool ReadFromBitStream(size_t alphabet_size, BrunsliBitReader* br,
                         Arena<HuffmanCode>* arena = nullptr);

  uint16_t ReadSymbol(BrunsliBitReader* br) const;

  std::vector<HuffmanCode> table_;
};

}  // namespace brunsli

#endif  // BRUNSLI_DEC_HUFFMAN_DECODE_H_
