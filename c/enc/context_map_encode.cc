// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library to encode the brunsli context map.

#include "./context_map_encode.h"

#include <algorithm>
#include <cstring> /* for memset */
#include <vector>

#include "../common/huffman_tree.h"
#include "../common/platform.h"
#include <brunsli/types.h>
#include "./write_bits.h"

namespace brunsli {

namespace {

static const int kCodeLengthCodes = 18;

// We use only the context map alphabet in brunsli, where the maximum alphabet
// size is 256 + 16 = 272. (We can have 256 clusters and 16 run length codes).
static const size_t kMaxAlphabetSize = 272;

void StoreVarLenUint8(size_t n, Storage* storage) {
  if (n == 0) {
    WriteBits(1, 0, storage);
  } else {
    WriteBits(1, 1, storage);
    size_t nbits = Log2FloorNonZero(n);
    WriteBits(3, nbits, storage);
    WriteBits(nbits, n - (1 << nbits), storage);
  }
}

void StoreHuffmanTreeOfHuffmanTreeToBitMask(const int num_codes,
                                            const uint8_t* code_length_bitdepth,
                                            Storage* storage) {
  static const uint8_t kStorageOrder[kCodeLengthCodes] = {
      1, 2, 3, 4, 0, 5, 17, 6, 16, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  // The bit lengths of the Huffman code over the code length alphabet
  // are compressed with the following static Huffman code:
  //   Symbol   Code
  //   ------   ----
  //   0          00
  //   1        1110
  //   2         110
  //   3          01
  //   4          10
  //   5        1111
  static const uint8_t kHuffmanBitLengthHuffmanCodeSymbols[6] = {0, 7, 3,
                                                                 2, 1, 15};
  static const uint8_t kHuffmanBitLengthHuffmanCodeBitLengths[6] = {2, 4, 3,
                                                                    2, 2, 4};

  // Throw away trailing zeros:
  size_t codes_to_store = kCodeLengthCodes;
  if (num_codes > 1) {
    for (; codes_to_store > 0; --codes_to_store) {
      if (code_length_bitdepth[kStorageOrder[codes_to_store - 1]] != 0) {
        break;
      }
    }
  }
  size_t skip_some = 0;  // skips none.
  if (code_length_bitdepth[kStorageOrder[0]] == 0 &&
      code_length_bitdepth[kStorageOrder[1]] == 0) {
    skip_some = 2;  // skips two.
    if (code_length_bitdepth[kStorageOrder[2]] == 0) {
      skip_some = 3;  // skips three.
    }
  }
  WriteBits(2, skip_some, storage);
  for (size_t i = skip_some; i < codes_to_store; ++i) {
    size_t l = code_length_bitdepth[kStorageOrder[i]];
    WriteBits(kHuffmanBitLengthHuffmanCodeBitLengths[l],
              kHuffmanBitLengthHuffmanCodeSymbols[l], storage);
  }
}

void StoreHuffmanTreeToBitMask(const size_t huffman_tree_size,
                               const uint8_t* huffman_tree,
                               const uint8_t* huffman_tree_extra_bits,
                               const uint8_t* code_length_bitdepth,
                               const uint16_t* code_length_bitdepth_symbols,
                               Storage* storage) {
  for (size_t i = 0; i < huffman_tree_size; ++i) {
    size_t ix = huffman_tree[i];
    WriteBits(code_length_bitdepth[ix], code_length_bitdepth_symbols[ix],
              storage);
    // Extra bits
    switch (ix) {
      case 16:
        WriteBits(2, huffman_tree_extra_bits[i], storage);
        break;
      case 17:
        WriteBits(3, huffman_tree_extra_bits[i], storage);
        break;
    }
  }
}

void StoreSimpleHuffmanTree(const uint8_t* depths, size_t symbols[4],
                            size_t num_symbols, size_t max_bits,
                            Storage* storage) {
  // value of 1 indicates a simple Huffman code
  WriteBits(2, 1, storage);
  WriteBits(2, num_symbols - 1, storage);  // NSYM - 1

  // Sort
  for (size_t i = 0; i < num_symbols; i++) {
    for (size_t j = i + 1; j < num_symbols; j++) {
      if (depths[symbols[j]] < depths[symbols[i]]) {
        std::swap(symbols[j], symbols[i]);
      }
    }
  }

  if (num_symbols == 2) {
    WriteBits(max_bits, symbols[0], storage);
    WriteBits(max_bits, symbols[1], storage);
  } else if (num_symbols == 3) {
    WriteBits(max_bits, symbols[0], storage);
    WriteBits(max_bits, symbols[1], storage);
    WriteBits(max_bits, symbols[2], storage);
  } else {
    WriteBits(max_bits, symbols[0], storage);
    WriteBits(max_bits, symbols[1], storage);
    WriteBits(max_bits, symbols[2], storage);
    WriteBits(max_bits, symbols[3], storage);
    // tree-select
    WriteBits(1, depths[symbols[0]] == 1 ? 1 : 0, storage);
  }
}

// num = alphabet size
// depths = symbol depths
void StoreHuffmanTree(const uint8_t* depths, size_t num, Storage* storage) {
  // Write the Huffman tree into the compact representation.
  BRUNSLI_DCHECK(num <= kMaxAlphabetSize);
  uint8_t huffman_tree[kMaxAlphabetSize];
  uint8_t huffman_tree_extra_bits[kMaxAlphabetSize];
  size_t huffman_tree_size = 0;
  WriteHuffmanTree(depths, num, &huffman_tree_size, huffman_tree,
                   huffman_tree_extra_bits);

  // Calculate the statistics of the Huffman tree in the compact representation.
  uint32_t huffman_tree_histogram[kCodeLengthCodes] = {0};
  for (size_t i = 0; i < huffman_tree_size; ++i) {
    ++huffman_tree_histogram[huffman_tree[i]];
  }

  int num_codes = 0;
  int code = 0;
  for (int i = 0; i < kCodeLengthCodes; ++i) {
    if (huffman_tree_histogram[i]) {
      if (num_codes == 0) {
        code = i;
        num_codes = 1;
      } else if (num_codes == 1) {
        num_codes = 2;
        break;
      }
    }
  }

  // Calculate another Huffman tree to use for compressing both the
  // earlier Huffman tree with.
  uint8_t code_length_bitdepth[kCodeLengthCodes] = {0};
  uint16_t code_length_bitdepth_symbols[kCodeLengthCodes] = {0};
  CreateHuffmanTree(&huffman_tree_histogram[0], kCodeLengthCodes, 5,
                    &code_length_bitdepth[0]);
  ConvertBitDepthsToSymbols(code_length_bitdepth, kCodeLengthCodes,
                            &code_length_bitdepth_symbols[0]);

  // Now, we have all the data, let's start storing it
  StoreHuffmanTreeOfHuffmanTreeToBitMask(num_codes, code_length_bitdepth,
                                         storage);

  if (num_codes == 1) {
    code_length_bitdepth[code] = 0;
  }

  // Store the real huffman tree now.
  StoreHuffmanTreeToBitMask(huffman_tree_size, huffman_tree,
                            huffman_tree_extra_bits, &code_length_bitdepth[0],
                            code_length_bitdepth_symbols, storage);
}

size_t IndexOf(const std::vector<uint32_t>& v, uint32_t value) {
  size_t i = 0;
  for (; i < v.size(); ++i) {
    if (v[i] == value) return i;
  }
  return i;
}

void MoveToFront(std::vector<uint32_t>* v, size_t index) {
  uint32_t value = (*v)[index];
  for (size_t i = index; i != 0; --i) {
    (*v)[i] = (*v)[i - 1];
  }
  (*v)[0] = value;
}

std::vector<uint32_t> MoveToFrontTransform(const std::vector<uint32_t>& v) {
  if (v.empty()) return v;
  uint32_t max_value = *std::max_element(v.begin(), v.end());
  std::vector<uint32_t> mtf(max_value + 1);
  for (uint32_t i = 0; i <= max_value; ++i) mtf[i] = i;
  std::vector<uint32_t> result(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    size_t index = IndexOf(mtf, v[i]);
    BRUNSLI_DCHECK(index < mtf.size());
    result[i] = static_cast<uint32_t>(index);
    MoveToFront(&mtf, index);
  }
  return result;
}

// Finds runs of zeros in v_in and replaces them with a prefix code of the run
// length plus extra bits in *v_out and *extra_bits. Non-zero values in v_in are
// shifted by *max_length_prefix. Will not create prefix codes bigger than the
// initial value of *max_run_length_prefix. The prefix code of run length L is
// simply Log2Floor(L) and the number of extra bits is the same as the prefix
// code.
void RunLengthCodeZeros(const std::vector<uint32_t>& v_in,
                        uint32_t* max_run_length_prefix,
                        std::vector<uint32_t>* v_out,
                        std::vector<uint32_t>* extra_bits) {
  size_t max_reps = 0;
  for (size_t i = 0; i < v_in.size();) {
    while (i < v_in.size() && v_in[i] != 0) ++i;
    size_t i0 = i;
    while (i < v_in.size() && v_in[i] == 0) ++i;
    max_reps = std::max(i - i0, max_reps);
  }
  uint32_t max_prefix = max_reps > 0 ? Log2FloorNonZero(max_reps) : 0;
  max_prefix = std::min(max_prefix, *max_run_length_prefix);
  *max_run_length_prefix = max_prefix;
  for (size_t i = 0; i < v_in.size();) {
    if (v_in[i] != 0) {
      v_out->push_back(v_in[i] + *max_run_length_prefix);
      extra_bits->push_back(0);
      ++i;
    } else {
      uint32_t reps = 1;
      for (size_t k = i + 1; k < v_in.size() && v_in[k] == 0; ++k) {
        ++reps;
      }
      i += reps;
      while (reps != 0) {
        if (reps < (2u << max_prefix)) {
          uint32_t run_length_prefix = Log2FloorNonZero(reps);
          v_out->push_back(run_length_prefix);
          extra_bits->push_back(reps - (1u << run_length_prefix));
          break;
        } else {
          v_out->push_back(max_prefix);
          extra_bits->push_back((1u << max_prefix) - 1u);
          reps -= (2u << max_prefix) - 1u;
        }
      }
    }
  }
}

}  // namespace

void BuildAndStoreHuffmanTree(const uint32_t* histogram, const size_t length,
                              uint8_t* depth, uint16_t* bits,
                              Storage* storage) {
  size_t count = 0;
  size_t s4[4] = {0};
  for (size_t i = 0; i < length; i++) {
    if (histogram[i]) {
      if (count < 4) {
        s4[count] = i;
      } else if (count > 4) {
        break;
      }
      count++;
    }
  }

  size_t max_bits_counter = length - 1;
  size_t max_bits = 0;
  while (max_bits_counter) {
    max_bits_counter >>= 1;
    ++max_bits;
  }

  if (count <= 1) {
    WriteBits(4, 1, storage);
    WriteBits(max_bits, s4[0], storage);
    return;
  }

  CreateHuffmanTree(histogram, length, 15, depth);
  ConvertBitDepthsToSymbols(depth, length, bits);

  if (count <= 4) {
    StoreSimpleHuffmanTree(depth, s4, count, max_bits, storage);
  } else {
    StoreHuffmanTree(depth, length, storage);
  }
}

void EncodeContextMap(const std::vector<uint32_t>& context_map,
                      size_t num_clusters, Storage* storage) {
  StoreVarLenUint8(num_clusters - 1, storage);

  if (num_clusters == 1) {
    return;
  }

  std::vector<uint32_t> transformed_symbols = MoveToFrontTransform(context_map);
  std::vector<uint32_t> rle_symbols;
  std::vector<uint32_t> extra_bits;
  uint32_t max_run_length_prefix = 6;
  RunLengthCodeZeros(transformed_symbols, &max_run_length_prefix, &rle_symbols,
                     &extra_bits);
  uint32_t symbol_histogram[kMaxAlphabetSize];
  memset(symbol_histogram, 0, sizeof(symbol_histogram));
  for (size_t i = 0; i < rle_symbols.size(); ++i) {
    ++symbol_histogram[rle_symbols[i]];
  }
  bool use_rle = max_run_length_prefix > 0;
  WriteBits(1, use_rle, storage);
  if (use_rle) {
    WriteBits(4, max_run_length_prefix - 1, storage);
  }
  uint8_t bit_depths[kMaxAlphabetSize];
  uint16_t bit_codes[kMaxAlphabetSize];
  memset(bit_depths, 0, sizeof(bit_depths));
  memset(bit_codes, 0, sizeof(bit_codes));
  BuildAndStoreHuffmanTree(symbol_histogram,
                           num_clusters + max_run_length_prefix, bit_depths,
                           bit_codes, storage);
  for (size_t i = 0; i < rle_symbols.size(); ++i) {
    WriteBits(bit_depths[rle_symbols[i]], bit_codes[rle_symbols[i]], storage);
    if (rle_symbols[i] > 0 && rle_symbols[i] <= max_run_length_prefix) {
      WriteBits(rle_symbols[i], extra_bits[i], storage);
    }
  }
  WriteBits(1, 1, storage);  // use move-to-front
}

}  // namespace brunsli
