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

#include "../common/constants.h"
#include "../common/platform.h"
#include <brunsli/types.h>
#include "./huffman_tree.h"
#include "./write_bits.h"

namespace brunsli {

namespace {

void StoreVarLenUint8(size_t n, Storage* storage) {
  if (n == 0) {
    WriteBits(1, 0, storage);
  } else {
    WriteBits(1, 1, storage);
    size_t nbits = Log2FloorNonZero(static_cast<uint32_t>(n));
    WriteBits(3, nbits, storage);
    WriteBits(nbits, n - (size_t(1) << nbits), storage);
  }
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
  uint32_t max_prefix =
      max_reps > 0 ? Log2FloorNonZero(static_cast<uint32_t>(max_reps)) : 0;
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
  uint32_t symbol_histogram[kMaxContextMapAlphabetSize];
  memset(symbol_histogram, 0, sizeof(symbol_histogram));
  for (size_t i = 0; i < rle_symbols.size(); ++i) {
    ++symbol_histogram[rle_symbols[i]];
  }
  bool use_rle = max_run_length_prefix > 0;
  WriteBits(1, use_rle, storage);
  if (use_rle) {
    WriteBits(4, max_run_length_prefix - 1, storage);
  }
  uint8_t bit_depths[kMaxContextMapAlphabetSize];
  uint16_t bit_codes[kMaxContextMapAlphabetSize];
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
