// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./context_map_decode.h"

#include <cstring>  /* for memset */
#include <vector>

#include "../common/platform.h"
#include <brunsli/types.h>
#include "./huffman_decode.h"

namespace brunsli {

namespace {

void MoveToFront(uint8_t* v, uint8_t index) {
  uint8_t value = v[index];
  uint8_t i = index;
  for (; i; --i) v[i] = v[i - 1];
  v[0] = value;
}

void InverseMoveToFrontTransform(uint8_t* v, int v_len) {
  uint8_t mtf[256];
  int i;
  for (i = 0; i < 256; ++i) {
    mtf[i] = (uint8_t)i;
  }
  for (i = 0; i < v_len; ++i) {
    uint8_t index = v[i];
    v[i] = mtf[index];
    if (index) MoveToFront(mtf, index);
  }
}

}  // namespace

bool DecodeContextMap(int num_h_trees, int context_map_size,
                      uint8_t* context_map, BrunsliBitReader* br) {
  if (num_h_trees <= 1) {
    memset(context_map, 0, (size_t)context_map_size);
    return true;
  }

  int max_run_length_prefix = 0;
  int use_rle_for_zeros = (int)BrunsliBitReaderReadBits(br, 1);
  if (use_rle_for_zeros) {
    max_run_length_prefix = (int)BrunsliBitReaderReadBits(br, 4) + 1;
  }
  std::vector<HuffmanCode> table(kMaxHuffmanTableSize);
  HuffmanDecodingData entropy;
  if (!entropy.ReadFromBitStream(num_h_trees + max_run_length_prefix, br)) {
    return false;
  }
  for (int i = 0; i < context_map_size;) {
    int code;
    if (!BrunsliBitReaderReadMoreInput(br)) {
      BRUNSLI_LOG_DEBUG() << "[DecodeContextMap] Unexpected end of input."
                          << BRUNSLI_ENDL();
      return false;
    }
    code = HuffmanDecoder::ReadSymbol(entropy, br);
    if (code == 0) {
      context_map[i] = 0;
      ++i;
    } else if (code <= max_run_length_prefix) {
      int reps = 1 + (1u << code) + (int)BrunsliBitReaderReadBits(br, code);
      while (--reps) {
        if (i >= context_map_size) {
          return false;
        }
        context_map[i] = 0;
        ++i;
      }
    } else {
      context_map[i] = (uint8_t)(code - max_run_length_prefix);
      ++i;
    }
  }
  if (BrunsliBitReaderReadBits(br, 1)) {
    InverseMoveToFrontTransform(context_map, context_map_size);
  }
  return true;
}

}  // namespace brunsli
