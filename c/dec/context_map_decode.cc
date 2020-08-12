// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./context_map_decode.h"

#include "../common/constants.h"
#include "../common/platform.h"
#include <brunsli/status.h>
#include <brunsli/types.h>
#include "./bit_reader.h"
#include "./huffman_decode.h"

namespace brunsli {

namespace {

void MoveToFront(uint8_t* v, uint8_t index) {
  uint8_t value = v[index];
  uint8_t i = index;
  for (; i; --i) v[i] = v[i - 1];
  v[0] = value;
}

void InverseMoveToFrontTransform(uint8_t* v, size_t v_len) {
  uint8_t mtf[256];
  for (size_t i = 0; i < 256; ++i) {
    mtf[i] = static_cast<uint8_t>(i);
  }
  for (size_t i = 0; i < v_len; ++i) {
    uint8_t index = v[i];
    v[i] = mtf[index];
    if (index) MoveToFront(mtf, index);
  }
}

}  // namespace

BrunsliStatus DecodeContextMap(const HuffmanDecodingData& entropy,
                               size_t max_run_length_prefix, size_t* index,
                               std::vector<uint8_t>* context_map,
                               BrunsliBitReader* br) {
  size_t& i = *index;
  uint8_t* map = context_map->data();
  const size_t length = context_map->size();
  while (i < length) {
    // Check there is enough deta for Huffman code, RLE and IMTF bit.
    if (!BrunsliBitReaderCanRead(br, 15 + max_run_length_prefix + 1)) {
      return BRUNSLI_NOT_ENOUGH_DATA;
    }
    uint32_t code = entropy.ReadSymbol(br);
    if (code == 0) {
      map[i] = 0;
      ++i;
    } else if (code <= max_run_length_prefix) {
      size_t reps = 1u + (1u << code) + (int)BrunsliBitReaderRead(br, code);
      while (--reps) {
        if (i >= length) return BRUNSLI_INVALID_BRN;
        map[i] = 0;
        ++i;
      }
    } else {
      map[i] = (uint8_t)(code - max_run_length_prefix);
      ++i;
    }
  }
  if (BrunsliBitReaderRead(br, 1)) {
    InverseMoveToFrontTransform(map, length);
  }
  return BrunsliBitReaderIsHealthy(br) ? BRUNSLI_OK : BRUNSLI_INVALID_BRN;
}

}  // namespace brunsli
