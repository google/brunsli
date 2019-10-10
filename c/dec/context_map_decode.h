// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_CONTEXT_MAP_DECODE_H_
#define BRUNSLI_DEC_CONTEXT_MAP_DECODE_H_

#include <brunsli/types.h>
#include "./bit_reader.h"

namespace brunsli {

// Reads the context map from the bit stream. The context map is an array of
// context_map_size histogram ids. The number of different histogram ids is
// given in num_h_trees.
bool DecodeContextMap(int num_h_trees, int context_map_size,
                      uint8_t* context_map, BrunsliBitReader* br);

}  // namespace brunsli

#endif  // BRUNSLI_DEC_CONTEXT_MAP_DECODE_H_
