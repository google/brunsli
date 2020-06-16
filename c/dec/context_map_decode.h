// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_CONTEXT_MAP_DECODE_H_
#define BRUNSLI_DEC_CONTEXT_MAP_DECODE_H_

#include <vector>

#include <brunsli/status.h>
#include <brunsli/types.h>

namespace brunsli {

struct BrunsliBitReader;
struct HuffmanDecodingData;

// Reads the context map from the bit stream using the provided entropy decoder.
BrunsliStatus DecodeContextMap(const HuffmanDecodingData& entropy,
                               size_t max_run_length_prefix, size_t* index,
                               std::vector<uint8_t>* context_map,
                               BrunsliBitReader* br);

}  // namespace brunsli

#endif  // BRUNSLI_DEC_CONTEXT_MAP_DECODE_H_
