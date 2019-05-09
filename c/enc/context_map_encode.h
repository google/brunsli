// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_ENC_CONTEXT_MAP_ENCODE_H_
#define BRUNSLI_ENC_CONTEXT_MAP_ENCODE_H_

#include <vector>

#include "../common/types.h"

namespace brunsli {

// Encodes the given context map to the bit stream. The number of different
// histogram ids is given by num_clusters.
void EncodeContextMap(const std::vector<uint32_t>& context_map,
                      size_t num_clusters,
                      size_t* storage_ix, uint8_t* storage);

}  // namespace brunsli

#endif  // BRUNSLI_ENC_CONTEXT_MAP_ENCODE_H_
