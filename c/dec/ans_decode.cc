// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./ans_decode.h"

#include <vector>

#include <brunsli/types.h>

namespace brunsli {

bool ANSDecodingData::Init(const std::vector<uint32_t>& counts) {
  size_t pos = 0;
  for (size_t i = 0; i < counts.size(); ++i) {
    for (size_t j = 0; j < counts[i]; ++j, ++pos) {
      map_[pos].symbol_ = static_cast<uint8_t>(i);
      map_[pos].freq_ = static_cast<uint16_t>(counts[i]);
      map_[pos].offset_ = static_cast<uint16_t>(j);
    }
  }
  return (pos == BRUNSLI_ANS_TAB_SIZE);
}

}  // namespace brunsli
