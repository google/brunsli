// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./lehmer_code.h"

#include <utility>
#include <vector>

namespace brunsli {

void ComputeLehmerCode(const uint32_t* sigma, const size_t len,
                       uint32_t* code) {
  std::vector<uint32_t> items(len);
  for (size_t i = 0; i < len; ++i) items[i] = static_cast<uint32_t>(i);
  for (size_t i = 0; i < len; ++i) {
    std::vector<uint32_t>::iterator it =
        std::find(items.begin(), items.end(), sigma[i]);
    BRUNSLI_DCHECK(it != items.end());
    code[i] = static_cast<uint32_t>(it - items.begin());
    items.erase(it);
  }
}

bool DecodeLehmerCode(const uint32_t* code, size_t len, uint32_t* sigma) {
  std::vector<uint32_t> items(len);
  for (size_t i = 0; i < len; ++i) items[i] = static_cast<uint32_t>(i);
  for (size_t i = 0; i < len; ++i) {
    uint32_t index = code[i];
    if (index >= items.size()) return false;
    const uint32_t value = items[index];
    items.erase(items.begin() + index);
    sigma[i] = value;
  }
  return true;
}

}  // namespace brunsli
