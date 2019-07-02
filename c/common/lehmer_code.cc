// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./lehmer_code.h"

#include <utility>
#include <vector>

namespace brunsli {

void ComputeLehmerCode(const int* sigma, const int len, int* code) {
  std::vector<int> items(len);
  for (int i = 0; i < len; ++i) items[i] = i;
  for (int i = 0; i < len; ++i) {
    std::vector<int>::iterator it =
      std::find(items.begin(), items.end(), sigma[i]);
    BRUNSLI_DCHECK(it != items.end());
    code[i] = it - items.begin();
    items.erase(it);
  }
}

bool DecodeLehmerCode(const int* code, int len, int* sigma) {
  std::vector<int> items(len);
  for (int i = 0; i < len; ++i) items[i] = i;
  for (int i = 0; i < len; ++i) {
    int index = code[i];
    if (index >= items.size() || index < 0) return false;
    const int value = items[index];
    items.erase(items.begin() + index);
    sigma[i] = value;
  }
  return true;
}

}  // namespace brunsli
