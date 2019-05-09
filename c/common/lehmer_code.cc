// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./lehmer_code.h"

#include <vector>

namespace brunsli {

int FindIndexAndRemove(int val, int* s, int len) {
  int idx = 0;
  for (int i = 0; i < len; ++i) {
    if (s[i] == val) {
      s[i] = -1;
      break;
    } else if (s[i] != -1) {
      ++idx;
    }
  }
  return idx;
}

void ComputeLehmerCode(const int* sigma, const int len, int* code) {
  std::vector<int> stdorder(len);
  for (int i = 0; i < len; ++i) {
    stdorder[i] = i;
  }
  for (int i = 0; i < len; ++i) {
    code[i] = FindIndexAndRemove(sigma[i], &stdorder[0], len);
  }
}

int FindValueAndRemove(int idx, int* s, int len) {
  int pos = 0;
  int val = 0;
  for (int i = 0; i < len; ++i) {
    if (s[i] == -1) continue;
    if (pos == idx) {
      val = s[i];
      s[i] = -1;
      break;
    }
    ++pos;
  }
  return val;
}

void DecodeLehmerCode(const int* code, int len, int* sigma) {
  std::vector<int> stdorder(len);
  for (int i = 0; i < len; ++i) {
    stdorder[i] = i;
  }
  for (int i = 0; i < len; ++i) {
    sigma[i] = FindValueAndRemove(code[i], &stdorder[0], len);
  }
}

}  // namespace brunsli
