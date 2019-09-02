// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library to compute the Lehmer code of a permutation and to reconstruct the
// permutation from its Lehmer code. For more details on Lehmer codes, see
// http://en.wikipedia.org/wiki/Lehmer_code

#ifndef BRUNSLI_COMMON_LEHMER_CODE_H_
#define BRUNSLI_COMMON_LEHMER_CODE_H_

#include <algorithm>
#include <utility>
#include <vector>

#include "./platform.h"
#include <brunsli/types.h>

namespace brunsli {

// Computes the Lehmer code of the permutation sigma[0..len) and puts the
// result into code[0..len).
void ComputeLehmerCode(const int* sigma, int len, int* code);

// Decodes the Lehmer code in code[0..len) and puts the resulting permutation
// into sigma[0..len).
bool DecodeLehmerCode(const int* code, int len, int* sigma);

// This class is an optimized Lehmer-like coder that takes the remaining
// number of possible values into account to reduce the bit usage.
// TODO: in worst case (always removing the first element), O(N^2)
// elements are moved; "Fenwick tree" is simple to implement and could reduce
// the complexity to O(N * log(N)).
class PermutationCoder {
 public:
  explicit PermutationCoder(std::vector<uint8_t> values)
      : values_(std::move(values)) {}
  // number of bits needed to represent the next code.
  int num_bits() const {
    size_t num_values = values_.size();
    BRUNSLI_DCHECK(num_values > 0);
    return num_values <= 1 ? 0 : (Log2FloorNonZero(num_values - 1) + 1);
  }

  // Removes (and return) the value coded by 'code'. Returns -1 in
  // case of error (invalid slot).
  int Remove(int code) {
    if (code >= values_.size() || code < 0) {
      return -1;
    }
    const int value = values_[code];
    values_.erase(values_.begin() + code);
    return value;
  }

  // Removes 'value' from the list and assign a code + number-of-bits
  // for it. Returns false if value could not be encoded.
  bool RemoveValue(uint8_t value, int* code, int* nbits) {
    std::vector<uint8_t>::iterator it =
        std::find(values_.begin(), values_.end(), value);
    if (it == values_.end()) {
      return false;  // invalid/non-existing value was passed.
    }
    *code = it - values_.begin();
    *nbits = num_bits();
    values_.erase(it);
    return true;
  }

 private:
  std::vector<uint8_t> values_;
};

}  // namespace brunsli

#endif  // BRUNSLI_COMMON_LEHMER_CODE_H_
