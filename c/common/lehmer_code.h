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
void ComputeLehmerCode(const uint32_t* sigma, size_t len, uint32_t* code);

// Decodes the Lehmer code in code[0..len) and puts the resulting permutation
// into sigma[0..len).
bool DecodeLehmerCode(const uint32_t* code, size_t len, uint32_t* sigma);

// This class is an optimized Lehmer-like coder that takes the remaining
// number of possible values into account to reduce the bit usage.
// TODO(eustas): in worst case (always removing the first element), O(N^2)
// elements are moved; "Fenwick tree" is simple to implement and could reduce
// the complexity to O(N * log(N)).
class PermutationCoder {
 public:
  PermutationCoder() {}

  void Init(std::vector<uint8_t> values) {
    values_ = std::move(values);
  }

  void Clear() {
    std::vector<uint8_t>().swap(values_);
  }

  // number of bits needed to represent the next code.
  int num_bits() const {
    uint32_t num_values = static_cast<uint32_t>(values_.size());
    BRUNSLI_DCHECK(num_values > 0);
    if (num_values <= 1) return 0;
    return static_cast<int>(Log2FloorNonZero(num_values - 1) + 1);
  }

  // Copy value at position 'code' and remove it. Returns false in
  // case of error (invalid slot).
  bool Remove(size_t code, uint8_t* value) {
    if (code >= values_.size()) {
      return false;
    }
    *value = values_[code];
    values_.erase(values_.begin() + code);
    return true;
  }

  // Removes 'value' from the list and assign a code + number-of-bits
  // for it. Returns false if value could not be encoded.
  bool RemoveValue(uint8_t value, int* code, int* nbits) {
    std::vector<uint8_t>::iterator it =
        std::find(values_.begin(), values_.end(), value);
    if (it == values_.end()) {
      return false;  // invalid/non-existing value was passed.
    }
    *code = static_cast<int>(it - values_.begin());
    *nbits = num_bits();
    values_.erase(it);
    return true;
  }

 private:
  std::vector<uint8_t> values_;
};

}  // namespace brunsli

#endif  // BRUNSLI_COMMON_LEHMER_CODE_H_
