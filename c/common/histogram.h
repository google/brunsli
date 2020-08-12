// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_COMMON_HISTOGRAM_H_
#define BRUNSLI_COMMON_HISTOGRAM_H_

#include <brunsli/types.h>

namespace brunsli {

// Returns the precision (number of bits) that should be used to store
// a histogram count such that Log2Floor(count) == logcount.
inline uint32_t GetPopulationCountPrecision(uint32_t logcount) {
  return (logcount + 1) >> 1;
}

}  // namespace brunsli

#endif  // BRUNSLI_COMMON_HISTOGRAM_H_
