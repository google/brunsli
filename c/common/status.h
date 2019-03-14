// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_COMMON_STATUS_H_
#define BRUNSLI_COMMON_STATUS_H_

namespace brunsli {

typedef enum {
  BRUNSLI_OK = 0,

  // Used if the input is not representable in the compressed brunsli format,
  // either because it is not a valid JPEG file or if some other limitation
  // is exceeded (e.g. absolute value of coefficients or number of Huffman
  // codes).
  BRUNSLI_NON_REPRESENTABLE,

  BRUNSLI_MEMORY_ERROR,
  BRUNSLI_INVALID_PARAM,

  BRUNSLI_COMPRESSION_ERROR,

  BRUNSLI_INVALID_BRN,
  BRUNSLI_DECOMPRESSION_ERROR,

  BRUNSLI_NOT_ENOUGH_DATA,
} BrunsliStatus;

}  // namespace brunsli

#endif  // BRUNSLI_COMMON_STATUS_H_
