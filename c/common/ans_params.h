// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Common parameters that are needed for both the ANS entropy encoding and
// decoding methods.

#ifndef BRUNSLI_COMMON_ANS_PARAMS_H_
#define BRUNSLI_COMMON_ANS_PARAMS_H_

namespace brunsli {

#define ANS_LOG_TAB_SIZE 10u
#define ANS_TAB_SIZE (1u << ANS_LOG_TAB_SIZE)
#define ANS_MAX_SYMBOLS 18
#define ANS_SIGNATURE 0x13u  // Initial state, used as CRC.

}  // namespace brunsli

#endif  // BRUNSLI_COMMON_ANS_PARAMS_H_
