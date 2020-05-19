// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_TESTS_TEST_UTILS_H_
#define BRUNSLI_TESTS_TEST_UTILS_H_

#include <vector>

namespace brunsli {

/**
 * Output callback for JPEGOutput.
 *
 * assert(data instanceof std::string)
 */
size_t StringOutputFunction(void* data, const uint8_t* buf, size_t count);

std::vector<uint8_t> GetSmallBrunsliFile();
const size_t kSmallBrunsliSignatuteSize = 6;
const size_t kSmallBrunsliHeaderSize = 10;

std::vector<uint8_t> GetFallbackBrunsliFile();

}  // namespace brunsli

#endif  // BRUNSLI_TESTS_TEST_UTILS_H_
