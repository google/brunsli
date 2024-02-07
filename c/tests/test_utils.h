// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_TESTS_TEST_UTILS_H_
#define BRUNSLI_TESTS_TEST_UTILS_H_

#include <cstdint>
#include <cstddef>
#include <string>
#include <tuple>
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

std::vector<std::tuple<std::vector<uint8_t>>> ParseMar(const void* data,
                                                             size_t size);

std::vector<uint8_t> ReadTestData(const std::string& filename);

}  // namespace brunsli

#if !defined(TEST)
#define TEST(A, B)    \
  class A##B##_Test { \
   private:           \
    void TestBody();  \
  };                  \
  void A##B##_Test::TestBody()
#endif

#if !defined(FUZZ_TEST)
struct FuzzTestSink {
  template<typename F>
  FuzzTestSink WithSeeds(F) {
    return *this;
  }
};
#define FUZZ_TEST(A, B) \
  const FuzzTestSink unused##A##B = FuzzTestSink()
#endif

#endif  // BRUNSLI_TESTS_TEST_UTILS_H_
