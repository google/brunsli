// Copyright (c) Google LLC 2020
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_OUPUT_CHUNK_H_
#define BRUNSLI_DEC_OUPUT_CHUNK_H_

#include <initializer_list>
#include <memory>
#include <vector>

#include <brunsli/types.h>

namespace brunsli {
namespace internal {
namespace dec {

/**
 * A chunk of output data.
 *
 * Data producer creates OutputChunks and adds them to the end output queue.
 * Once control flow leaves the producer code, it is considered that chunk of
 * data is final and can not be changed; to underline this fact |next| is a
 * const-pointer.
 *
 * Data consumer removes OutputChunks from the beginning of the output queue.
 * It is possible to consume OutputChunks partially, by updating |next| and
 * |len|.
 *
 * There are 2 types of output chunks:
 *  - owning: actual data is stored in |buffer| field; producer fills data after
 *    the instance it created; it is legal to reduce |len| to show that not all
 *    the capacity of |buffer| is used
 *  - non-owning: represents the data stored (owned) somewhere else
 */
struct OutputChunk {
  // Non-owning
  template<typename Bytes>
  OutputChunk(Bytes& bytes) : len(bytes.size()) {
    // Deal both with const qualifier and data type.
    const void* src = bytes.data();
    next = reinterpret_cast<const uint8_t*>(src);
  }

  // Non-owning
  OutputChunk(const uint8_t* data, size_t size) : next(data), len(size) {}

  // Owning
  explicit OutputChunk(size_t size = 0) {
    buffer.reset(new std::vector<uint8_t>(size));
    next = buffer->data();
    len = size;
  }

  // Owning
  OutputChunk(std::initializer_list<uint8_t> bytes) {
    buffer.reset(new std::vector<uint8_t>(bytes));
    next = buffer->data();
    len = bytes.size();
  }

  const uint8_t* next;
  size_t len;
  std::unique_ptr<std::vector<uint8_t>> buffer;
};

}  // namespace dec
}  // namespace internal
}  // namespace brunsli

#endif  // BRUNSLI_DEC_OUPUT_CHUNK_H_
