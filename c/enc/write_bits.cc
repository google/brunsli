// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "../common/platform.h"
#include "./write_bits.h"

namespace brunsli {

Storage::Storage(uint8_t* data, size_t length)
    : data(data), length(length), pos(0) {
  BRUNSLI_CHECK(length > 0);
  data[0] = 0;
}

void Storage::AppendBytes(const uint8_t* data, size_t len) {
  BRUNSLI_DCHECK((pos & 7) == 0);
  BRUNSLI_DCHECK(GetBytesUsed() + len <= length);
  memcpy(this->data + (pos >> 3), data, len);
  pos += 8 * len;
}

Storage::~Storage() {
  BRUNSLI_CHECK(GetBytesUsed() <= length);
}

}  // namespace brunsli
