// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./write_bits.h"

#include <brunsli/types.h>

#include "../common/platform.h"

namespace brunsli {

Storage::Storage(uint8_t* data, size_t length)
    : data(data), length(length), pos(0) {
  // 8 bytes are necessary for efficient writing (unaligned uint64_t store).
  // Extra 8 bytes are useful to avoid frequent IsHealthy() checks.
  BRUNSLI_CHECK(length > 16);
  data[0] = 0;
  max_pos = (length >= 16) ? ((length - 16) << 3) : 0;
}

bool Storage::AppendBytes(const uint8_t* src, size_t len) {
  BRUNSLI_DCHECK((pos & 7) == 0);
  if (length - GetBytesUsed() < len) return false;
  memcpy(data + (pos >> 3), src, len);
  pos += 8 * len;
  return true;
}

Storage::~Storage() { BRUNSLI_CHECK(GetBytesUsed() <= length); }

}  // namespace brunsli
