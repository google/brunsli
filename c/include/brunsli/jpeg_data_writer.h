// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Functions for writing a JPEGData object into a jpeg byte stream.

#ifndef BRUNSLI_DEC_JPEG_DATA_WRITER_H_
#define BRUNSLI_DEC_JPEG_DATA_WRITER_H_

#include <brunsli/jpeg_data.h>
#include <brunsli/types.h>

namespace brunsli {

// Function pointer type used to write len bytes into buf. Returns the
// number of bytes written.
typedef size_t (*JPEGOutputHook)(void* data, const uint8_t* buf, size_t len);

// Output callback function with associated data.
struct JPEGOutput {
  JPEGOutput(JPEGOutputHook cb, void* data) : cb(cb), data(data) {}
  bool Write(const uint8_t* buf, size_t len) const {
    if (len == 0) return true;
    size_t bytes_written = cb(data, buf, len);
    return (bytes_written == len);
  }

 private:
  JPEGOutputHook cb;
  void* data;
};

bool WriteJpeg(const JPEGData& jpg, JPEGOutput out);

}  // namespace brunsli

#endif  // BRUNSLI_DEC_JPEG_DATA_WRITER_H_
