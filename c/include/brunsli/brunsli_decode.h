// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Functions for reading a brunsli byte stream into a JPEGData object and
// converting a brunsli byte stream to a jpeg byte stream.

#ifndef BRUNSLI_DEC_BRUNSLI_DECODE_H_
#define BRUNSLI_DEC_BRUNSLI_DECODE_H_

#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>

namespace brunsli {

typedef enum {
  BRUNSLI_READ_HEADER,   // only signature and header
  BRUNSLI_READ_SIZES,    // header and histogram sections
  BRUNSLI_READ_ALL,      // everything
} BrunsliReadMode;

struct BrunsliAuxData;

// Parses the brunsli byte stream contained in data[0 ... len) and fills in *jpg
// with the parsed information.
// The *jpg object is valid only as long as the input data is valid.
// Returns BRUNSLI_OK, unless the data is not valid brunsli byte stream, or is
// truncated.
BrunsliStatus BrunsliDecodeJpeg(const uint8_t* data, const size_t len,
                                BrunsliReadMode mode,
                                JPEGData* jpg,
                                BrunsliAuxData* aux);

}  // namespace brunsli

#endif  // BRUNSLI_DEC_BRUNSLI_DECODE_H_
