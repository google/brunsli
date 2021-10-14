// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Functions for writing a JPEGData object into a brunsli byte stream and
// converting a jpeg byte stream to a brunsli byte stream.

#ifndef BRUNSLI_ENC_BRUNSLI_ENCODE_H_
#define BRUNSLI_ENC_BRUNSLI_ENCODE_H_

#include <brunsli/jpeg_data.h>
#include <brunsli/types.h>

namespace brunsli {

// Returns an upper bound on the size of the buffer needed to encode the given
// jpg data in brunsli format.
size_t GetMaximumBrunsliEncodedSize(const JPEGData& jpg);

#if defined(BRUNSLI_EXTRA_API)
// Returns the estimated peak memory usage (in bytes) of BrunsliEncodeJpeg().
size_t EstimateBrunsliEncodePeakMemoryUsage(size_t jpg_size,
                                            const JPEGData& jpg);
#endif  // defined(BRUNSLI_EXTRA_API)

// Encodes the given jpg to the buffer data[0 ... *len) in brunsli format and
// sets *len to the encoded size. Returns false on buffer overflow or invalid
// jpg data.
bool BrunsliEncodeJpeg(const JPEGData& jpg, uint8_t* data, size_t* len);

// Return the storage size needed to store raw jpg data in bypass mode.
size_t GetBrunsliBypassSize(size_t jpg_size);

// Bypass mode: store the JPEG data directly into brunsli format. *len should
// contain the maximum storage space available. Upon return, *len is updated to
// the actual size used in 'data'.
bool BrunsliEncodeJpegBypass(const uint8_t* jpg_data, size_t jpg_data_len,
                             uint8_t* data, size_t* len);

}  // namespace brunsli

#endif  // BRUNSLI_ENC_BRUNSLI_ENCODE_H_
