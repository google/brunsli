// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Functions for reading a brunsli byte stream into a JPEGData object and
// converting a brunsli byte stream to a jpeg byte stream.

#ifndef BRUNSLI_DEC_BRUNSLI_DECODE_H_
#define BRUNSLI_DEC_BRUNSLI_DECODE_H_

#include <memory>
#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>

namespace brunsli {

namespace internal {
namespace dec {
class State;
}  // namespace dec
}  // namespace internal

// Parses the brunsli byte stream contained in data[0 ... len) and fills in *jpg
// with the parsed information.
// The *jpg object is valid only as long as the input data is valid.
// Returns BRUNSLI_OK, unless the data is not valid brunsli byte stream, or is
// truncated.
BrunsliStatus BrunsliDecodeJpeg(const uint8_t* data, size_t len, JPEGData* jpg);

/* Check if data looks like Brunsli stream.
 * Currently, only 6 byte signature is compared
 * (i.e. if |len| < 6, result is always "false").
 */
bool IsBrunsli(const uint8_t* data, size_t len);

// Returns the estimated peak memory usage (in bytes) of the BrunsliDecodeJpeg
// function. If parsing is failed, then result is 0.
size_t BrunsliEstimateDecoderPeakMemoryUsage(const uint8_t* data, size_t len);

class BrunsliDecoder {
 public:
  BrunsliDecoder();
  ~BrunsliDecoder();

  enum Status {
    NEEDS_MORE_INPUT,
    NEEDS_MORE_OUTPUT,
    ERROR,
    DONE,
  };

  Status Decode(size_t* available_in, const uint8_t** next_in,
                size_t* available_out, uint8_t** next_out);

 private:
  std::unique_ptr<JPEGData> jpg_;
  std::unique_ptr<::brunsli::internal::dec::State> state_;
};

}  // namespace brunsli

#endif  // BRUNSLI_DEC_BRUNSLI_DECODE_H_
