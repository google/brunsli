// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_BRUNSLI_AUX_DATA_H_
#define BRUNSLI_DEC_BRUNSLI_AUX_DATA_H_

#include <brunsli/types.h>

namespace brunsli {

struct BrunsliAuxData {
  BrunsliAuxData() : num_brunsli_contexts(0),
                     num_brunsli_histograms(0) {}

  // These fields are filled in by BrunsliDecodeJpeg() in
  // BRUNSLI_READ_SIZES mode.
  int num_brunsli_contexts;
  int num_brunsli_histograms;
};

}  // namespace brunsli

#endif  // BRUNSLI_DEC_BRUNSLI_AUX_DATA_H_
