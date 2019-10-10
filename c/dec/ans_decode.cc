// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./ans_decode.h"

#include <vector>

#include "./histogram_decode.h"

namespace brunsli {

namespace {

bool ANSBuildMapTable(const int counts[], int alphabet_size,
                      ANSSymbolInfo map[ANS_TAB_SIZE]) {
  int i;
  int pos = 0;
  for (i = 0; i < alphabet_size; ++i) {
    int j;
    for (j = 0; j < counts[i]; ++j, ++pos) {
      map[pos].symbol_ = i;
      map[pos].freq_ = counts[i];
      map[pos].offset_ = j;
    }
  }
  return (pos == ANS_TAB_SIZE);
}

}  // namespace

bool ANSDecodingData::ReadFromBitStream(size_t alphabet_size,
                                        BrunsliBitReader* br) {
  std::vector<int> counts(alphabet_size);
  return (ReadHistogram(ANS_LOG_TAB_SIZE, alphabet_size, &counts[0], br) &&
          ANSBuildMapTable(&counts[0], alphabet_size, map_));
}

}  // namespace brunsli
