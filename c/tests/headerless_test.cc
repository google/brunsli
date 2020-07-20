// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <string>

#include "gtest/gtest.h"
#include "../common/constants.h"
#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>
#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>
#include "../dec/state.h"
#include "./test_utils.h"

namespace brunsli {

using ::brunsli::internal::dec::Stage;
using ::brunsli::internal::dec::State;

TEST(DecodeTest, TestHeaderless) {
  JPEGData jpg;
  // Vanilla context modelling.
  jpg.version = 0;
  jpg.width = 16;
  jpg.height = 16;
  jpg.components.resize(3);
  jpg.max_h_samp_factor = 1;
  jpg.max_v_samp_factor = 1;
  jpg.MCU_rows = 2;
  jpg.MCU_cols = 2;
  for (size_t c = 0; c < 3; ++c) {
    jpg.components[c].v_samp_factor = 1;
    jpg.components[c].h_samp_factor = 1;
  }
  ASSERT_TRUE(::brunsli::internal::dec::UpdateSubsamplingDerivatives(&jpg));

  std::vector<uint8_t> src = GetSmallBrunsliFile();
  size_t src_offset = kSmallBrunsliSignatuteSize + kSmallBrunsliHeaderSize;

  State state;
  // Vanilla context modelling.
  state.use_legacy_context_model = true;
  state.stage = Stage::SECTION;
  state.tags_met = (1 << kBrunsliSignatureTag) | (1 << kBrunsliHeaderTag) |
                   (1 << kBrunsliOriginalJpgTag);
  state.data = src.data() + src_offset;
  state.len = src.size() - src_offset;

  ::brunsli::internal::dec::PrepareMeta(&jpg, &state);

  ASSERT_EQ(BRUNSLI_OK, ::brunsli::internal::dec::ProcessJpeg(&state, &jpg));
}

}  // namespace brunsli
