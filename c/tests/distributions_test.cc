// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "../common/distributions.h"

#include <vector>

#include "gtest/gtest.h"
#include <brunsli/types.h>
#include "../dec/arith_decode.h"
#include "../dec/brunsli_input.h"

namespace brunsli {

namespace {

// Stripped version of DataStream from brunsli_encoder.cc
class DataStream {
 public:
  DataStream() : low_(0), high_(~0) {}

  void Flush() {
    words.push_back(high_ >> 16);
    words.push_back(high_ & 0xffff);
  }

  void AddBit(Prob* p, int bit) {
    const uint8_t prob = p->get_proba();
    p->Add(bit);
    const uint32_t diff = high_ - low_;
    const uint32_t split = low_ + (((uint64_t)diff * prob) >> 8);
    if (bit) {
      low_ = split + 1;
    } else {
      high_ = split;
    }
    if (((low_ ^ high_) >> 16) == 0) {
      words.push_back(high_ >> 16);
      low_ <<= 16;
      high_ <<= 16;
      high_ |= 0xffff;
    }
  }

  uint32_t low_;
  uint32_t high_;
  std::vector<uint16_t> words;
};

void Roundtrip(const std::vector<int>& text) {
  DataStream ds;
  std::vector<Prob> p_out(4);
  for (size_t i = 0; i < text.size(); ++i) {
    ds.AddBit(&p_out[i & 3], text[i]);
  }
  ds.Flush();

  WordSource in(reinterpret_cast<const uint8_t*>(ds.words.data()),
                ds.words.size() * sizeof(uint16_t), true);
  BinaryArithmeticDecoder bad;
  bad.Init(&in);
  std::vector<Prob> p_in(4);
  for (size_t i = 0; i < text.size(); ++i) {
    int bit = bad.ReadBit(p_in[i & 3].get_proba(), &in);
    p_in[i & 3].Add(bit);
    ASSERT_EQ(text[i], bit);
  }
}

}  // namespace

TEST(Distributions, Worst) {
  DataStream ds;
  for (size_t i = 0; i < 640; ++i) {
    Prob p;
    p.Init(0);
    ds.AddBit(&p, 0);
  }
  ds.Flush();
  // This means that every AddBit emits the whole word (16 bits).
  // That is far from ideal (see "Bad" test case).
  ASSERT_EQ(642, ds.words.size());
}

// TODO(eustas): explore if this could happen for valid JPEG input; also, if
//               this could be triggered by Brunsli input (CPU-zip-bomb)
TEST(Distributions, Perfect) {
  DataStream ds;
  {
    Prob p;
    p.Init(255);
    ds.AddBit(&p, 1);
  }
  // This is a backflip of "Worst" case - 4 billions bits could be stored
  // in 2 words...
  for (size_t i = 0; i < 0xFF0000 - 1; ++i) {
    Prob p;
    p.Init(0);
    ds.AddBit(&p, 1);
  }
  ASSERT_EQ(0, ds.words.size());
  {
    Prob p;
    p.Init(0);
    ds.AddBit(&p, 1);
  }
  ASSERT_EQ(1, ds.words.size());
}

TEST(Distributions, Bad) {
  DataStream ds;
  for (size_t i = 0; i < 640; ++i) {
    Prob p;
    p.Init(255);
    ds.AddBit(&p, 1);
  }
  ds.Flush();
  // Unlike "Worst" just one byte (amortized) is used per AddBit.
  ASSERT_EQ(322, ds.words.size());
}

TEST(Distributions, Good) {
  DataStream ds;
  for (size_t i = 0; i < 2833; ++i) {
    Prob p;
    p.Init(255);
    ds.AddBit(&p, 0);
  }
  ASSERT_EQ(0, ds.words.size());
  {
    Prob p;
    p.Init(255);
    ds.AddBit(&p, 0);
  }
  ASSERT_EQ(1, ds.words.size());
}

TEST(Distributions, TowardsZero) {
  Prob p;
  p.Init(255);
  for (size_t i = 0; i < 10000; ++i) {
    p.Add(0);
    ASSERT_EQ(255, p.get_proba()) << i;
  }
}

TEST(Distributions, TowardsOne) {
  Prob p;
  p.Init(0);
  for (size_t i = 0; i < 10000; ++i) {
    p.Add(1);
    ASSERT_EQ(0, p.get_proba()) << i;
  }
}

TEST(Distributions, Extremes) {
  std::vector<int> text;
  for (size_t i = 0; i < 50000; ++i) {
    text.push_back(i >= 25000);
  }
  for (size_t i = 0; i < 4; ++i) {
    text[24000 + i] ^= 1;
    text[49000 + i] ^= 1;
  }
  for (size_t i = 0; i < 16; ++i) {
    text[24531 + i] ^= 1;
    text[49531 + i] ^= 1;
  }
  Roundtrip(text);
}

}  // namespace brunsli
