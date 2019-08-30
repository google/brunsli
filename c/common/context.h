// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_COMMON_CONTEXT_H_
#define BRUNSLI_COMMON_CONTEXT_H_

#include <vector>

#include "./distributions.h"
#include <brunsli/jpeg_data.h>
#include "./platform.h"
#include <brunsli/types.h>

namespace brunsli {

static const int kMaxAverageContext = 8;
static const int kNumAvrgContexts = kMaxAverageContext + 1;
static const int kInterSymbolContexts = 63;

static const uint8_t kNonzeroBuckets[64] = {
   0,   1,   2,   3,   4,   4,   5,   5,   5,   6,   6,   6,   6,   7,   7,   7,
   7,   7,   7,   7,   7,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,   8,
   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,  10,  10,  10,
  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,
};
// kNonzeroBuckets[i] < kNumNonzeroBuckets
static const uint8_t kNumNonzeroBuckets = 11;

static const int kNumSchemes = 7;

static const uint8_t kFreqContext[kNumSchemes][64] = {
  { 0, },

  { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  0,  0, },

  { 0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  3,
    3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,
    3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  1,  1,  1, },

  { 0,  1,  1,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,
    5,  5,  5,  5,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,  6,
    6,  6,  6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,  7,  7,
    7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  2,  2,  2, },

  { 0,  1,  2,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8,  8,  8,
    9,  9,  9,  9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12,
   13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14,
   15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, },

  { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
   16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23,
   24, 24, 24, 24, 25, 25, 25, 25, 26, 26, 26, 26, 27, 27, 27, 27,
   28, 28, 28, 28, 29, 29, 29, 29, 30, 30, 30, 30, 31, 31, 31, 31, },

  { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
   16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
   32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
   48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, },
};

static const uint16_t kNumNonzeroContext[kNumSchemes][64] = {
 { 0,   1,   1,   2,   2,   2,   3,   3,   3,   3,   4,   4,   4,   4,   4,   4,
   5,   5,   5,   5,   5,   5,   5,   5,   6,   6,   6,   6,   6,   6,   6,   6,
   6,   6,   6,   6,   6,   6,   6,   6,   7,   7,   7,   7,   7,   7,   7,   7,
   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7,   7
 },
 { 0,   2,   2,   4,   4,   4,   6,   6,   6,   6,   8,   8,  8,   8,   8,   8,
  10,  10,  10,  10,  10,  10,  10,  10,  12,  12,  12,  12,  12,  12,  12,  12,
  12,  12,  12,  12,  12,  12,  12,  12,  14,  14,  14,  14,  14,  14,  14,  14,
  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14,  14
 },
 { 0,   4,   4,   8,   8,   8,  12,  12,  12,  12,  16,  16,  16,  16,  16,  16,
  20,  20,  20,  20,  20,  20,  20,  20,  24,  24,  24,  24,  24,  24,  24,  24,
  24,  24,  24,  24,  24,  24,  24,  24,  28,  28,  28,  28,  28,  28,  28,  28,
  28,  28,  28,  28,  28,  28,  28,  28,  28,  28,  28,  28,  28,  28,  28,  28
 },
 { 0,   8,   8,  16,  16,  16,  24,  24,  24,  24,  32,  32,  32,  32,  32,  32,
  40,  40,  40,  40,  40,  40,  40,  40,  48,  48,  48,  48,  48,  48,  48,  48,
  48,  48,  48,  48,  48,  48,  48,  48,  55,  55,  55,  55,  55,  55,  55,  55,
  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55,  55
 },
 { 0,  16,  16,  32,  32,  32,  48,  48,  48,  48,  64,  64,  64,  64,  64,  64,
  80,  80,  80,  80,  80,  80,  80,  80,  95,  95,  95,  95,  95,  95,  95,  95,
  95,  95,  95,  95,  95,  95,  95,  95, 109, 109, 109, 109, 109, 109, 109, 109,
 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109
 },
 { 0,  32,  32,  64,  64,  64,  96,  96,  96,  96, 127, 127, 127, 127, 127, 127,
 157, 157, 157, 157, 157, 157, 157, 157, 185, 185, 185, 185, 185, 185, 185, 185,
 185, 185, 185, 185, 185, 185, 185, 185, 211, 211, 211, 211, 211, 211, 211, 211,
 211, 211, 211, 211, 211, 211, 211, 211, 211, 211, 211, 211, 211, 211, 211, 211
 },
 { 0,  64,  64, 127, 127, 127, 188, 188, 188, 188, 246, 246, 246, 246, 246, 246,
 300, 300, 300, 300, 300, 300, 300, 300, 348, 348, 348, 348, 348, 348, 348, 348,
 348, 348, 348, 348, 348, 348, 348, 348, 388, 388, 388, 388, 388, 388, 388, 388,
 388, 388, 388, 388, 388, 388, 388, 388, 388, 388, 388, 388, 388, 388, 388, 388
 }
};

static const uint16_t kNumNonzeroContextSkip[kNumSchemes] = {
  8, 15, 31, 61, 120, 231, 412
};

inline int ZeroDensityContext(int nonzeros_left, int k, int bits) {
  return kNumNonzeroContext[bits][nonzeros_left] + kFreqContext[bits][k];
}

// Returns the context for the absolute value of the prediction error of
// the next DC coefficient in column x, using the one row size ringbuffer of
// previous absolute prediction errors in vals.
inline int WeightedAverageContextDC(const int* vals, int x) {
  // Since vals is a ringbuffer, vals[x] and vals[x + 1] refer to the
  // previous row.
  int sum = 1 + vals[x - 2] + vals[x - 1] + vals[x] + vals[x + 1];
  if ((sum >> kMaxAverageContext) != 0) {
    return kMaxAverageContext;
  }
  return Log2FloorNonZero(sum);
}

inline int WeightedAverageContext(const int* vals, int prev_row_delta) {
  int sum = 4 +
      vals[0] +
      (vals[-kDCTBlockSize] + vals[prev_row_delta]) * 2 +
      vals[-2 * kDCTBlockSize] +
      vals[prev_row_delta - kDCTBlockSize] +
      vals[prev_row_delta + kDCTBlockSize];
  if ((sum >> (kMaxAverageContext + 2)) != 0) {
    return kMaxAverageContext;
  }
  return Log2FloorNonZero(sum) - 2;
}

static const int kACPredictPrecisionBits = 13;
static const int kACPredictPrecision = 1 << kACPredictPrecisionBits;

void ComputeACPredictMultipliers(const int* quant,
                                 int* mult_row,
                                 int* mult_col);

// Returns a context value from the AC prediction in the range
// [-kMaxAverageContext, kMaxAverageContext]
inline int ACPredictContext(int64_t p_orig) {
  uint64_t p;
  static const int64_t kMaxPred = (1 << (kMaxAverageContext + 1)) - 1;
  if (p_orig >= 0) {
    p = static_cast<uint64_t>(p_orig) << 1;
    if (p > kMaxPred) {
      return kMaxAverageContext;
    }
    return Log2FloorNonZero(++p);
  }
  p = static_cast<uint64_t>(-p_orig) << 1;
  if (p > kMaxPred) {
    return -kMaxAverageContext;
  } else {
    return -Log2FloorNonZero(++p);
  }
}

inline int ACPredictContextCol(const coeff_t* prev,
                               const coeff_t* cur,
                               const int* mult) {
  coeff_t terms[8];
  terms[0] = 0;
  terms[1] = cur[1] + prev[1];
  terms[2] = cur[2] - prev[2];
  terms[3] = cur[3] + prev[3];
  terms[4] = cur[4] - prev[4];
  terms[5] = cur[5] + prev[5];
  terms[6] = cur[6] - prev[6];
  terms[7] = cur[7] + prev[7];
  int64_t delta = terms[0] * static_cast<int64_t>(mult[0]) +
                  terms[1] * static_cast<int64_t>(mult[1]) +
                  terms[2] * static_cast<int64_t>(mult[2]) +
                  terms[3] * static_cast<int64_t>(mult[3]) +
                  terms[4] * static_cast<int64_t>(mult[4]) +
                  terms[5] * static_cast<int64_t>(mult[5]) +
                  terms[6] * static_cast<int64_t>(mult[6]) +
                  terms[7] * static_cast<int64_t>(mult[7]);
  return ACPredictContext(prev[0] - delta / kACPredictPrecision);
}

inline int ACPredictContextRow(const coeff_t* prev,
                               const coeff_t* cur,
                               const int* mult) {
  coeff_t terms[8];
  terms[0] = 0;
  terms[1] = cur[8] + prev[8];
  terms[2] = cur[16] - prev[16];
  terms[3] = cur[24] + prev[24];
  terms[4] = cur[32] - prev[32];
  terms[5] = cur[40] + prev[40];
  terms[6] = cur[48] - prev[48];
  terms[7] = cur[56] + prev[56];
  int64_t delta = terms[0] * static_cast<int64_t>(mult[0]) +
                  terms[1] * static_cast<int64_t>(mult[1]) +
                  terms[2] * static_cast<int64_t>(mult[2]) +
                  terms[3] * static_cast<int64_t>(mult[3]) +
                  terms[4] * static_cast<int64_t>(mult[4]) +
                  terms[5] * static_cast<int64_t>(mult[5]) +
                  terms[6] * static_cast<int64_t>(mult[6]) +
                  terms[7] * static_cast<int64_t>(mult[7]);
  return ACPredictContext(prev[0] - delta / kACPredictPrecision);
}

inline int NumNonzerosContext(const int* prev, int x, int y) {
  return (y == 0 ? (prev[x - 1] >> 1) :
          x == 0 ? (prev[0] >> 1) :
          (prev[x - 1] + prev[x] + 1) >> 2);
}

// Context for the emptyness of a block is the number of non-empty blocks in the
// previous and up neighborhood (blocks beyond the border are assumed
// non-empty).
static const int kNumIsEmptyBlockContexts = 3;
inline int IsEmptyBlockContext(const int* prev, int x) {
  return prev[x - 1] + prev[x];
}

// Holds all encoding/decoding state for an image component that is needed to
// switch between components during interleaved encoding/decoding.
struct ComponentStateDC {
  ComponentStateDC() :
      width(0),
      is_empty_block_prob(kNumIsEmptyBlockContexts),
      sign_prob(9),
      first_extra_bit_prob(10) {
    InitAll();
  }

  void SetWidth(int w) {
    width = w;
    prev_is_nonempty.resize(w + 1, 1);
    prev_abs_coeff.resize(w + 3);
    prev_sign.resize(w + 1);
  }

  int width;
  Prob is_zero_prob;
  std::vector<Prob> is_empty_block_prob;
  std::vector<Prob> sign_prob;
  std::vector<Prob> first_extra_bit_prob;
  std::vector<int> prev_is_nonempty;
  std::vector<int> prev_abs_coeff;
  std::vector<int> prev_sign;

 protected:
  void InitAll();
};

struct ComponentState {
  ComponentState() :
      width(0),
      is_zero_prob(kNumNonzeroBuckets * kDCTBlockSize),
      sign_prob((2 * kMaxAverageContext + 1) * kDCTBlockSize),
      first_extra_bit_prob(10 * kDCTBlockSize) {
    InitAll();
  }

  void SetWidth(int w) {
    width = w;
    prev_is_nonempty.resize(w + 1, 1);
    prev_num_nonzeros.resize(w + 1);
    prev_abs_coeff.resize(kDCTBlockSize * 2 * (w + 3));
    prev_sign.resize(kDCTBlockSize * (w + 1));
  }

  // Returns the size of the object after constructor and SetWidth(w).
  // Used in estimating peak heap memory usage of the brunsli codec.
  static size_t SizeInBytes(int w) {
    return (4 + (10 + 3 * w) * kDCTBlockSize + 2 * w) * sizeof(int) +
        ((kNumNonzeroBuckets + 2 * kMaxAverageContext + 11) * kDCTBlockSize +
         32 * kInterSymbolContexts) * sizeof(Prob);
  }

  int width;
  int context_offset;
  int order[kDCTBlockSize];
  int mult_row[kDCTBlockSize];
  // mult_col is transposed for more effective ACPredictContextRow execution.
  int mult_col[kDCTBlockSize];
  std::vector<Prob> is_zero_prob;
  std::vector<Prob> sign_prob;
  Prob num_nonzero_prob[32][kInterSymbolContexts];
  std::vector<Prob> first_extra_bit_prob;
  std::vector<int> prev_is_nonempty;
  std::vector<int> prev_num_nonzeros;
  std::vector<int> prev_abs_coeff;
  std::vector<int> prev_sign;

 protected:
  void InitAll();
};

}  // namespace brunsli

#endif  // BRUNSLI_COMMON_CONTEXT_H_
