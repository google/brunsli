// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Library of cumulative distribution functions that can be used for arithmetic
// coding. For the semantics of these classes, see the comment on EncodeSymbol()
// in //third_party/brunsli/enc/arith_encode.h

#ifndef BRUNSLI_COMMON_DISTRIBUTIONS_H_
#define BRUNSLI_COMMON_DISTRIBUTIONS_H_

#include <brunsli/types.h>

namespace brunsli {

// An adaptive binary distribution with 8-bit precision.
struct Prob {
  static const int kInitProb = 134;
  static const int kInitProbCount = 3;
  Prob()
      : prob8(kInitProb),
        total(kInitProbCount),
        count(kInitProb * kInitProbCount) {}

  void Init(int probability) {
    prob8 = probability;
    total = kInitProbCount;
    count = kInitProbCount * probability;
  }

  void Add(int val) {
    // Python >>> [0, 0, 0] + [((1 << 17))/x for x in range(3,281)]
    static const uint16_t divlut[281] = {
      0, 0, 0, 43690, 32768, 26214, 21845, 18724, 16384, 14563,
      13107, 11915, 10922, 10082, 9362, 8738, 8192, 7710, 7281, 6898,
      6553, 6241, 5957, 5698, 5461, 5242, 5041, 4854, 4681, 4519,
      4369, 4228, 4096, 3971, 3855, 3744, 3640, 3542, 3449, 3360,
      3276, 3196, 3120, 3048, 2978, 2912, 2849, 2788, 2730, 2674,
      2621, 2570, 2520, 2473, 2427, 2383, 2340, 2299, 2259, 2221,
      2184, 2148, 2114, 2080, 2048, 2016, 1985, 1956, 1927, 1899,
      1872, 1846, 1820, 1795, 1771, 1747, 1724, 1702, 1680, 1659,
      1638, 1618, 1598, 1579, 1560, 1542, 1524, 1506, 1489, 1472,
      1456, 1440, 1424, 1409, 1394, 1379, 1365, 1351, 1337, 1323,
      1310, 1297, 1285, 1272, 1260, 1248, 1236, 1224, 1213, 1202,
      1191, 1180, 1170, 1159, 1149, 1139, 1129, 1120, 1110, 1101,
      1092, 1083, 1074, 1065, 1057, 1048, 1040, 1032, 1024, 1016,
      1008, 1000, 992, 985, 978, 970, 963, 956, 949, 942,
      936, 929, 923, 916, 910, 903, 897, 891, 885, 879,
      873, 868, 862, 856, 851, 845, 840, 834, 829, 824,
      819, 814, 809, 804, 799, 794, 789, 784, 780, 775,
      771, 766, 762, 757, 753, 748, 744, 740, 736, 732,
      728, 724, 720, 716, 712, 708, 704, 700, 697, 693,
      689, 686, 682, 679, 675, 672, 668, 665, 661, 658,
      655, 652, 648, 645, 642, 639, 636, 633, 630, 627,
      624, 621, 618, 615, 612, 609, 606, 604, 601, 598,
      595, 593, 590, 587, 585, 582, 579, 577, 574, 572,
      569, 567, 564, 562, 560, 557, 555, 553, 550, 548,
      546, 543, 541, 539, 537, 534, 532, 530, 528, 526,
      524, 522, 520, 518, 516, 514, 512, 510, 508, 506,
      504, 502, 500, 498, 496, 494, 492, 490, 489, 487,
      485, 483, 481, 480, 478, 476, 474, 473, 471, 469,
      468
    };
    ++total;
    if (val == 0) {
      count += 256;
    } else {
      ++count;
    }
    prob8 = ((uint32_t)(divlut[total]) * count) >> 17;
    static const int forget_rate = 254;
    if (total == forget_rate) {
      count >>= 1;
      total = forget_rate >> 1;
    }
  }
  uint8_t get_proba() const { return prob8; }

 private:
  uint8_t prob8;
  uint8_t total;
  uint16_t count;
};

}  // namespace brunsli

#endif  // BRUNSLI_COMMON_DISTRIBUTIONS_H_
