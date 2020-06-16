// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "../enc/huffman_tree.h"

#include <vector>

#include "gtest/gtest.h"
#include "../common/platform.h"
#include <brunsli/types.h>

namespace brunsli {

TEST(HuffmanTree, Simple) {
  std::vector<uint32_t> histogram;
  histogram.push_back(1);
  histogram.push_back(2);
  histogram.push_back(4);
  histogram.push_back(8);
  histogram.push_back(16);
  histogram.push_back(32);
  std::vector<uint8_t> depth(histogram.size());
  CreateHuffmanTree(&histogram[0], histogram.size(), 6, &depth[0]);
  // The static casts to int are only to help ASSERT to show the value
  // in integer, not as a char.
  ASSERT_EQ(1u, depth[5]);
  ASSERT_EQ(2u, depth[4]);
  ASSERT_EQ(3u, depth[3]);
  ASSERT_EQ(4u, depth[2]);
  ASSERT_EQ(5u, depth[1]);
  ASSERT_EQ(5u, depth[0]);
}

TEST(HuffmanTree, Limited) {
  std::vector<uint32_t> histogram;
  histogram.push_back(1);
  histogram.push_back(2);
  histogram.push_back(4);
  histogram.push_back(8);
  histogram.push_back(16);
  histogram.push_back(32);
  histogram.push_back(64);
  histogram.push_back(128);
  // Check that by default the bit lengths exceed the limit we use later.
  {
    std::vector<uint8_t> depth(histogram.size());
    CreateHuffmanTree(&histogram[0], histogram.size(), 16, &depth[0]);
    EXPECT_LT(5u, depth[1]);
    EXPECT_LT(5u, depth[0]);
  }
  // Check that a limited length is obeyed.
  {
    std::vector<uint8_t> depth(histogram.size());
    CreateHuffmanTree(&histogram[0], histogram.size(), 5, &depth[0]);
    EXPECT_EQ(5u, depth[1]);
    EXPECT_EQ(5u, depth[0]);
  }
}

TEST(HuffmanTree, SimpleLiteral) {
  std::vector<uint32_t> histogram;
  histogram.push_back(0);
  histogram.push_back(10);
  histogram.push_back(0);
  histogram.push_back(0);
  std::vector<uint8_t> depth(histogram.size());
  CreateHuffmanTree(&histogram[0], histogram.size(), 10, &depth[0]);
  EXPECT_EQ(1u, depth[1]);
}

TEST(HuffmanTree, StableBitDepth) {
  std::vector<uint32_t> histogram;
  histogram.push_back(1);
  histogram.push_back(1);
  histogram.push_back(1);
  histogram.push_back(1);
  histogram.push_back(1);
  std::vector<uint8_t> depth(histogram.size());
  CreateHuffmanTree(&histogram[0], histogram.size(), 10, &depth[0]);
  EXPECT_EQ(2, depth[0]);
  EXPECT_EQ(2, depth[1]);
  EXPECT_EQ(2, depth[2]);
  EXPECT_EQ(3, depth[3]);
  EXPECT_EQ(3, depth[4]);
}

TEST(HuffmanTree, ConvertBitDepthsToSymbols) {
  std::vector<uint8_t> depth;
  depth.push_back(1);
  depth.push_back(0);
  depth.push_back(2);
  depth.push_back(2);
  std::vector<uint16_t> bits(depth.size());
  ConvertBitDepthsToSymbols(&depth[0], depth.size(), bits.data());
  ASSERT_EQ(0, bits[0]);
  ASSERT_EQ(0, bits[1]);
  ASSERT_EQ(1, bits[2]);
  ASSERT_EQ(3, bits[3]);
}

TEST(HuffmanTree, ConvertBitDepthsToSymbols2) {
  std::vector<uint8_t> depth;

  depth.push_back(0);
  depth.push_back(0);
  depth.push_back(3);
  depth.push_back(3);
  depth.push_back(3);
  depth.push_back(1);

  std::vector<uint16_t> bits(depth.size());
  ConvertBitDepthsToSymbols(&depth[0], depth.size(), bits.data());

  ASSERT_EQ(0, bits[0]);
  ASSERT_EQ(0, bits[1]);
  ASSERT_EQ(1, bits[2]);
  ASSERT_EQ(5, bits[3]);
  ASSERT_EQ(3, bits[4]);
  ASSERT_EQ(0, bits[5]);
}

TEST(HuffmanTree, WriteHuffmanTree) {
  std::vector<uint8_t> depth;
  depth.push_back(1);
  depth.push_back(0);  // This streak of 0s should be represented by a single
  depth.push_back(0);  // marker.
  depth.push_back(2);
  std::vector<uint8_t> tree(depth.size());
  std::vector<uint8_t> extra_bits_data(depth.size());
  size_t tree_size = 0;
  WriteHuffmanTree(depth.data(), 4, &tree_size, &tree[0], &extra_bits_data[0]);
  ASSERT_EQ(4u, tree_size);
}

TEST(HuffmanTree, WriteHuffmanTreeSparse) {
  std::vector<uint8_t> depth(286 + 30);
  // Some literals and lengths:
  depth[0] = 1;
  depth[256] = 2;
  depth[257] = 3;
  depth[265] = 3;

  // Some backward references:
  depth[286 + 0] = 1;
  depth[286 + 3] = 1;

  std::vector<uint8_t> tree(depth.size());
  std::vector<uint8_t> extra_bits_data(depth.size());
  size_t tree_size = 0;
  WriteHuffmanTree(depth.data(), 286 + 30, &tree_size, &tree[0],
                   &extra_bits_data[0]);
  ASSERT_EQ(14u, tree_size);

  ASSERT_EQ(1u, tree[0]);
  ASSERT_EQ(17u, tree[1]);
  ASSERT_EQ(17u, tree[2]);
  ASSERT_EQ(17u, tree[3]);
  EXPECT_EQ(2u, extra_bits_data[1]);  // 2 + 1 = 3
  EXPECT_EQ(6u, extra_bits_data[2]);  // 3 * 8 + 6 + 1 = 31
  EXPECT_EQ(4u, extra_bits_data[3]);  // 31 * 8 + 4 + 3 = 255

  ASSERT_EQ(2u, tree[4]);
  ASSERT_EQ(3u, tree[5]);
  ASSERT_EQ(17u, tree[6]);
  ASSERT_EQ(7u, extra_bits_data[6] + 3u);
  ASSERT_EQ(3u, tree[7]);
  ASSERT_EQ(17u, tree[8]);
  ASSERT_EQ(17u, tree[9]);
  EXPECT_EQ(1u, extra_bits_data[8]);  // 1 + 1 = 2
  EXPECT_EQ(1u, extra_bits_data[9]);  // 2 * 8 + 1 + 3 = 20
  ASSERT_EQ(1u, tree[10]);
  ASSERT_EQ(0u, tree[11]);
  ASSERT_EQ(0u, tree[12]);
  ASSERT_EQ(1u, tree[13]);

  // Fine so far, but let's take a Huffman tree of the Huffman tree.

  // RFC 1951. 3.2.7. Compression with dynamic Huffman codes (BTYPE=10)
  // Calculate the statistics of the Huffman tree in deflate-representation.
  const size_t kCodeLengthCodes = 19;
  std::vector<uint32_t> huffman_tree_histogram(kCodeLengthCodes);
  for (size_t i = 0; i < tree.size(); ++i) {
    ++huffman_tree_histogram[tree[i]];
    BRUNSLI_LOG_DEBUG() << "huffman_tree: " << i << " "
                        << static_cast<int>(tree[i]) << BRUNSLI_ENDL();
  }
  std::vector<uint8_t> code_length_bitdepth(huffman_tree_histogram.size());
  std::vector<uint16_t> code_length_bitdepth_symbols(
      huffman_tree_histogram.size());
  CreateHuffmanTree(&huffman_tree_histogram[0],
                    huffman_tree_histogram.size(),
                    7, &code_length_bitdepth[0]);
  for (size_t i = 0; i < code_length_bitdepth.size(); ++i) {
    BRUNSLI_LOG_DEBUG() << "code_length_bitdepth " << i << " "
                        << static_cast<int>(code_length_bitdepth[i])
                        << BRUNSLI_ENDL();
  }

  for (size_t i = 0; i < code_length_bitdepth.size(); ++i) {
    ASSERT_GT(8, code_length_bitdepth[i]);
  }
  ConvertBitDepthsToSymbols(&code_length_bitdepth[0],
                            code_length_bitdepth.size(),
                            code_length_bitdepth_symbols.data());
  for (size_t i = 0; i < code_length_bitdepth_symbols.size(); ++i) {
    BRUNSLI_LOG_DEBUG() << "code_length_symbols " << i << " "
                        << code_length_bitdepth_symbols[i] << BRUNSLI_ENDL();
  }
}

TEST(HuffmanTree, WriteHuffmanTreeManyZeros) {
  std::vector<uint8_t> depth;
  depth.push_back(1);
  for (size_t i = 0; i < 148; ++i) {
    depth.push_back(0);
  }
  depth.push_back(2);
  std::vector<uint8_t> tree(depth.size());
  std::vector<uint8_t> extra_bits_data(depth.size());
  size_t tree_size = 0;
  WriteHuffmanTree(depth.data(), 150, &tree_size, &tree[0],
                   &extra_bits_data[0]);
  ASSERT_EQ(17u, tree[1]);
  ASSERT_EQ(17u, tree[2]);
  ASSERT_EQ(17u, tree[3]);
  ASSERT_EQ(5u, tree_size);
}

TEST(HuffmanTree, WriteHuffmanTreeShortStipeOfNonZeros) {
  std::vector<uint8_t> depth;
  depth.push_back(9);
  for (size_t i = 0; i < 5; ++i) {
    depth.push_back(8);
  }
  depth.push_back(10);
  std::vector<uint8_t> tree(depth.size());
  std::vector<uint8_t> extra_bits_data(depth.size());
  size_t tree_size = 0;
  WriteHuffmanTree(depth.data(), 7, &tree_size, &tree[0], &extra_bits_data[0]);
  EXPECT_GE(7u, tree_size);
  EXPECT_EQ(9, tree[0]);
  EXPECT_EQ(8, tree[1]);
  EXPECT_EQ(8, tree[2]);
  EXPECT_EQ(8, tree[3]);
  EXPECT_EQ(8, tree[4]);
  EXPECT_EQ(8, tree[5]);
  EXPECT_EQ(10, tree[6]);
}

TEST(HuffmanTree, WriteHuffmanTreeManyNonZeros) {
  std::vector<uint8_t> depth;
  depth.push_back(9);
  for (size_t i = 0; i < 200; ++i) {
    depth.push_back(8);
  }
  depth.push_back(10);
  std::vector<uint8_t> tree(depth.size());
  std::vector<uint8_t> extra_bits_data(depth.size());
  size_t tree_size = 0;
  WriteHuffmanTree(depth.data(), 202, &tree_size, &tree[0],
                   &extra_bits_data[0]);
  ASSERT_GE(60u, tree_size);  // Is now 60, the actual value is repeated
  // too many times. Ideally, this should be 32. This leads to bloat of
  // 20 bytes per block or so in some rare cases.
  ASSERT_EQ(9, tree[0]);
  ASSERT_EQ(8, tree[1]);
  ASSERT_EQ(16, tree[2]);
  ASSERT_EQ(10, tree[tree_size - 1]);
}

}  // namespace brunsli
