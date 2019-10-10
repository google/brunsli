// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./huffman_decode.h"

#include <cstring>  /* for memset */
#include <vector>

#include <brunsli/types.h>
#include "./bit_reader.h"
#include "./huffman_table.h"

namespace brunsli {

static const int kCodeLengthCodes = 18;
static const uint8_t kCodeLengthCodeOrder[kCodeLengthCodes] = {
    1, 2, 3, 4, 0, 5, 17, 6, 16, 7, 8, 9, 10, 11, 12, 13, 14, 15,
};
static const uint8_t kDefaultCodeLength = 8;
static const uint8_t kCodeLengthRepeatCode = 16;

int ReadHuffmanCodeLengths(const uint8_t* code_length_code_lengths,
                           int num_symbols, uint8_t* code_lengths,
                           BrunsliBitReader* br) {
  int symbol = 0;
  uint8_t prev_code_len = kDefaultCodeLength;
  int repeat = 0;
  uint8_t repeat_code_len = 0;
  int space = 32768;
  HuffmanCode table[32];

  uint16_t counts[16] = {0};
  for (int i = 0; i < kCodeLengthCodes; ++i) {
    ++counts[code_length_code_lengths[i]];
  }
  if (!BuildHuffmanTable(table, 5, code_length_code_lengths, kCodeLengthCodes,
                         &counts[0])) {
    return 0;
  }

  while (symbol < num_symbols && space > 0) {
    const HuffmanCode* p = table;
    uint8_t code_len;
    if (!BrunsliBitReaderReadMoreInput(br)) {
      return 0;
    }
    BrunsliBitReaderFillWindow(br, 5);
    p += (br->val_ >> br->bit_pos_) & 31;
    br->bit_pos_ += p->bits;
    code_len = (uint8_t)p->value;
    if (code_len < kCodeLengthRepeatCode) {
      repeat = 0;
      code_lengths[symbol++] = code_len;
      if (code_len != 0) {
        prev_code_len = code_len;
        space -= 32768u >> code_len;
      }
    } else {
      const int extra_bits = code_len - 14;
      int old_repeat;
      int repeat_delta;
      uint8_t new_len = 0;
      if (code_len == kCodeLengthRepeatCode) {
        new_len = prev_code_len;
      }
      if (repeat_code_len != new_len) {
        repeat = 0;
        repeat_code_len = new_len;
      }
      old_repeat = repeat;
      if (repeat > 0) {
        repeat -= 2;
        repeat <<= extra_bits;
      }
      repeat += (int)BrunsliBitReaderReadBits(br, extra_bits) + 3;
      repeat_delta = repeat - old_repeat;
      if (symbol + repeat_delta > num_symbols) {
        return 0;
      }
      memset(&code_lengths[symbol], repeat_code_len, (size_t)repeat_delta);
      symbol += repeat_delta;
      if (repeat_code_len != 0) {
        space -= repeat_delta << (15 - repeat_code_len);
      }
    }
  }
  if (space != 0) {
    return 0;
  }
  memset(&code_lengths[symbol], 0, (size_t)(num_symbols - symbol));
  return 1;
}

bool HuffmanDecodingData::ReadFromBitStream(int alphabet_size,
                                            BrunsliBitReader* br) {
  int ok = 1;
  int table_size = 0;
  int simple_code_or_skip;

  std::vector<uint8_t> code_lengths(alphabet_size, 0);
  if (!BrunsliBitReaderReadMoreInput(br)) {
    return false;
  }
  /* simple_code_or_skip is used as follows:
     1 for simple code;
     0 for no skipping, 2 skips 2 code lengths, 3 skips 3 code lengths */
  simple_code_or_skip = (int)BrunsliBitReaderReadBits(br, 2);
  if (simple_code_or_skip == 1) {
    /* Read symbols, codes & code lengths directly. */
    int i;
    int max_bits_counter = alphabet_size - 1;
    int max_bits = 0;
    int symbols[4] = {0};
    const int num_symbols = (int)BrunsliBitReaderReadBits(br, 2) + 1;
    while (max_bits_counter) {
      max_bits_counter >>= 1;
      ++max_bits;
    }
    for (i = 0; i < num_symbols; ++i) {
      symbols[i] = (int)BrunsliBitReaderReadBits(br, max_bits) % alphabet_size;
      code_lengths[symbols[i]] = 2;
    }
    code_lengths[symbols[0]] = 1;
    switch (num_symbols) {
      case 1:
        break;
      case 3:
        ok = ((symbols[0] != symbols[1]) && (symbols[0] != symbols[2]) &&
              (symbols[1] != symbols[2]));
        break;
      case 2:
        ok = (symbols[0] != symbols[1]);
        code_lengths[symbols[1]] = 1;
        break;
      case 4:
        ok = ((symbols[0] != symbols[1]) && (symbols[0] != symbols[2]) &&
              (symbols[0] != symbols[3]) && (symbols[1] != symbols[2]) &&
              (symbols[1] != symbols[3]) && (symbols[2] != symbols[3]));
        if (BrunsliBitReaderReadBits(br, 1)) {
          code_lengths[symbols[2]] = 3;
          code_lengths[symbols[3]] = 3;
        } else {
          code_lengths[symbols[0]] = 2;
        }
        break;
      default:
        // Unreachable.
        return false;
    }
  } else { /* Decode Huffman-coded code lengths. */
    int i;
    uint8_t code_length_code_lengths[kCodeLengthCodes] = {0};
    int space = 32;
    int num_codes = 0;
    /* Static Huffman code for the code length code lengths */
    static const HuffmanCode huff[16] = {
        {2, 0}, {2, 4}, {2, 3}, {3, 2}, {2, 0}, {2, 4}, {2, 3}, {4, 1},
        {2, 0}, {2, 4}, {2, 3}, {3, 2}, {2, 0}, {2, 4}, {2, 3}, {4, 5},
    };
    for (i = simple_code_or_skip; i < kCodeLengthCodes && space > 0; ++i) {
      const int code_len_idx = kCodeLengthCodeOrder[i];
      const HuffmanCode* p = huff;
      uint8_t v;
      BrunsliBitReaderFillWindow(br, 4);
      p += (br->val_ >> br->bit_pos_) & 15;
      br->bit_pos_ += p->bits;
      v = (uint8_t)p->value;
      code_length_code_lengths[code_len_idx] = v;
      if (v != 0) {
        space -= (32u >> v);
        ++num_codes;
      }
    }
    ok = (num_codes == 1 || space == 0) &&
         ReadHuffmanCodeLengths(code_length_code_lengths, alphabet_size,
                                &code_lengths[0], br);
  }
  uint16_t counts[16] = {0};
  for (int i = 0; i < alphabet_size; ++i) {
    ++counts[code_lengths[i]];
  }
  if (ok) {
    table_size = BuildHuffmanTable(&table_[0], kHuffmanTableBits,
                                   &code_lengths[0], alphabet_size, &counts[0]);
  }
  return (table_size > 0);
}

}  // namespace brunsli
