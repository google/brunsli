// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_JPEG_BIT_WRITER_H_
#define BRUNSLI_DEC_JPEG_BIT_WRITER_H_

#include <memory>

#include <brunsli/types.h>

namespace brunsli {

// Returns non-zero if and only if x has a zero byte, i.e. one of
// x & 0xff, x & 0xff00, ..., x & 0xff00000000000000 is zero.
inline uint64_t HasZeroByte(uint64_t x) {
  return (x - 0x0101010101010101ULL) & ~x & 0x8080808080808080ULL;
}

// Handles the packing of bits into output bytes.
struct BitWriter {
  explicit BitWriter(size_t length)
      : len(length),
        data(new uint8_t[len]),
        pos(0),
        put_buffer(0),
        put_bits(64),
        overflow(false),
        invalid_write(false) {}

  void WriteBits(int nbits, uint64_t bits) {
    if (nbits == 0) {
      invalid_write = true;
      return;
    }
    put_bits -= nbits;
    put_buffer |= (bits << put_bits);
    if (put_bits <= 16) {
      // At this point we are ready to emit the most significant 6 bytes of
      // put_buffer_ to the output.
      // The JPEG format requires that after every 0xff byte in the entropy
      // coded section, there is a zero byte, therefore we first check if any of
      // the 6 most significant bytes of put_buffer_ is 0xff.
      if (HasZeroByte(~put_buffer | 0xffff)) {
        // We have a 0xff byte somewhere, examine each byte and append a zero
        // byte if necessary.
        EmitByte((put_buffer >> 56) & 0xff);
        EmitByte((put_buffer >> 48) & 0xff);
        EmitByte((put_buffer >> 40) & 0xff);
        EmitByte((put_buffer >> 32) & 0xff);
        EmitByte((put_buffer >> 24) & 0xff);
        EmitByte((put_buffer >> 16) & 0xff);
      } else if (pos + 6 < len) {
        // We don't have any 0xff bytes, output all 6 bytes without checking.
        data[pos] = (put_buffer >> 56) & 0xff;
        data[pos + 1] = (put_buffer >> 48) & 0xff;
        data[pos + 2] = (put_buffer >> 40) & 0xff;
        data[pos + 3] = (put_buffer >> 32) & 0xff;
        data[pos + 4] = (put_buffer >> 24) & 0xff;
        data[pos + 5] = (put_buffer >> 16) & 0xff;
        pos += 6;
      } else {
        overflow = true;
      }
      put_buffer <<= 48;
      put_bits += 48;
    }
  }

  // Writes the given byte to the output, writes an extra zero if byte is 0xff.
  void EmitByte(int byte) {
    if (pos < len) {
      data[pos++] = byte;
    } else {
      overflow = true;
    }
    if (byte == 0xff) {
      EmitByte(0);
    }
  }

  void EmitMarker(int marker) {
    if (pos + 1 < len) {
      data[pos++] = 0xff;
      data[pos++] = marker;
    } else {
      overflow = true;
    }
  }

  void JumpToByteBoundary(uint8_t pad_pattern) {
    while (put_bits <= 56) {
      int c = (put_buffer >> 56) & 0xff;
      EmitByte(c);
      put_buffer <<= 8;
      put_bits += 8;
    }
    if (put_bits < 64) {
      int pad_mask = 0xFFu >> (64 - put_bits);
      int c = ((put_buffer >> 56) & ~pad_mask) | pad_pattern;
      EmitByte(c);
    }
    put_buffer = 0;
    put_bits = 64;
  }

  size_t len;
  std::unique_ptr<uint8_t[]> data;
  int pos;
  uint64_t put_buffer;
  int put_bits;
  bool overflow;
  bool invalid_write;
};

}  // namespace brunsli

#endif  // BRUNSLI_DEC_JPEG_BIT_WRITER_H_
