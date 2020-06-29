// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_JPEG_BIT_WRITER_H_
#define BRUNSLI_DEC_JPEG_BIT_WRITER_H_

#include <deque>
#include <memory>

#include "../common/platform.h"
#include <brunsli/types.h>
#include "./output_chunk.h"

namespace brunsli {
namespace internal {
namespace dec {

// Returns non-zero if and only if x has a zero byte, i.e. one of
// x & 0xff, x & 0xff00, ..., x & 0xff00000000000000 is zero.
inline uint64_t HasZeroByte(uint64_t x) {
  return (x - 0x0101010101010101ULL) & ~x & 0x8080808080808080ULL;
}

// Handles the packing of bits into output bytes.
struct BitWriter {
  BitWriter() : chunk(0) {}

  void Init(std::deque<OutputChunk>* output_queue) {
    output = output_queue;
    chunk = OutputChunk(kChunkSize);
    pos = 0;
    put_buffer = 0;
    put_bits = 64;
    healthy = true;
    data = chunk.buffer->data();
  }

  void WriteBits(int nbits, uint64_t bits) {
    // This is an optimization; if everything goes well,
    // then |nbits| is positive; if non-existing Huffman symbol is going to be
    // encoded, its length should be zero; later encoder could check the
    // "health" of BitWriter.
    if (nbits == 0) {
      healthy = false;
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
        Reserve(12);
        // We have a 0xff byte somewhere, examine each byte and append a zero
        // byte if necessary.
        EmitByte((put_buffer >> 56) & 0xff);
        EmitByte((put_buffer >> 48) & 0xff);
        EmitByte((put_buffer >> 40) & 0xff);
        EmitByte((put_buffer >> 32) & 0xff);
        EmitByte((put_buffer >> 24) & 0xff);
        EmitByte((put_buffer >> 16) & 0xff);
      } else {
        Reserve(6);
        // We don't have any 0xff bytes, output all 6 bytes without checking.
        data[pos] = (put_buffer >> 56) & 0xff;
        data[pos + 1] = (put_buffer >> 48) & 0xff;
        data[pos + 2] = (put_buffer >> 40) & 0xff;
        data[pos + 3] = (put_buffer >> 32) & 0xff;
        data[pos + 4] = (put_buffer >> 24) & 0xff;
        data[pos + 5] = (put_buffer >> 16) & 0xff;
        pos += 6;
      }
      put_buffer <<= 48;
      put_bits += 48;
    }
  }

  void EmitMarker(int marker) {
    Reserve(2);
    BRUNSLI_DCHECK(marker != 0xFF);
    data[pos++] = 0xFF;
    data[pos++] = marker;
  }

  bool JumpToByteBoundary(const int** pad_bits, const int* pad_bits_end) {
    size_t n_bits = put_bits & 7u;
    uint8_t pad_pattern;
    if (*pad_bits == nullptr) {
      pad_pattern = (1u << n_bits) - 1;
    } else {
      pad_pattern = 0;
      const int* src = *pad_bits;
      // TODO(eustas): bitwise reading looks insanely ineffective...
      while (n_bits--) {
        pad_pattern <<= 1;
        if (src >= pad_bits_end) return false;
        // TODO(eustas): DCHECK *src == {0, 1}
        pad_pattern |= !!*(src++);
      }
      *pad_bits = src;
    }

    Reserve(16);

    while (put_bits <= 56) {
      int c = (put_buffer >> 56) & 0xFF;
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

    return true;
  }

  bool IsHealthy() { return healthy; }

  void Finish() {
    if (pos > 0) {
      chunk.len = pos;
      output->emplace_back(std::move(chunk));
      chunk = OutputChunk(nullptr, 0);
      data = nullptr;
      pos = 0;
    }
  }

 private:
  /**
   * Writes the given byte to the output, writes an extra zero if byte is 0xFF.
   *
   * This method is "careless" - caller must make sure that there is enough
   * space in the output buffer. Emits up to 2 bytes to buffer.
   */
  void EmitByte(int byte) {
    data[pos++] = byte;
    if (byte == 0xFF) data[pos++] = 0;
  }

  void inline Reserve(size_t n_bytes) {
    if (BRUNSLI_PREDICT_FALSE((pos + n_bytes) > kChunkSize)) {
      chunk.len = pos;
      output->emplace_back(std::move(chunk));
      chunk = OutputChunk(kChunkSize);
      data = chunk.buffer->data();
      pos = 0;
    }
  }

  static const size_t kChunkSize = 16384;

  std::deque<OutputChunk>* output;
  OutputChunk chunk;
  uint8_t* data;
  size_t pos;
  uint64_t put_buffer;
  int put_bits;
  bool healthy;
};

}  // namespace dec
}  // namespace internal
}  // namespace brunsli

#endif  // BRUNSLI_DEC_JPEG_BIT_WRITER_H_
