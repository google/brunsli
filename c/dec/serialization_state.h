// Copyright (c) Google LLC 2020
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_SERIALIZATION_STATE_H_
#define BRUNSLI_DEC_SERIALIZATION_STATE_H_

#include <deque>
#include <vector>

#include "./jpeg_bit_writer.h"
#include <brunsli/jpeg_data.h>
#include <brunsli/types.h>
#include "./output_chunk.h"

namespace brunsli {

struct HuffmanCodeTable {
  int depth[256];
  int code[256];
};

namespace internal {
namespace dec {

// Maximum number of correction bits to buffer.
const int kJPEGMaxCorrectionBits = 1u << 16;

// Holds data that is buffered between 8x8 blocks in progressive mode.
class DCTCodingState {
 public:
  DCTCodingState() {
    refinement_bits_.reserve(kJPEGMaxCorrectionBits);
  }

  void Init() {
    eob_run_ = 0;
    cur_ac_huff_ = nullptr;
    refinement_bits_.clear();
  }

  // Emit all buffered data to the bit stream using the given Huffman code and
  // bit writer.
  void Flush(BitWriter* bw) {
    if (eob_run_ > 0) {
      int nbits = Log2FloorNonZero(eob_run_);
      int symbol = nbits << 4u;
      bw->WriteBits(cur_ac_huff_->depth[symbol], cur_ac_huff_->code[symbol]);
      if (nbits > 0) {
        bw->WriteBits(nbits, eob_run_ & ((1 << nbits) - 1));
      }
      eob_run_ = 0;
    }
    for (size_t i = 0; i < refinement_bits_.size(); ++i) {
      bw->WriteBits(1, refinement_bits_[i]);
    }
    refinement_bits_.clear();
  }

  // Buffer some more data at the end-of-band (the last non-zero or newly
  // non-zero coefficient within the [Ss, Se] spectral band).
  void BufferEndOfBand(const HuffmanCodeTable& ac_huff,
                       const std::vector<int>* new_bits, BitWriter* bw) {
    if (eob_run_ == 0) {
      cur_ac_huff_ = &ac_huff;
    }
    ++eob_run_;
    if (new_bits) {
      refinement_bits_.insert(refinement_bits_.end(), new_bits->begin(),
                              new_bits->end());
    }
    if (eob_run_ == 0x7fff ||
        refinement_bits_.size() > kJPEGMaxCorrectionBits - kDCTBlockSize + 1) {
      Flush(bw);
    }
  }

 private:
  // The run length of end-of-band symbols in a progressive scan.
  int eob_run_;
  // The huffman table to be used when flushing the state.
  const HuffmanCodeTable* cur_ac_huff_;
  // The sequence of currently buffered refinement bits for a successive
  // approximation scan (one where Ah > 0).
  std::vector<int> refinement_bits_;
};

struct EncodeScanState {
  enum Stage {
    HEAD,
    BODY
  };

  Stage stage = HEAD;

  int mcu_y;
  BitWriter bw;
  coeff_t last_dc_coeff[kMaxComponents] = {0};
  int restarts_to_go;
  int next_restart_marker;
  int block_scan_index;
  DCTCodingState coding_state;
  size_t extra_zero_runs_pos;
  int next_extra_zero_run_index;
};

struct SerializationState {
  enum Stage {
    INIT,
    SERIALIZE_SECTION,
    DONE,
    ERROR,
  };

  Stage stage = INIT;

  std::deque<OutputChunk> output_queue;

  size_t section_index = 0;
  int dht_index = 0;
  int dqt_index = 0;
  int app_index = 0;
  int com_index = 0;
  int data_index = 0;
  int scan_index = 0;
  std::vector<HuffmanCodeTable> dc_huff_table;
  std::vector<HuffmanCodeTable> ac_huff_table;
  const int* pad_bits = nullptr;
  const int* pad_bits_end = nullptr;
  bool seen_dri_marker = false;
  bool is_progressive = false;

  EncodeScanState scan_state;
};

}  // namespace dec
}  // namespace internal
}  // namespace brunsli

#endif  // BRUNSLI_DEC_SERIALIZATION_STATE_H_
