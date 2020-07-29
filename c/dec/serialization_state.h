// Copyright (c) Google LLC 2020
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef BRUNSLI_DEC_SERIALIZATION_STATE_H_
#define BRUNSLI_DEC_SERIALIZATION_STATE_H_

#include <deque>
#include <vector>

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

// Handles the packing of bits into output bytes.
struct BitWriter {
  bool healthy;
  std::deque<OutputChunk>* output;
  OutputChunk chunk;
  uint8_t* data;
  size_t pos;
  uint64_t put_buffer;
  int put_bits;
};

// Holds data that is buffered between 8x8 blocks in progressive mode.
struct DCTCodingState {
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
  size_t next_reset_point_pos;
  int next_reset_point;
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
