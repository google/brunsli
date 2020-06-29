// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <brunsli/jpeg_data_writer.h>

#include <cstdlib>
#include <cstring>  /* for memset, memcpy */
#include <deque>
#include <string>
#include <vector>

#include "../common/constants.h"
#include <brunsli/jpeg_data.h>
#include "../common/platform.h"
#include <brunsli/types.h>
#include "./jpeg_bit_writer.h"
#include "./serialization_state.h"
#include "./state.h"
#include "./state_internal.h"

namespace brunsli {

using ::brunsli::internal::dec::BitWriter;
using ::brunsli::internal::dec::DCTCodingState;
using ::brunsli::internal::dec::EncodeScanState;
using ::brunsli::internal::dec::OutputChunk;
using ::brunsli::internal::dec::SerializationState;
using ::brunsli::internal::dec::SerializationStatus;
using ::brunsli::internal::dec::SerializeJpeg;
using ::brunsli::internal::dec::Stage;
using ::brunsli::internal::dec::State;

namespace {

const int kJpegPrecision = 8;

// Returns ceil(a/b).
inline int DivCeil(int a, int b) { return (a + b - 1) / b; }

bool BuildHuffmanCodeTable(const JPEGHuffmanCode& huff,
                           HuffmanCodeTable* table) {
  int huff_code[kJpegHuffmanAlphabetSize];
  // +1 for a sentinel element.
  uint32_t huff_size[kJpegHuffmanAlphabetSize + 1];
  int p = 0;
  for (size_t l = 1; l <= kJpegHuffmanMaxBitLength; ++l) {
    int i = huff.counts[l];
    if (p + i > kJpegHuffmanAlphabetSize + 1) {
      return false;
    }
    while (i--) huff_size[p++] = l;
  }

  if (p == 0) {
    return true;
  }

  // Reuse sentinel element.
  int last_p = p - 1;
  huff_size[last_p] = 0;

  int code = 0;
  uint32_t si = huff_size[0];
  p = 0;
  while (huff_size[p]) {
    while ((huff_size[p]) == si) {
      huff_code[p++] = code;
      code++;
    }
    code <<= 1;
    si++;
  }
  for (p = 0; p < last_p; p++) {
    int i = huff.values[p];
    table->depth[i] = huff_size[p];
    table->code[i] = huff_code[p];
  }
  return true;
}

bool EncodeSOI(SerializationState* state) {
  state->output_queue.push_back(OutputChunk({0xFF, 0xD8}));
  return true;
}

bool EncodeEOI(const JPEGData& jpg, SerializationState* state) {
  state->output_queue.push_back(OutputChunk({0xFF, 0xD9}));
  state->output_queue.emplace_back(jpg.tail_data);
  return true;
}

bool EncodeSOF(const JPEGData& jpg, uint8_t marker, SerializationState* state) {
  if (marker <= 0xC2) state->is_progressive = (marker == 0xC2);

  const size_t n_comps = jpg.components.size();
  const size_t marker_len = 8 + 3 * n_comps;
  state->output_queue.emplace_back(marker_len + 2);
  uint8_t* data = state->output_queue.back().buffer->data();
  size_t pos = 0;
  data[pos++] = 0xFF;
  data[pos++] = marker;
  data[pos++] = marker_len >> 8u;
  data[pos++] = marker_len & 0xFFu;
  data[pos++] = kJpegPrecision;
  data[pos++] = jpg.height >> 8u;
  data[pos++] = jpg.height & 0xFFu;
  data[pos++] = jpg.width >> 8u;
  data[pos++] = jpg.width & 0xFFu;
  data[pos++] = n_comps;
  for (size_t i = 0; i < n_comps; ++i) {
    data[pos++] = jpg.components[i].id;
    data[pos++] = ((jpg.components[i].h_samp_factor << 4u) |
                   (jpg.components[i].v_samp_factor));
    const size_t quant_idx = jpg.components[i].quant_idx;
    if (quant_idx >= jpg.quant.size()) return false;
    data[pos++] = jpg.quant[quant_idx].index;
  }
  return true;
}

bool EncodeSOS(const JPEGData& jpg, const JPEGScanInfo& scan_info,
               SerializationState* state) {
  const size_t n_scans = scan_info.components.size();
  const size_t marker_len = 6 + 2 * n_scans;
  state->output_queue.emplace_back(marker_len + 2);
  uint8_t* data = state->output_queue.back().buffer->data();
  size_t pos = 0;
  data[pos++] = 0xFF;
  data[pos++] = 0xDA;
  data[pos++] = marker_len >> 8u;
  data[pos++] = marker_len & 0xFFu;
  data[pos++] = n_scans;
  for (size_t i = 0; i < n_scans; ++i) {
    const JPEGComponentScanInfo& si = scan_info.components[i];
    if (si.comp_idx >= jpg.components.size()) return false;
    data[pos++] = jpg.components[si.comp_idx].id;
    data[pos++] = (si.dc_tbl_idx << 4u) + si.ac_tbl_idx;
  }
  data[pos++] = scan_info.Ss;
  data[pos++] = scan_info.Se;
  data[pos++] = ((scan_info.Ah << 4u) | (scan_info.Al));
  return true;
}

bool EncodeDHT(const JPEGData& jpg, SerializationState* state) {
  const std::vector<JPEGHuffmanCode>& huffman_code = jpg.huffman_code;

  size_t marker_len = 2;
  for (size_t i = state->dht_index; i < huffman_code.size(); ++i) {
    const JPEGHuffmanCode& huff = huffman_code[i];
    marker_len += kJpegHuffmanMaxBitLength;
    for (size_t j = 0; j < huff.counts.size(); ++j) {
      marker_len += huff.counts[j];
    }
    if (huff.is_last) break;
  }
  state->output_queue.emplace_back(marker_len + 2);
  uint8_t* data = state->output_queue.back().buffer->data();
  size_t pos = 0;
  data[pos++] = 0xFF;
  data[pos++] = 0xC4;
  data[pos++] = marker_len >> 8u;
  data[pos++] = marker_len & 0xFFu;
  while (true) {
    const size_t huffman_code_index = state->dht_index++;
    if (huffman_code_index >= huffman_code.size()) {
      return false;
    }
    const JPEGHuffmanCode& huff = huffman_code[huffman_code_index];
    size_t index = huff.slot_id;
    HuffmanCodeTable* huff_table;
    if (index & 0x10) {
      index -= 0x10;
      huff_table = &state->ac_huff_table[index];
    } else {
      huff_table = &state->dc_huff_table[index];
    }
    // TODO(eustas): cache
    // TODO(eustas): set up non-existing symbols
    if (!BuildHuffmanCodeTable(huff, huff_table)) {
      return false;
    }
    size_t total_count = 0;
    size_t max_length = 0;
    for (size_t i = 0; i < huff.counts.size(); ++i) {
      if (huff.counts[i] != 0) {
        max_length = i;
      }
      total_count += huff.counts[i];
    }
    --total_count;
    data[pos++] = huff.slot_id;
    for (size_t i = 1; i <= kJpegHuffmanMaxBitLength; ++i) {
      data[pos++] = (i == max_length ? huff.counts[i] - 1 : huff.counts[i]);
    }
    for (size_t i = 0; i < total_count; ++i) {
      data[pos++] = huff.values[i];
    }
    if (huff.is_last) break;
  }
  return true;
}

bool EncodeDQT(const JPEGData& jpg, SerializationState* state) {
  int marker_len = 2;
  for (size_t i = state->dqt_index; i < jpg.quant.size(); ++i) {
    const JPEGQuantTable& table = jpg.quant[i];
    marker_len += 1 + (table.precision ? 2 : 1) * kDCTBlockSize;
    if (table.is_last) break;
  }
  state->output_queue.emplace_back(marker_len + 2);
  uint8_t* data = state->output_queue.back().buffer->data();
  size_t pos = 0;
  data[pos++] = 0xFF;
  data[pos++] = 0xDB;
  data[pos++] = marker_len >> 8u;
  data[pos++] = marker_len & 0xFFu;
  while (true) {
    const size_t idx = state->dqt_index++;
    if (idx >= jpg.quant.size()) {
      return false;  // corrupt input
    }
    const JPEGQuantTable& table = jpg.quant[idx];
    data[pos++] = (table.precision << 4u) + table.index;
    for (size_t i = 0; i < kDCTBlockSize; ++i) {
      int val_idx = kJPEGNaturalOrder[i];
      int val = table.values[val_idx];
      if (table.precision) {
        data[pos++] = val >> 8u;
      }
      data[pos++] = val & 0xFFu;
    }
    if (table.is_last) break;
  }
  return true;
}

bool EncodeDRI(const JPEGData& jpg, SerializationState* state) {
  state->seen_dri_marker = true;
  OutputChunk dri_marker = {0xFF,
                            0xDD,
                            0,
                            4,
                            static_cast<uint8_t>(jpg.restart_interval >> 8),
                            static_cast<uint8_t>(jpg.restart_interval & 0xFF)};
  state->output_queue.push_back(std::move(dri_marker));
  return true;
}

bool EncodeRestart(uint8_t marker, SerializationState* state) {
  state->output_queue.push_back(OutputChunk({0xFF, marker}));
  return true;
}

bool EncodeAPP(const JPEGData& jpg, uint8_t marker, SerializationState* state) {
  // TODO(eustas): check that marker corresponds to payload?
  (void)marker;

  size_t app_index = state->app_index++;
  if (app_index >= jpg.app_data.size()) return false;
  state->output_queue.push_back(OutputChunk({0xFF}));
  state->output_queue.emplace_back(jpg.app_data[app_index]);
  return true;
}

bool EncodeCOM(const JPEGData& jpg, SerializationState* state) {
  size_t com_index = state->com_index++;
  if (com_index >= jpg.com_data.size()) return false;
  state->output_queue.push_back(OutputChunk({0xFF}));
  state->output_queue.emplace_back(jpg.com_data[com_index]);
  return true;
}

bool EncodeInterMarkerData(const JPEGData& jpg, SerializationState* state) {
  size_t index = state->data_index++;
  if (index >= jpg.inter_marker_data.size()) return false;
  state->output_queue.emplace_back(jpg.inter_marker_data[index]);
  return true;
}

bool EncodeDCTBlockSequential(const coeff_t* coeffs,
                              const HuffmanCodeTable& dc_huff,
                              const HuffmanCodeTable& ac_huff,
                              int num_zero_runs, coeff_t* last_dc_coeff,
                              BitWriter* bw) {
  coeff_t temp2;
  coeff_t temp;
  temp2 = coeffs[0];
  temp = temp2 - *last_dc_coeff;
  *last_dc_coeff = temp2;
  temp2 = temp;
  if (temp < 0) {
    temp = -temp;
    temp2--;
  }
  int dc_nbits = (temp == 0) ? 0 : (Log2FloorNonZero(temp) + 1);
  bw->WriteBits(dc_huff.depth[dc_nbits], dc_huff.code[dc_nbits]);
  if (dc_nbits > 0) {
    bw->WriteBits(dc_nbits, temp2 & ((1u << dc_nbits) - 1));
  }
  int r = 0;
  for (int k = 1; k < 64; ++k) {
    if ((temp = coeffs[kJPEGNaturalOrder[k]]) == 0) {
      r++;
      continue;
    }
    if (temp < 0) {
      temp = -temp;
      temp2 = ~temp;
    } else {
      temp2 = temp;
    }
    while (r > 15) {
      bw->WriteBits(ac_huff.depth[0xf0], ac_huff.code[0xf0]);
      r -= 16;
    }
    int ac_nbits = Log2FloorNonZero(temp) + 1;
    int symbol = (r << 4u) + ac_nbits;
    bw->WriteBits(ac_huff.depth[symbol], ac_huff.code[symbol]);
    bw->WriteBits(ac_nbits, temp2 & ((1 << ac_nbits) - 1));
    r = 0;
  }
  for (int i = 0; i < num_zero_runs; ++i) {
    bw->WriteBits(ac_huff.depth[0xf0], ac_huff.code[0xf0]);
    r -= 16;
  }
  if (r > 0) {
    bw->WriteBits(ac_huff.depth[0], ac_huff.code[0]);
  }
  return true;
}

bool EncodeDCTBlockProgressive(const coeff_t* coeffs,
                               const HuffmanCodeTable& dc_huff,
                               const HuffmanCodeTable& ac_huff, int Ss, int Se,
                               int Al, int num_zero_runs,
                               DCTCodingState* coding_state,
                               coeff_t* last_dc_coeff, BitWriter* bw) {
  bool eob_run_allowed = Ss > 0;
  coeff_t temp2;
  coeff_t temp;
  if (Ss == 0) {
    temp2 = coeffs[0] >> Al;
    temp = temp2 - *last_dc_coeff;
    *last_dc_coeff = temp2;
    temp2 = temp;
    if (temp < 0) {
      temp = -temp;
      temp2--;
    }
    int nbits = (temp == 0) ? 0 : (Log2FloorNonZero(temp) + 1);
    bw->WriteBits(dc_huff.depth[nbits], dc_huff.code[nbits]);
    if (nbits > 0) {
      bw->WriteBits(nbits, temp2 & ((1 << nbits) - 1));
    }
    ++Ss;
  }
  if (Ss > Se) {
    return true;
  }
  int r = 0;
  for (int k = Ss; k <= Se; ++k) {
    if ((temp = coeffs[kJPEGNaturalOrder[k]]) == 0) {
      r++;
      continue;
    }
    if (temp < 0) {
      temp = -temp;
      temp >>= Al;
      temp2 = ~temp;
    } else {
      temp >>= Al;
      temp2 = temp;
    }
    if (temp == 0) {
      r++;
      continue;
    }
    coding_state->Flush(bw);
    while (r > 15) {
      bw->WriteBits(ac_huff.depth[0xf0], ac_huff.code[0xf0]);
      r -= 16;
    }
    int nbits = Log2FloorNonZero(temp) + 1;
    int symbol = (r << 4u) + nbits;
    bw->WriteBits(ac_huff.depth[symbol], ac_huff.code[symbol]);
    bw->WriteBits(nbits, temp2 & ((1 << nbits) - 1));
    r = 0;
  }
  if (num_zero_runs > 0) {
    coding_state->Flush(bw);
    for (int i = 0; i < num_zero_runs; ++i) {
      bw->WriteBits(ac_huff.depth[0xf0], ac_huff.code[0xf0]);
      r -= 16;
    }
  }
  if (r > 0) {
    coding_state->BufferEndOfBand(ac_huff, nullptr, bw);
    if (!eob_run_allowed) {
      coding_state->Flush(bw);
    }
  }
  return true;
}

bool EncodeRefinementBits(const coeff_t* coeffs,
                          const HuffmanCodeTable& ac_huff, int Ss, int Se,
                          int Al, DCTCodingState* coding_state, BitWriter* bw) {
  bool eob_run_allowed = Ss > 0;
  if (Ss == 0) {
    // Emit next bit of DC component.
    bw->WriteBits(1, (coeffs[0] >> Al) & 1);
    ++Ss;
  }
  if (Ss > Se) {
    return true;
  }
  int abs_values[kDCTBlockSize];
  int eob = 0;
  for (int k = Ss; k <= Se; k++) {
    const coeff_t abs_val = std::abs(coeffs[kJPEGNaturalOrder[k]]);
    abs_values[k] = abs_val >> Al;
    if (abs_values[k] == 1) {
      eob = k;
    }
  }
  int r = 0;
  std::vector<int> refinement_bits;
  refinement_bits.reserve(kDCTBlockSize);
  for (int k = Ss; k <= Se; k++) {
    if (abs_values[k] == 0) {
      r++;
      continue;
    }
    while (r > 15 && k <= eob) {
      coding_state->Flush(bw);
      bw->WriteBits(ac_huff.depth[0xf0], ac_huff.code[0xf0]);
      r -= 16;
      for (int bit : refinement_bits) {
        bw->WriteBits(1, bit);
      }
      refinement_bits.clear();
    }
    if (abs_values[k] > 1) {
      refinement_bits.push_back(abs_values[k] & 1u);
      continue;
    }
    coding_state->Flush(bw);
    int symbol = (r << 4u) + 1;
    int new_non_zero_bit = (coeffs[kJPEGNaturalOrder[k]] < 0) ? 0 : 1;
    bw->WriteBits(ac_huff.depth[symbol], ac_huff.code[symbol]);
    bw->WriteBits(1, new_non_zero_bit);
    for (int bit : refinement_bits) {
      bw->WriteBits(1, bit);
    }
    refinement_bits.clear();
    r = 0;
  }
  if (r > 0 || !refinement_bits.empty()) {
    coding_state->BufferEndOfBand(ac_huff, &refinement_bits, bw);
    if (!eob_run_allowed) {
      coding_state->Flush(bw);
    }
  }
  return true;
}

SerializationStatus EncodeScan(const JPEGData& jpg, const State& parsing_state,
                               SerializationState* state) {
  const JPEGScanInfo& scan_info = jpg.scan_info[state->scan_index];
  EncodeScanState& ss = state->scan_state;

  const int restart_interval =
      state->seen_dri_marker ? jpg.restart_interval : 0;

  const auto get_next_extra_zero_run_index = [&ss, &scan_info]() {
    if (ss.extra_zero_runs_pos < scan_info.extra_zero_runs.size()) {
      return scan_info.extra_zero_runs[ss.extra_zero_runs_pos].block_idx;
    } else {
      return -1;
    }
  };

  if (ss.stage == EncodeScanState::HEAD) {
    if (!EncodeSOS(jpg, scan_info, state)) return SerializationStatus::ERROR;
    ss.bw.Init(&state->output_queue);
    ss.coding_state.Init();
    ss.restarts_to_go = restart_interval;
    ss.next_restart_marker = 0;
    ss.block_scan_index = 0;
    ss.extra_zero_runs_pos = 0;
    ss.next_extra_zero_run_index = get_next_extra_zero_run_index();
    ss.mcu_y = 0;
    memset(ss.last_dc_coeff, 0, sizeof(ss.last_dc_coeff));
    ss.stage = EncodeScanState::BODY;
  }
  BitWriter* bw = &ss.bw;
  DCTCodingState* coding_state = &ss.coding_state;

  BRUNSLI_DCHECK(ss.stage == EncodeScanState::BODY);

  // "Non-interleaved" means color data comes in separate scans, in other words
  // each scan can contain only one color component.
  const bool is_interleaved = (scan_info.components.size() > 1);
  const JPEGComponent& base_component =
      jpg.components[scan_info.components[0].comp_idx];
  // h_group / v_group act as numerators for converting number of blocks to
  // number of MCU. In interleaved mode it is 1, so MCU is represented with
  // max_*_samp_factor blocks. In non-interleaved mode we choose numerator to
  // be the samping factor, consequently MCU is always represented with single
  // block.
  const int h_group = is_interleaved ? 1 : base_component.h_samp_factor;
  const int v_group = is_interleaved ? 1 : base_component.v_samp_factor;
  const int MCUs_per_row =
      DivCeil(jpg.width * h_group, 8 * jpg.max_h_samp_factor);
  const int MCU_rows = DivCeil(jpg.height * v_group, 8 * jpg.max_v_samp_factor);
  const bool is_progressive = state->is_progressive;
  const int Al = is_progressive ? scan_info.Al : 0;
  const int Ah = is_progressive ? scan_info.Ah : 0;
  const int Ss = is_progressive ? scan_info.Ss : 0;
  const int Se = is_progressive ? scan_info.Se : 63;
  const bool need_sequential =
      !is_progressive || (Ah == 0 && Al == 0 && Ss == 0 && Se == 63);

  // DC-only is defined by [0..0] spectral range.
  const bool want_ac = ((Ss != 0) || (Se != 0));
  const bool complete_ac = (parsing_state.stage == Stage::DONE);
  const bool has_ac =
      complete_ac || HasSection(&parsing_state, kBrunsliACDataTag);
  if (want_ac && !has_ac) return SerializationStatus::NEEDS_MORE_INPUT;

  // |has_ac| implies |complete_dc| but not vice versa; for the sake of
  // simplicity we pretend they are equal, because they are separated by just a
  // few bytes of input.
  const bool complete_dc = has_ac;
  const bool complete = want_ac ? complete_ac : complete_dc;
  // When "incomplete" |ac_dc| tracks information about current ("incomplete")
  // band parsing progress.
  const int last_mcu_y =
      complete ? MCU_rows : parsing_state.internal->ac_dc.next_mcu_y * v_group;

  for (; ss.mcu_y < last_mcu_y; ++ss.mcu_y) {
    for (int mcu_x = 0; mcu_x < MCUs_per_row; ++mcu_x) {
      // Possibly emit a restart marker.
      if (restart_interval > 0 && ss.restarts_to_go == 0) {
        coding_state->Flush(bw);
        if (!bw->JumpToByteBoundary(&state->pad_bits, state->pad_bits_end)) {
          return SerializationStatus::ERROR;
        }
        bw->EmitMarker(0xD0 + ss.next_restart_marker);
        ss.next_restart_marker += 1;
        ss.next_restart_marker &= 0x7;
        ss.restarts_to_go = restart_interval;
        memset(ss.last_dc_coeff, 0, sizeof(ss.last_dc_coeff));
      }
      // Encode one MCU
      for (size_t i = 0; i < scan_info.components.size(); ++i) {
        const JPEGComponentScanInfo& si = scan_info.components[i];
        const JPEGComponent& c = jpg.components[si.comp_idx];
        const HuffmanCodeTable& dc_huff = state->dc_huff_table[si.dc_tbl_idx];
        const HuffmanCodeTable& ac_huff = state->ac_huff_table[si.ac_tbl_idx];
        int n_blocks_y = is_interleaved ? c.v_samp_factor : 1;
        int n_blocks_x = is_interleaved ? c.h_samp_factor : 1;
        for (int iy = 0; iy < n_blocks_y; ++iy) {
          for (int ix = 0; ix < n_blocks_x; ++ix) {
            int block_y = ss.mcu_y * n_blocks_y + iy;
            int block_x = mcu_x * n_blocks_x + ix;
            int block_idx = block_y * c.width_in_blocks + block_x;
            if (scan_info.reset_points.find(ss.block_scan_index) !=
                scan_info.reset_points.end()) {
              coding_state->Flush(bw);
            }
            int num_zero_runs = 0;
            if (ss.block_scan_index == ss.next_extra_zero_run_index) {
              num_zero_runs = scan_info.extra_zero_runs[ss.extra_zero_runs_pos]
                                  .num_extra_zero_runs;
              ++ss.extra_zero_runs_pos;
              ss.next_extra_zero_run_index = get_next_extra_zero_run_index();
            }
            const coeff_t* coeffs = &c.coeffs[block_idx << 6];
            bool ok;
            if (need_sequential) {
              ok = EncodeDCTBlockSequential(coeffs, dc_huff, ac_huff,
                                            num_zero_runs,
                                            ss.last_dc_coeff + si.comp_idx, bw);
            } else if (Ah == 0) {
              ok = EncodeDCTBlockProgressive(
                  coeffs, dc_huff, ac_huff, Ss, Se, Al, num_zero_runs,
                  coding_state, ss.last_dc_coeff + si.comp_idx, bw);
            } else {
              ok = EncodeRefinementBits(coeffs, ac_huff, Ss, Se, Al,
                                        coding_state, bw);
            }
            if (!ok) return SerializationStatus::ERROR;
            ++ss.block_scan_index;
          }
        }
      }
      --ss.restarts_to_go;
    }
  }
  if (ss.mcu_y < MCU_rows) {
    if (!bw->IsHealthy()) return SerializationStatus::ERROR;
    return SerializationStatus::NEEDS_MORE_INPUT;
  }
  coding_state->Flush(bw);
  if (!bw->JumpToByteBoundary(&state->pad_bits, state->pad_bits_end)) {
    return SerializationStatus::ERROR;
  }
  bw->Finish();
  ss.stage = EncodeScanState::HEAD;
  state->scan_index++;
  if (!bw->IsHealthy()) return SerializationStatus::ERROR;

  return SerializationStatus::DONE;
}

SerializationStatus SerializeSection(uint8_t marker, const State& parsing_state,
                                     SerializationState* state,
                                     const JPEGData& jpg) {
  const auto to_status = [] (bool result) {
    return result ? SerializationStatus::DONE : SerializationStatus::ERROR;
  };
  // TODO(eustas): add and use marker enum
  switch (marker) {
    case 0xC0:
    case 0xC1:
    case 0xC2:
    case 0xC9:
    case 0xCA:
      return to_status(EncodeSOF(jpg, marker, state));

    case 0xC4:
      return to_status(EncodeDHT(jpg, state));

    case 0xD0:
    case 0xD1:
    case 0xD2:
    case 0xD3:
    case 0xD4:
    case 0xD5:
    case 0xD6:
    case 0xD7:
      return to_status(EncodeRestart(marker, state));

    case 0xD9:
      return to_status(EncodeEOI(jpg, state));

    case 0xDA:
      return EncodeScan(jpg, parsing_state, state);

    case 0xDB:
      return to_status(EncodeDQT(jpg, state));

    case 0xDD:
      return to_status(EncodeDRI(jpg, state));

    case 0xE0:
    case 0xE1:
    case 0xE2:
    case 0xE3:
    case 0xE4:
    case 0xE5:
    case 0xE6:
    case 0xE7:
    case 0xE8:
    case 0xE9:
    case 0xEA:
    case 0xEB:
    case 0xEC:
    case 0xED:
    case 0xEE:
    case 0xEF:
      return to_status(EncodeAPP(jpg, marker, state));

    case 0xFE:
      return to_status(EncodeCOM(jpg, state));

    case 0xFF:
      return to_status(EncodeInterMarkerData(jpg, state));

    default:
      return SerializationStatus::ERROR;
  }
}

void PushOutput(std::deque<OutputChunk>* in, size_t* available_out,
                uint8_t** next_out) {
  while (*available_out > 0) {
    // No more data.
    if (in->empty()) return;
    OutputChunk& chunk = in->front();
    size_t to_copy = std::min(*available_out, chunk.len);
    if (to_copy > 0) {
      memcpy(*next_out, chunk.next, to_copy);
      *next_out += to_copy;
      *available_out -= to_copy;
      chunk.next += to_copy;
      chunk.len -= to_copy;
    }
    if (chunk.len == 0) in->pop_front();
  }
}

}  // namespace

// Adaptor for old API users. Will be removed once new API will support proper
// streaming serialization.
bool WriteJpeg(const JPEGData& jpg, JPEGOutput out) {
  State state;
  state.stage = Stage::DONE;
  std::vector<uint8_t> buffer(16384);
  while (true) {
    uint8_t* next_out = buffer.data();
    size_t available_out = buffer.size();
    SerializationStatus status =
        SerializeJpeg(&state, jpg, &available_out, &next_out);
    if (status != SerializationStatus::DONE &&
        status != SerializationStatus::NEEDS_MORE_OUTPUT) {
      return false;
    }
    size_t to_write = buffer.size() - available_out;
    if (!out.Write(buffer.data(), to_write)) return false;
    if (status == SerializationStatus::DONE) return true;
  }
}

namespace internal {
namespace dec {
SerializationStatus SerializeJpeg(State* state, const JPEGData& jpg,
                                  size_t* available_out, uint8_t** next_out) {
  SerializationState& ss = state->internal->serialization;

  const auto maybe_push_output = [&] () {
    if (ss.stage != SerializationState::ERROR) {
      PushOutput(&ss.output_queue, available_out, next_out);
    }
  };

  // Push remaining output from prevoius session.
  maybe_push_output();

  while (true) {
    switch (ss.stage) {
      case SerializationState::INIT: {
        // If parsing is complete, serialization is possible.
        bool can_start_serialization = (state->stage == Stage::DONE);
        // Parsing of AC/DC has started; i.e. quant/huffman/metadata is ready
        // to be used.
        if (HasSection(state, kBrunsliDCDataTag) ||
            HasSection(state, kBrunsliACDataTag)) {
          can_start_serialization = true;
        }
        if (!can_start_serialization) {
          return SerializationStatus::NEEDS_MORE_INPUT;
        }
        // JpegBypass is a very simple / special case.
        if (jpg.version == 1) {
          if (jpg.original_jpg == nullptr) {
            ss.stage = SerializationState::ERROR;
            break;
          }
          // TODO(eustas): investigate if bad things can happen when complete
          //               file is passed to parser, but it is impossible to
          //               push complete output.
          ss.output_queue.emplace_back(jpg.original_jpg, jpg.original_jpg_size);
          ss.stage = SerializationState::DONE;
          break;
        }

        // Valid Brunsli requires, at least, 0xD9 marker.
        // This might happen on corrupted stream, or on unconditioned JPEGData.
        // TODO(eustas): check D9 in the only one and is the last one.
        if (jpg.marker_order.empty()) {
          ss.stage = SerializationState::ERROR;
          break;
        }

        ss.dc_huff_table.resize(kMaxHuffmanTables);
        ss.ac_huff_table.resize(kMaxHuffmanTables);
        if (jpg.has_zero_padding_bit) {
          ss.pad_bits = jpg.padding_bits.data();
          ss.pad_bits_end = ss.pad_bits + jpg.padding_bits.size();
        }

        EncodeSOI(&ss);
        maybe_push_output();
        ss.stage = SerializationState::SERIALIZE_SECTION;
        break;
      }

      case SerializationState::SERIALIZE_SECTION: {
        if (ss.section_index >= jpg.marker_order.size()) {
          ss.stage = SerializationState::DONE;
          break;
        }
        uint8_t marker = jpg.marker_order[ss.section_index];
        SerializationStatus status = SerializeSection(marker, *state, &ss, jpg);
        if (status == SerializationStatus::ERROR) {
          BRUNSLI_LOG_DEBUG() << "Failed to encode marker " << std::hex
                              << marker << BRUNSLI_ENDL();
          ss.stage = SerializationState::ERROR;
          break;
        }
        maybe_push_output();
        if (status == SerializationStatus::NEEDS_MORE_INPUT) {
          return SerializationStatus::NEEDS_MORE_INPUT;
        } else if (status != SerializationStatus::DONE) {
          BRUNSLI_DCHECK(false);
          ss.stage = SerializationState::ERROR;
          break;
        }
        ++ss.section_index;
        break;
      }

      case SerializationState::DONE: {
        if (!ss.output_queue.empty()) {
          return SerializationStatus::NEEDS_MORE_OUTPUT;
        } else {
          return SerializationStatus::DONE;
        }
      }

      default:
        return SerializationStatus::ERROR;
    }
  }
}
}  // namespace dec
}  // namespace internal

}  // namespace brunsli
