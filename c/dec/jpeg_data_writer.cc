// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <brunsli/jpeg_data_writer.h>

#include <cstdlib>
#include <cstring>  /* for memset, memcpy */
#include <string>
#include <vector>

#include "../common/huffman_tree.h"
#include "../common/platform.h"
#include <brunsli/types.h>
#include "./jpeg_bit_writer.h"

namespace brunsli {

namespace {

const int kJpegPrecision = 8;

// Maximum number of correction bits to buffer.
const int kJPEGMaxCorrectionBits = 1u << 16;

// Returns ceil(a/b).
inline int DivCeil(int a, int b) { return (a + b - 1) / b; }

struct HuffmanCodeTable {
  int depth[256];
  int code[256];
};

bool BuildHuffmanCodeTable(const JPEGHuffmanCode& huff,
                           HuffmanCodeTable* table) {
  int huff_code[kJpegHuffmanAlphabetSize];
  // +1 for a sentinel element.
  int huff_size[kJpegHuffmanAlphabetSize + 1];
  int p = 0;
  for (int l = 1; l <= kJpegHuffmanMaxBitLength; ++l) {
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
  int si = huff_size[0];
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

// Writes len bytes from buf, using the out callback.
inline bool JPEGWrite(JPEGOutput out, const uint8_t* buf, size_t len) {
  static const size_t kBlockSize = 1u << 30u;
  size_t pos = 0;
  while (len - pos > kBlockSize) {
    if (!out.Write(buf + pos, kBlockSize)) {
      return false;
    }
    pos += kBlockSize;
  }
  return out.Write(buf + pos, len - pos);
}

// Writes a string using the out callback.
inline bool JPEGWrite(JPEGOutput out, const std::string& s) {
  const uint8_t* data = reinterpret_cast<const uint8_t*>(&s[0]);
  return JPEGWrite(out, data, s.size());
}

bool EncodeSOF(const JPEGData& jpg, uint8_t marker, JPEGOutput out) {
  const size_t n_comps = jpg.components.size();
  const size_t marker_len = 8 + 3 * n_comps;
  std::vector<uint8_t> data(marker_len + 2);
  size_t pos = 0;
  data[pos++] = 0xff;
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
    if (quant_idx >= jpg.quant.size()) {
      return false;
    }
    data[pos++] = jpg.quant[quant_idx].index;
  }
  return JPEGWrite(out, &data[0], pos);
}

bool EncodeSOS(const JPEGData& jpg, const JPEGScanInfo& scan_info,
               JPEGOutput out) {
  const size_t n_scans = scan_info.components.size();
  const size_t marker_len = 6 + 2 * n_scans;
  std::vector<uint8_t> data(marker_len + 2);
  size_t pos = 0;
  data[pos++] = 0xFF;
  data[pos++] = 0xDA;
  data[pos++] = marker_len >> 8u;
  data[pos++] = marker_len & 0xFFu;
  data[pos++] = n_scans;
  for (size_t i = 0; i < n_scans; ++i) {
    const JPEGComponentScanInfo& si = scan_info.components[i];
    if (si.comp_idx >= jpg.components.size()) {
      return false;
    }
    data[pos++] = jpg.components[si.comp_idx].id;
    data[pos++] = (si.dc_tbl_idx << 4u) + si.ac_tbl_idx;
  }
  data[pos++] = scan_info.Ss;
  data[pos++] = scan_info.Se;
  data[pos++] = ((scan_info.Ah << 4u) | (scan_info.Al));
  return JPEGWrite(out, &data[0], pos);
}

bool EncodeDHT(const std::vector<JPEGHuffmanCode>& huffman_code, int* dht_index,
               JPEGOutput out, std::vector<HuffmanCodeTable>* dc_huff_table,
               std::vector<HuffmanCodeTable>* ac_huff_table) {
  size_t marker_len = 2;
  for (size_t i = *dht_index; i < huffman_code.size(); ++i) {
    const JPEGHuffmanCode& huff = huffman_code[i];
    marker_len += kJpegHuffmanMaxBitLength;
    for (size_t j = 0; j < huff.counts.size(); ++j) {
      marker_len += huff.counts[j];
    }
    if (huff.is_last) break;
  }
  std::vector<uint8_t> data(marker_len + 2);
  size_t pos = 0;
  data[pos++] = 0xff;
  data[pos++] = 0xc4;
  data[pos++] = marker_len >> 8u;
  data[pos++] = marker_len & 0xFFu;
  while (true) {
    const size_t huffman_code_index = (*dht_index)++;
    if (huffman_code_index >= huffman_code.size()) {
      return false;
    }
    const JPEGHuffmanCode& huff = huffman_code[huffman_code_index];
    size_t index = huff.slot_id;
    HuffmanCodeTable* huff_table;
    if (index & 0x10) {
      index -= 0x10;
      huff_table = &(*ac_huff_table)[index];
    } else {
      huff_table = &(*dc_huff_table)[index];
    }
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
  return JPEGWrite(out, &data[0], pos);
}

bool EncodeDQT(const JPEGData& jpg, int* dqt_index, JPEGOutput out) {
  int marker_len = 2;
  for (size_t i = *dqt_index; i < jpg.quant.size(); ++i) {
    const JPEGQuantTable& table = jpg.quant[i];
    marker_len += 1 + (table.precision ? 2 : 1) * kDCTBlockSize;
    if (table.is_last) break;
  }
  std::vector<uint8_t> data(marker_len + 2);
  size_t pos = 0;
  data[pos++] = 0xFF;
  data[pos++] = 0xDB;
  data[pos++] = marker_len >> 8u;
  data[pos++] = marker_len & 0xFFu;
  while (true) {
    const size_t idx = (*dqt_index)++;
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
  return JPEGWrite(out, &data[0], pos);
}

bool EncodeDRI(int restart_interval, JPEGOutput out) {
  uint8_t data[6] = {0xFF, 0xDD, 0, 4};
  data[4] = restart_interval >> 8u;
  data[5] = restart_interval & 0xFFu;
  return JPEGWrite(out, data, sizeof(data));
}

bool EncodeAPP(const JPEGData& jpg, size_t app_index, JPEGOutput out) {
  if (app_index >= jpg.app_data.size()) {
    return false;
  }
  uint8_t data[1] = {0xff};
  return (JPEGWrite(out, data, sizeof(data)) &&
          JPEGWrite(out, jpg.app_data[app_index]));
}

bool EncodeCOM(const JPEGData& jpg, size_t com_index, JPEGOutput out) {
  if (com_index >= jpg.com_data.size()) {
    return false;
  }
  uint8_t data[2] = {0xff, 0xfe};
  return (JPEGWrite(out, data, sizeof(data)) &&
          JPEGWrite(out, jpg.com_data[com_index]));
}

bool EncodeInterMarkerData(const JPEGData& jpg, size_t index, JPEGOutput out) {
  if (index >= jpg.inter_marker_data.size()) {
    return false;
  }
  return JPEGWrite(out, jpg.inter_marker_data[index]);
}

// Holds data that is buffered between 8x8 blocks in progressive mode.
class DCTCodingState {
 public:
  DCTCodingState() : eob_run_(0), cur_ac_huff_(nullptr) {
    refinement_bits_.reserve(kJPEGMaxCorrectionBits);
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

bool GetNextPadPattern(const int** pad_bits, const int* pad_bits_end, int nbits,
                       uint8_t* pad_pattern) {
  // TODO: DCHECK pad_bits < 8
  if (*pad_bits == nullptr) {
    *pad_pattern = (1u << nbits) - 1;
    return true;
  }
  uint8_t p = 0;
  const int* src = *pad_bits;
  // TODO: bitwise reading looks insanely ineffective...
  while (nbits--) {
    p <<= 1;
    if (src >= pad_bits_end) {
      return false;
    }
    // TODO: DCHECK *src == {0, 1}
    p |= *(src++);
  }
  *pad_bits = src;
  *pad_pattern = p;
  return true;
}

bool EncodeScan(const JPEGData& jpg, const JPEGScanInfo& scan_info,
                bool is_progressive,
                const std::vector<HuffmanCodeTable>& dc_huff_table,
                const std::vector<HuffmanCodeTable>& ac_huff_table,
                const int restart_interval, const int** pad_bits,
                const int* pad_bits_end, JPEGOutput out) {
  if (!EncodeSOS(jpg, scan_info, out)) {
    return false;
  }
  bool is_interleaved = (scan_info.components.size() > 1);
  int MCUs_per_row;
  int MCU_rows;
  if (is_interleaved) {
    MCUs_per_row = DivCeil(jpg.width, jpg.max_h_samp_factor * 8);
    MCU_rows = DivCeil(jpg.height, jpg.max_v_samp_factor * 8);
  } else {
    const JPEGComponent& c = jpg.components[scan_info.components[0].comp_idx];
    MCUs_per_row =
        DivCeil(jpg.width * c.h_samp_factor, 8 * jpg.max_h_samp_factor);
    MCU_rows = DivCeil(jpg.height * c.v_samp_factor, 8 * jpg.max_v_samp_factor);
  }
  coeff_t last_dc_coeff[kMaxComponents] = {0};
  BitWriter bw(1 << 17);
  int restarts_to_go = restart_interval;
  int next_restart_marker = 0;
  int block_scan_index = 0;
  int extra_zero_runs_pos = 0;
  int next_extra_zero_run_index = scan_info.extra_zero_runs.empty()
                                      ? -1
                                      : scan_info.extra_zero_runs[0].block_idx;
  DCTCodingState coding_state;
  const int Al = is_progressive ? scan_info.Al : 0;
  const int Ah = is_progressive ? scan_info.Ah : 0;
  const int Ss = is_progressive ? scan_info.Ss : 0;
  const int Se = is_progressive ? scan_info.Se : 63;
  const bool need_sequential =
      !is_progressive || (Ah == 0 && Al == 0 && Ss == 0 && Se == 63);
  for (int mcu_y = 0; mcu_y < MCU_rows; ++mcu_y) {
    for (int mcu_x = 0; mcu_x < MCUs_per_row; ++mcu_x) {
      // Possibly emit a restart marker.
      if (restart_interval > 0 && restarts_to_go == 0) {
        coding_state.Flush(&bw);
        const int n_pad_bits = bw.put_bits & 7u;
        uint8_t pad_pattern;
        if (!GetNextPadPattern(pad_bits, pad_bits_end, n_pad_bits,
                               &pad_pattern)) {
          return false;
        }
        bw.JumpToByteBoundary(pad_pattern);
        bw.EmitMarker(0xd0 + next_restart_marker);
        next_restart_marker += 1;
        next_restart_marker &= 0x7;
        restarts_to_go = restart_interval;
        memset(last_dc_coeff, 0, sizeof(last_dc_coeff));
      }
      // Encode one MCU
      for (size_t i = 0; i < scan_info.components.size(); ++i) {
        const JPEGComponentScanInfo& si = scan_info.components[i];
        const JPEGComponent& c = jpg.components[si.comp_idx];
        const HuffmanCodeTable& dc_huff = dc_huff_table[si.dc_tbl_idx];
        const HuffmanCodeTable& ac_huff = ac_huff_table[si.ac_tbl_idx];
        int n_blocks_y = is_interleaved ? c.v_samp_factor : 1;
        int n_blocks_x = is_interleaved ? c.h_samp_factor : 1;
        for (int iy = 0; iy < n_blocks_y; ++iy) {
          for (int ix = 0; ix < n_blocks_x; ++ix) {
            int block_y = mcu_y * n_blocks_y + iy;
            int block_x = mcu_x * n_blocks_x + ix;
            int block_idx = block_y * c.width_in_blocks + block_x;
            if (scan_info.reset_points.find(block_scan_index) !=
                scan_info.reset_points.end()) {
              coding_state.Flush(&bw);
            }
            int num_zero_runs = 0;
            if (block_scan_index == next_extra_zero_run_index) {
              num_zero_runs = scan_info.extra_zero_runs[extra_zero_runs_pos]
                                  .num_extra_zero_runs;
              ++extra_zero_runs_pos;
              next_extra_zero_run_index =
                  (extra_zero_runs_pos < scan_info.extra_zero_runs.size())
                      ? scan_info.extra_zero_runs[extra_zero_runs_pos].block_idx
                      : -1;
            }
            const coeff_t* coeffs = &c.coeffs[block_idx << 6];
            if (need_sequential) {
              if (!EncodeDCTBlockSequential(coeffs, dc_huff, ac_huff,
                                            num_zero_runs,
                                            &last_dc_coeff[si.comp_idx], &bw)) {
                return false;
              }
            } else if (Ah == 0) {
              if (!EncodeDCTBlockProgressive(
                      coeffs, dc_huff, ac_huff, Ss, Se, Al, num_zero_runs,
                      &coding_state, &last_dc_coeff[si.comp_idx], &bw)) {
                return false;
              }
            } else {
              if (!EncodeRefinementBits(coeffs, ac_huff, Ss, Se, Al,
                                        &coding_state, &bw)) {
                return false;
              }
            }
            ++block_scan_index;
          }
        }
      }
      --restarts_to_go;
      if (bw.pos > (1 << 16)) {
        if (!JPEGWrite(out, bw.data.get(), bw.pos)) {
          return false;
        }
        bw.pos = 0;
      }
    }
  }
  coding_state.Flush(&bw);
  const int n_pad_bits = bw.put_bits & 7u;
  uint8_t pad_pattern;
  if (!GetNextPadPattern(pad_bits, pad_bits_end, n_pad_bits, &pad_pattern)) {
    return false;
  }
  bw.JumpToByteBoundary(pad_pattern);
  return !bw.overflow && !bw.invalid_write &&
         JPEGWrite(out, bw.data.get(), bw.pos);
}

}  // namespace

bool WriteJpegBypass(const JPEGData& jpg, JPEGOutput out) {
  if (jpg.version != 1 || jpg.original_jpg == nullptr) {
    return false;
  }
  return JPEGWrite(out, jpg.original_jpg, jpg.original_jpg_size);
}

bool WriteJpeg(const JPEGData& jpg, JPEGOutput out) {
  if (jpg.version == 1) {
    return WriteJpegBypass(jpg, out);
  }

  const std::vector<uint8_t>& marker_order = jpg.marker_order;
  // Valid Brunsli requires, at least, 0xD9 marker.
  // This might happen on corrupted stream, or on unconditioned JPEGData.
  if (marker_order.empty()) {
    return false;
  }

  static const uint8_t kSOIMarker[2] = {0xff, 0xd8};
  if (!JPEGWrite(out, kSOIMarker, sizeof(kSOIMarker))) {
    return false;
  }
  int dht_index = 0;
  int dqt_index = 0;
  int app_index = 0;
  int com_index = 0;
  int data_index = 0;
  int scan_index = 0;
  std::vector<HuffmanCodeTable> dc_huff_table(kMaxHuffmanTables);
  std::vector<HuffmanCodeTable> ac_huff_table(kMaxHuffmanTables);

  const std::vector<JPEGScanInfo>& scan_info = jpg.scan_info;
  const std::vector<JPEGHuffmanCode>& huffman_code = jpg.huffman_code;
  const int* pad_bits = nullptr;
  const int* pad_bits_end = nullptr;
  if (jpg.has_zero_padding_bit) {
    pad_bits = jpg.padding_bits.data();
    pad_bits_end = pad_bits + jpg.padding_bits.size();
  }
  bool seen_dri_marker = false;
  // progressive/sequential will be decided based on SOF marker
  bool is_progressive = false;

  for (size_t i = 0; i < marker_order.size(); ++i) {
    uint8_t marker = marker_order[i];
    bool ok = true;
    switch (marker) {
      case 0xc0:
      case 0xc1:
      case 0xc2:
        is_progressive = (marker == 0xc2);
        ok = EncodeSOF(jpg, marker, out);
        break;
      case 0xc9:
      case 0xca:
        ok = EncodeSOF(jpg, marker, out);
        break;
      case 0xc4:
        ok = EncodeDHT(huffman_code, &dht_index, out, &dc_huff_table,
                       &ac_huff_table);
        break;
      case 0xd0:
      case 0xd1:
      case 0xd2:
      case 0xd3:
      case 0xd4:
      case 0xd5:
      case 0xd6:
      case 0xd7: {
        uint8_t marker_data[2] = {0xff, marker};
        ok = JPEGWrite(out, marker_data, sizeof(marker_data));
      } break;
      case 0xd9:
        // Found end marker.
        break;
      case 0xda:
        ok = EncodeScan(jpg, scan_info[scan_index++], is_progressive,
                        dc_huff_table, ac_huff_table,
                        seen_dri_marker ? jpg.restart_interval : 0, &pad_bits,
                        pad_bits_end, out);
        break;
      case 0xdb:
        ok = EncodeDQT(jpg, &dqt_index, out);
        break;
      case 0xdd:
        seen_dri_marker = true;
        ok = EncodeDRI(jpg.restart_interval, out);
        break;
      case 0xe0:
      case 0xe1:
      case 0xe2:
      case 0xe3:
      case 0xe4:
      case 0xe5:
      case 0xe6:
      case 0xe7:
      case 0xe8:
      case 0xe9:
      case 0xea:
      case 0xeb:
      case 0xec:
      case 0xed:
      case 0xee:
      case 0xef:
        // TODO: check that marker corresponds to payload?
        ok = EncodeAPP(jpg, app_index++, out);
        break;
      case 0xfe:
        ok = EncodeCOM(jpg, com_index++, out);
        break;
      case 0xff:
        ok = EncodeInterMarkerData(jpg, data_index++, out);
        break;
      default:
        ok = false;
        break;
    }
    if (!ok) {
      BRUNSLI_LOG_DEBUG() << "Failed to encode marker " << std::hex << marker
                          << BRUNSLI_ENDL();
      return false;
    }
  }
  static const uint8_t kEOIMarker[2] = {0xff, 0xd9};
  return (JPEGWrite(out, kEOIMarker, sizeof(kEOIMarker)) &&
          JPEGWrite(out, jpg.tail_data));
}

}  // namespace brunsli
