// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <brunsli/jpeg_data_reader.h>

#include <algorithm>
#include <string>
#include <vector>

#include "../common/constants.h"
#include <brunsli/jpeg_data.h>
#include "../common/platform.h"
#include <brunsli/types.h>
#include "./jpeg_huffman_decode.h"

namespace brunsli {

namespace {

// Macros for commonly used error conditions.

#define VERIFY_LEN(n)                                                          \
  if (*pos + (n) > len) {                                                      \
    BRUNSLI_LOG_INFO() << "Unexpected end of input:"                           \
                       << " pos=" << *pos << " need=" << (n) << " len=" << len \
                       << BRUNSLI_ENDL();                                      \
    jpg->error = JPEGReadError::UNEXPECTED_EOF;                                \
    return false;                                                              \
  }

#define VERIFY_INPUT(var, low, high, code)                                     \
  if (var < low || var > high) {                                               \
    BRUNSLI_LOG_INFO() << "Invalid " << #var << ": " << var << BRUNSLI_ENDL(); \
    jpg->error = JPEGReadError::INVALID_##code;                                \
    return false;                                                              \
  }

#define VERIFY_MARKER_END()                                                   \
  if (start_pos + marker_len != *pos) {                                       \
    BRUNSLI_LOG_INFO() << "Invalid marker length:"                            \
                       << " declared=" << marker_len                          \
                       << " actual=" << (*pos - start_pos) << BRUNSLI_ENDL(); \
    jpg->error = JPEGReadError::WRONG_MARKER_SIZE;                            \
    return false;                                                             \
  }

#define EXPECT_MARKER()                                                       \
  if (pos + 2 > len || data[pos] != 0xff) {                                   \
    BRUNSLI_LOG_INFO() << "Marker byte (0xff) expected,"                      \
                       << " found: " << (pos < len ? data[pos] : 0)           \
                       << " pos=" << pos << " len=" << len << BRUNSLI_ENDL(); \
    jpg->error = JPEGReadError::MARKER_BYTE_NOT_FOUND;                        \
    return false;                                                             \
  }

// Returns ceil(a/b).
inline int DivCeil(int a, int b) {
  return (a + b - 1) / b;
}

inline int ReadUint8(const uint8_t* data, size_t* pos) {
  return data[(*pos)++];
}

inline int ReadUint16(const uint8_t* data, size_t* pos) {
  int v = (data[*pos] << 8) + data[*pos + 1];
  *pos += 2;
  return v;
}

// Reads the Start of Frame (SOF) marker segment and fills in *jpg with the
// parsed data.
bool ProcessSOF(const uint8_t* data, const size_t len,
                JpegReadMode mode, size_t* pos, JPEGData* jpg) {
  if (jpg->width != 0) {
    BRUNSLI_LOG_INFO() << "Duplicate SOF marker." << BRUNSLI_ENDL();
    jpg->error = JPEGReadError::DUPLICATE_SOF;
    return false;
  }
  const size_t start_pos = *pos;
  VERIFY_LEN(8);
  size_t marker_len = ReadUint16(data, pos);
  int precision = ReadUint8(data, pos);
  int height = ReadUint16(data, pos);
  int width = ReadUint16(data, pos);
  int num_components = ReadUint8(data, pos);
  VERIFY_INPUT(precision, 8, 8, PRECISION);
  VERIFY_INPUT(height, 1, kMaxDimPixels, HEIGHT);
  VERIFY_INPUT(width, 1, kMaxDimPixels, WIDTH);
  VERIFY_INPUT(num_components, 1, kMaxComponents, NUMCOMP);
  VERIFY_LEN(3 * num_components);
  jpg->height = height;
  jpg->width = width;
  jpg->components.resize(num_components);

  // Read sampling factors and quant table index for each component.
  std::vector<bool> ids_seen(256, false);
  for (int i = 0; i < jpg->components.size(); ++i) {
    const int id = ReadUint8(data, pos);
    if (ids_seen[id]) {  // (cf. section B.2.2, syntax of Ci)
      BRUNSLI_LOG_INFO() << "Duplicate ID " << id << " in SOF."
                         << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::DUPLICATE_COMPONENT_ID;
      return false;
    }
    ids_seen[id] = true;
    jpg->components[i].id = id;
    int factor = ReadUint8(data, pos);
    int h_samp_factor = factor >> 4;
    int v_samp_factor = factor & 0xf;
    VERIFY_INPUT(h_samp_factor, 1, kBrunsliMaxSampling, SAMP_FACTOR);
    VERIFY_INPUT(v_samp_factor, 1, kBrunsliMaxSampling, SAMP_FACTOR);
    jpg->components[i].h_samp_factor = h_samp_factor;
    jpg->components[i].v_samp_factor = v_samp_factor;
    jpg->components[i].quant_idx = ReadUint8(data, pos);
    jpg->components[i].max_block_index.resize(v_samp_factor);
    jpg->max_h_samp_factor = std::max(jpg->max_h_samp_factor, h_samp_factor);
    jpg->max_v_samp_factor = std::max(jpg->max_v_samp_factor, v_samp_factor);
  }

  // We have checked above that none of the sampling factors are 0, so the max
  // sampling factors can not be 0.
  jpg->MCU_rows = DivCeil(jpg->height, jpg->max_v_samp_factor * 8);
  jpg->MCU_cols = DivCeil(jpg->width, jpg->max_h_samp_factor * 8);
  // Compute the block dimensions for each component.
  for (int i = 0; i < jpg->components.size(); ++i) {
    JPEGComponent* c = &jpg->components[i];
    if (jpg->max_h_samp_factor % c->h_samp_factor != 0 ||
        jpg->max_v_samp_factor % c->v_samp_factor != 0) {
      BRUNSLI_LOG_INFO() << "Non-integral subsampling ratios."
                         << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::INVALID_SAMPLING_FACTORS;
      return false;
    }
    c->width_in_blocks = jpg->MCU_cols * c->h_samp_factor;
    c->height_in_blocks = jpg->MCU_rows * c->v_samp_factor;
    const uint64_t num_blocks =
        static_cast<uint64_t>(c->width_in_blocks) * c->height_in_blocks;
    if (num_blocks > kBrunsliMaxNumBlocks) {
      BRUNSLI_LOG_INFO() << "Image too large." << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::IMAGE_TOO_LARGE;
      return false;
    }
    c->num_blocks = static_cast<int>(num_blocks);
    if (mode == JPEG_READ_ALL) {
      c->coeffs.resize(c->num_blocks * kDCTBlockSize);
    }
  }
  VERIFY_MARKER_END();
  return true;
}

// Reads the Start of Scan (SOS) marker segment and fills in *scan_info with the
// parsed data.
bool ProcessSOS(const uint8_t* data, const size_t len, size_t* pos,
                JPEGData* jpg) {
  const size_t start_pos = *pos;
  VERIFY_LEN(3);
  size_t marker_len = ReadUint16(data, pos);
  int comps_in_scan = ReadUint8(data, pos);
  VERIFY_INPUT(comps_in_scan, 1, jpg->components.size(), COMPS_IN_SCAN);

  JPEGScanInfo scan_info;
  scan_info.components.resize(comps_in_scan);
  VERIFY_LEN(2 * comps_in_scan);
  std::vector<bool> ids_seen(256, false);
  for (int i = 0; i < comps_in_scan; ++i) {
    int id = ReadUint8(data, pos);
    if (ids_seen[id]) {  // (cf. section B.2.3, regarding CSj)
      BRUNSLI_LOG_INFO() << "Duplicate ID " << id << " in SOS."
                         << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::DUPLICATE_COMPONENT_ID;
      return false;
    }
    ids_seen[id] = true;
    bool found_index = false;
    for (int j = 0; j < jpg->components.size(); ++j) {
      if (jpg->components[j].id == id) {
        scan_info.components[i].comp_idx = j;
        found_index = true;
      }
    }
    if (!found_index) {
      BRUNSLI_LOG_INFO() << "SOS marker: Could not find component with id "
                         << id << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::COMPONENT_NOT_FOUND;
      return false;
    }
    int c = ReadUint8(data, pos);
    int dc_tbl_idx = c >> 4;
    int ac_tbl_idx = c & 0xf;
    VERIFY_INPUT(dc_tbl_idx, 0, 3, HUFFMAN_INDEX);
    VERIFY_INPUT(ac_tbl_idx, 0, 3, HUFFMAN_INDEX);
    scan_info.components[i].dc_tbl_idx = dc_tbl_idx;
    scan_info.components[i].ac_tbl_idx = ac_tbl_idx;
  }
  VERIFY_LEN(3);
  scan_info.Ss = ReadUint8(data, pos);
  scan_info.Se = ReadUint8(data, pos);
  VERIFY_INPUT(scan_info.Ss, 0, 63, START_OF_SCAN);
  VERIFY_INPUT(scan_info.Se, scan_info.Ss, 63, END_OF_SCAN);
  int c = ReadUint8(data, pos);
  scan_info.Ah = c >> 4;
  scan_info.Al = c & 0xf;
  if (scan_info.Ah != 0 && scan_info.Al != scan_info.Ah - 1) {
    // section G.1.1.1.2 : Successive approximation control only improves
    // by one bit at a time. But it's not always respected, so we just issue
    // a warning.
    BRUNSLI_LOG_WARNING() << "Invalid progressive parameters: "
                          << " Al = " << scan_info.Al
                          << " Ah = " << scan_info.Ah << BRUNSLI_ENDL();
  }
  // Check that all the Huffman tables needed for this scan are defined.
  for (int i = 0; i < comps_in_scan; ++i) {
    bool found_dc_table = false;
    bool found_ac_table = false;
    for (int j = 0; j < jpg->huffman_code.size(); ++j) {
      int slot_id = jpg->huffman_code[j].slot_id;
      if (slot_id == scan_info.components[i].dc_tbl_idx) {
        found_dc_table = true;
      } else if (slot_id == scan_info.components[i].ac_tbl_idx + 16) {
        found_ac_table = true;
      }
    }
    if (scan_info.Ss == 0 && !found_dc_table) {
      BRUNSLI_LOG_INFO() << "SOS marker: Could not find DC Huffman table with"
                         << " index " << scan_info.components[i].dc_tbl_idx
                         << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::HUFFMAN_TABLE_NOT_FOUND;
      return false;
    }
    if (scan_info.Se > 0 && !found_ac_table) {
      BRUNSLI_LOG_INFO() << "SOS marker: Could not find AC Huffman table with"
                         << " index " << scan_info.components[i].ac_tbl_idx
                         << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::HUFFMAN_TABLE_NOT_FOUND;
      return false;
    }
  }
  jpg->scan_info.push_back(scan_info);
  VERIFY_MARKER_END();
  return true;
}

// Reads the Define Huffman Table (DHT) marker segment and fills in *jpg with
// the parsed data. Builds the Huffman decoding table in either dc_huff_lut or
// ac_huff_lut, depending on the type and solt_id of Huffman code being read.
bool ProcessDHT(const uint8_t* data, const size_t len,
                JpegReadMode mode,
                std::vector<HuffmanTableEntry>* dc_huff_lut,
                std::vector<HuffmanTableEntry>* ac_huff_lut,
                size_t* pos,
                JPEGData* jpg) {
  const size_t start_pos = *pos;
  VERIFY_LEN(2);
  size_t marker_len = ReadUint16(data, pos);
  if (marker_len == 2) {
    BRUNSLI_LOG_INFO() << "DHT marker: no Huffman table found"
                       << BRUNSLI_ENDL();
    jpg->error = JPEGReadError::EMPTY_DHT;
    return false;
  }
  while (*pos < start_pos + marker_len) {
    VERIFY_LEN(1 + kJpegHuffmanMaxBitLength);
    JPEGHuffmanCode huff;
    huff.slot_id = ReadUint8(data, pos);
    int huffman_index = huff.slot_id;
    int is_ac_table = (huff.slot_id & 0x10) != 0;
    HuffmanTableEntry* huff_lut;
    if (is_ac_table) {
      huffman_index -= 0x10;
      VERIFY_INPUT(huffman_index, 0, 3, HUFFMAN_INDEX);
      huff_lut = &(*ac_huff_lut)[huffman_index * kJpegHuffmanLutSize];
    } else {
      VERIFY_INPUT(huffman_index, 0, 3, HUFFMAN_INDEX);
      huff_lut = &(*dc_huff_lut)[huffman_index * kJpegHuffmanLutSize];
    }
    huff.counts[0] = 0;
    int total_count = 0;
    int space = 1 << kJpegHuffmanMaxBitLength;
    int max_depth = 1;
    for (int i = 1; i <= kJpegHuffmanMaxBitLength; ++i) {
      int count = ReadUint8(data, pos);
      if (count != 0) {
        max_depth = i;
      }
      huff.counts[i] = count;
      total_count += count;
      space -= count * (1 << (kJpegHuffmanMaxBitLength - i));
    }
    if (is_ac_table) {
      VERIFY_INPUT(total_count, 0, kJpegHuffmanAlphabetSize, HUFFMAN_CODE);
    } else {
      VERIFY_INPUT(total_count, 0, kJpegDCAlphabetSize, HUFFMAN_CODE);
    }
    VERIFY_LEN(total_count);
    std::vector<bool> values_seen(256, false);
    for (int i = 0; i < total_count; ++i) {
      uint8_t value = ReadUint8(data, pos);
      if (!is_ac_table) {
        VERIFY_INPUT(value, 0, kJpegDCAlphabetSize - 1, HUFFMAN_CODE);
      }
      if (values_seen[value]) {
        BRUNSLI_LOG_INFO() << "Duplicate Huffman code value " << value
                           << BRUNSLI_ENDL();
        jpg->error = JPEGReadError::INVALID_HUFFMAN_CODE;
        return false;
      }
      values_seen[value] = true;
      huff.values[i] = value;
    }
    // Add an invalid symbol that will have the all 1 code.
    ++huff.counts[max_depth];
    huff.values[total_count] = kJpegHuffmanAlphabetSize;
    space -= (1 << (kJpegHuffmanMaxBitLength - max_depth));
    if (space < 0) {
      BRUNSLI_LOG_INFO() << "Invalid Huffman code lengths." << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::INVALID_HUFFMAN_CODE;
      return false;
    } else if (space > 0 && huff_lut[0].value != 0xffff) {
      // Re-initialize the values to an invalid symbol so that we can recognize
      // it when reading the bit stream using a Huffman code with space > 0.
      for (int i = 0; i < kJpegHuffmanLutSize; ++i) {
        huff_lut[i].bits = 0;
        huff_lut[i].value = 0xffff;
      }
    }
    huff.is_last = (*pos == start_pos + marker_len);
    if (mode == JPEG_READ_ALL) {
      BuildJpegHuffmanTable(&huff.counts[0], &huff.values[0], huff_lut);
    }
    jpg->huffman_code.push_back(huff);
  }
  VERIFY_MARKER_END();
  return true;
}

// Reads the Define Quantization Table (DQT) marker segment and fills in *jpg
// with the parsed data.
bool ProcessDQT(const uint8_t* data, const size_t len, size_t* pos,
                JPEGData* jpg) {
  const size_t start_pos = *pos;
  VERIFY_LEN(2);
  size_t marker_len = ReadUint16(data, pos);
  if (marker_len == 2) {
    BRUNSLI_LOG_INFO() << "DQT marker: no quantization table found"
                       << BRUNSLI_ENDL();
    jpg->error = JPEGReadError::EMPTY_DQT;
    return false;
  }
  while (*pos < start_pos + marker_len && jpg->quant.size() < kMaxQuantTables) {
    VERIFY_LEN(1);
    int quant_table_index = ReadUint8(data, pos);
    int quant_table_precision = quant_table_index >> 4;
    VERIFY_INPUT(quant_table_precision, 0, 1, QUANT_TBL_PRECISION);
    quant_table_index &= 0xf;
    VERIFY_INPUT(quant_table_index, 0, 3, QUANT_TBL_INDEX);
    VERIFY_LEN((quant_table_precision + 1) * kDCTBlockSize);
    JPEGQuantTable table;
    table.index = quant_table_index;
    table.precision = quant_table_precision;
    for (int i = 0; i < kDCTBlockSize; ++i) {
      int quant_val = quant_table_precision ?
          ReadUint16(data, pos) :
          ReadUint8(data, pos);
      VERIFY_INPUT(quant_val, 1, 65535, QUANT_VAL);
      table.values[kJPEGNaturalOrder[i]] = quant_val;
    }
    table.is_last = (*pos == start_pos + marker_len);
    jpg->quant.push_back(table);
  }
  VERIFY_MARKER_END();
  return true;
}

// Reads the DRI marker and saves the restart interval into *jpg.
bool ProcessDRI(const uint8_t* data, const size_t len, size_t* pos,
                bool* found_dri, JPEGData* jpg) {
  if (*found_dri) {
    BRUNSLI_LOG_INFO() << "Duplicate DRI marker." << BRUNSLI_ENDL();
    jpg->error = JPEGReadError::DUPLICATE_DRI;
    return false;
  }
  *found_dri = true;
  const size_t start_pos = *pos;
  VERIFY_LEN(4);
  size_t marker_len = ReadUint16(data, pos);
  int restart_interval = ReadUint16(data, pos);
  jpg->restart_interval = restart_interval;
  VERIFY_MARKER_END();
  return true;
}

// Saves the APP marker segment as a string to *jpg.
bool ProcessAPP(const uint8_t* data, const size_t len, size_t* pos,
                JPEGData* jpg) {
  VERIFY_LEN(2);
  size_t marker_len = ReadUint16(data, pos);
  VERIFY_INPUT(marker_len, 2, 65535, MARKER_LEN);
  VERIFY_LEN(marker_len - 2);
  // Save the marker type together with the app data.
  std::string app_str(reinterpret_cast<const char*>(&data[*pos - 3]),
                      marker_len + 1);
  *pos += marker_len - 2;
  jpg->app_data.push_back(app_str);
  return true;
}

// Saves the COM marker segment as a string to *jpg.
bool ProcessCOM(const uint8_t* data, const size_t len, size_t* pos,
                JPEGData* jpg) {
  VERIFY_LEN(2);
  size_t marker_len = ReadUint16(data, pos);
  VERIFY_INPUT(marker_len, 2, 65535, MARKER_LEN);
  VERIFY_LEN(marker_len - 2);
  std::string com_str(reinterpret_cast<const char*>(&data[*pos - 2]),
                      marker_len);
  *pos += marker_len - 2;
  jpg->com_data.push_back(com_str);
  return true;
}

// Helper structure to read bits from the entropy coded data segment.
struct BitReaderState {
  BitReaderState(const uint8_t* data, const size_t len, size_t pos)
      : data_(data), len_(len) {
    Reset(pos);
  }

  void Reset(size_t pos) {
    pos_ = pos;
    val_ = 0;
    bits_left_ = 0;
    next_marker_pos_ = len_ - 2;
    FillBitWindow();
  }

  // Returns the next byte and skips the 0xff/0x00 escape sequences.
  uint8_t GetNextByte() {
    if (pos_ >= next_marker_pos_) {
      ++pos_;
      return 0;
    }
    uint8_t c = data_[pos_++];
    if (c == 0xff) {
      uint8_t escape = data_[pos_];
      if (escape == 0) {
        ++pos_;
      } else {
        // 0xff was followed by a non-zero byte, which means that we found the
        // start of the next marker segment.
        next_marker_pos_ = pos_ - 1;
      }
    }
    return c;
  }

  void FillBitWindow() {
    if (bits_left_ <= 16) {
      while (bits_left_ <= 56) {
        val_ <<= 8;
        val_ |= (uint64_t)GetNextByte();
        bits_left_ += 8;
      }
    }
  }

  int ReadBits(int nbits) {
    FillBitWindow();
    uint64_t val = (val_ >> (bits_left_ - nbits)) & ((1ULL << nbits) - 1);
    bits_left_ -= nbits;
    return val;
  }

  // Sets *pos to the next stream position where parsing should continue.
  // Enqueue the padding bits seen (0 or 1).
  // Returns false if there is inconsistent or invalid padding or the stream
  // ended too early.
  bool FinishStream(JPEGData* jpg, size_t* pos) {
    int npadbits = bits_left_ & 7;
    if (npadbits > 0) {
      uint64_t padmask = (1ULL << npadbits) - 1;
      uint64_t padbits = (val_ >> (bits_left_ - npadbits)) & padmask;
      if (padbits != padmask) {
        jpg->has_zero_padding_bit = true;
      }
      for (int i = npadbits - 1; i >= 0; --i) {
        jpg->padding_bits.push_back((padbits >> i) & 1);
      }
    }
    // Give back some bytes that we did not use.
    int unused_bytes_left = bits_left_ >> 3;
    while (unused_bytes_left-- > 0) {
      --pos_;
      // If we give back a 0 byte, we need to check if it was a 0xff/0x00 escape
      // sequence, and if yes, we need to give back one more byte.
      if (pos_ < next_marker_pos_ &&
          data_[pos_] == 0 && data_[pos_ - 1] == 0xff) {
        --pos_;
      }
    }
    if (pos_ > next_marker_pos_) {
      // Data ran out before the scan was complete.
      BRUNSLI_LOG_INFO() << "Unexpected end of scan." << BRUNSLI_ENDL();
      return false;
    }
    *pos = pos_;
    return true;
  }

  // get truncated block position in bs
  bool get_pos(JPEGData* jpg, size_t* pos) {
    // Give back some bytes that we did not use.
    int unused_bytes_left = bits_left_ >> 3;
    uint64_t npos = pos_;
    while (unused_bytes_left-- > 0) {
      --npos;
      // If we give back a 0 byte, we need to check if it was a 0xff/0x00 escape
      // sequence, and if yes, we need to give back one more byte.
      if (npos < next_marker_pos_ &&
          data_[npos] == 0 && data_[npos - 1] == 0xff) {
        --npos;
      }
    }
    if (npos > next_marker_pos_) {
      // Data ran out before the scan was complete.
      //BRUNSLI_LOG_INFO() << "get pos Unexpected end of scan." << BRUNSLI_ENDL();
      return false;
    }
    *pos = npos;
    return true;
  }

  const uint8_t* data_;
  const size_t len_;
  size_t pos_;
  uint64_t val_;
  int bits_left_;
  size_t next_marker_pos_;
};

// Returns the next Huffman-coded symbol.
int ReadSymbol(const HuffmanTableEntry* table, BitReaderState* br) {
  int nbits;
  br->FillBitWindow();
  int val = (br->val_ >> (br->bits_left_ - 8)) & 0xff;
  table += val;
  nbits = table->bits - 8;
  if (nbits > 0) {
    br->bits_left_ -= 8;
    table += table->value;
    val = (br->val_ >> (br->bits_left_ - nbits)) & ((1 << nbits) - 1);
    table += val;
  }
  br->bits_left_ -= table->bits;
  return table->value;
}

// Returns the DC diff or AC value for extra bits value x and prefix code s.
// See Tables F.1 and F.2 of the spec.
int HuffExtend(int x, int s) {
  return (x < (1 << (s - 1)) ? x + ((-1) << s ) + 1 : x);
}

// Decodes one 8x8 block of DCT coefficients from the bit stream.
bool DecodeDCTBlock(const HuffmanTableEntry* dc_huff,
                    const HuffmanTableEntry* ac_huff,
                    int Ss, int Se, int Al,
                    int* eobrun,
                    bool* reset_state,
                    int* num_zero_runs,
                    BitReaderState* br,
                    JPEGData* jpg,
                    coeff_t* last_dc_coeff,
                    coeff_t* coeffs) {
  int s;
  int r;
  bool eobrun_allowed = Ss > 0;
  if (Ss == 0) {
    s = ReadSymbol(dc_huff, br);
    if (s >= kJpegDCAlphabetSize) {
      BRUNSLI_LOG_INFO() << "Invalid Huffman symbol " << s
                         << " for DC coefficient." << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::INVALID_SYMBOL;
      return false;
    }
    if (s > 0) {
      r = br->ReadBits(s);
      s = HuffExtend(r, s);
    }
    s += *last_dc_coeff;
    const int dc_coeff = s << Al;
    coeffs[0] = dc_coeff;
    if (dc_coeff != coeffs[0]) {
      BRUNSLI_LOG_INFO() << "Invalid DC coefficient " << dc_coeff
                         << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::NON_REPRESENTABLE_DC_COEFF;
      return false;
    }
    *last_dc_coeff = s;
    ++Ss;
  }
  if (Ss > Se) {
    return true;
  }
  if (*eobrun > 0) {
    --(*eobrun);
    return true;
  }
  *num_zero_runs = 0;
  for (int k = Ss; k <= Se; k++) {
    s = ReadSymbol(ac_huff, br);
    if (s >= kJpegHuffmanAlphabetSize) {
      BRUNSLI_LOG_INFO() << "Invalid Huffman symbol " << s
                         << " for AC coefficient " << k << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::INVALID_SYMBOL;
      return false;
    }
    r = s >> 4;
    s &= 15;
    if (s > 0) {
      k += r;
      if (k > Se) {
        BRUNSLI_LOG_INFO() << "Out-of-band coefficient " << k << " band was "
                           << Ss << "-" << Se << BRUNSLI_ENDL();
        jpg->error = JPEGReadError::OUT_OF_BAND_COEFF;
        return false;
      }
      if (s + Al >= kJpegDCAlphabetSize) {
        BRUNSLI_LOG_INFO() << "Out of range AC coefficient value: s = " << s
                           << " Al = " << Al << " k = " << k << BRUNSLI_ENDL();
        jpg->error = JPEGReadError::NON_REPRESENTABLE_AC_COEFF;
        return false;
      }
      r = br->ReadBits(s);
      s = HuffExtend(r, s);
      coeffs[kJPEGNaturalOrder[k]] = s << Al;
      *num_zero_runs = 0;
    } else if (r == 15) {
      k += 15;
      ++(*num_zero_runs);
    } else {
      if (eobrun_allowed && k == Ss && *eobrun == 0) {
        // We have two end-of-block runs right after each other, so we signal
        // the jpeg encoder to force a state reset at this point.
        *reset_state = true;
      }
      *eobrun = 1 << r;
      if (r > 0) {
        if (!eobrun_allowed) {
          BRUNSLI_LOG_INFO()
              << "End-of-block run crossing DC coeff." << BRUNSLI_ENDL();
          jpg->error = JPEGReadError::EOB_RUN_TOO_LONG;
          return false;
        }
        *eobrun += br->ReadBits(r);
      }
      break;
    }
  }
  --(*eobrun);
  return true;
}

bool RefineDCTBlock(const HuffmanTableEntry* ac_huff,
                    int Ss, int Se, int Al,
                    int* eobrun,
                    bool* reset_state,
                    BitReaderState* br,
                    JPEGData* jpg,
                    coeff_t* coeffs) {
  bool eobrun_allowed = Ss > 0;
  if (Ss == 0) {
    int s = br->ReadBits(1);
    coeff_t dc_coeff = coeffs[0];
    dc_coeff |= s << Al;
    coeffs[0] = dc_coeff;
    ++Ss;
  }
  if (Ss > Se) {
    return true;
  }
  int p1 = 1 << Al;
  int m1 = (-1) << Al;
  int k = Ss;
  int r;
  int s;
  bool in_zero_run = false;
  if (*eobrun <= 0) {
    for (; k <= Se; k++) {
      s = ReadSymbol(ac_huff, br);
      if (s >= kJpegHuffmanAlphabetSize) {
        BRUNSLI_LOG_INFO() << "Invalid Huffman symbol " << s
                           << " for AC coefficient " << k << BRUNSLI_ENDL();
        jpg->error = JPEGReadError::INVALID_SYMBOL;
        return false;
      }
      r = s >> 4;
      s &= 15;
      if (s) {
        if (s != 1) {
          BRUNSLI_LOG_INFO() << "Invalid Huffman symbol " << s
                             << " for AC coefficient " << k << BRUNSLI_ENDL();
          jpg->error = JPEGReadError::INVALID_SYMBOL;
          return false;
        }
        s = br->ReadBits(1) ? p1 : m1;
        in_zero_run = false;
      } else {
        if (r != 15) {
          if (eobrun_allowed && k == Ss && *eobrun == 0) {
            // We have two end-of-block runs right after each other, so we
            // signal the jpeg encoder to force a state reset at this point.
            *reset_state = true;
          }
          *eobrun = 1 << r;
          if (r > 0) {
            if (!eobrun_allowed) {
              BRUNSLI_LOG_INFO()
                  << "End-of-block run crossing DC coeff." << BRUNSLI_ENDL();
              jpg->error = JPEGReadError::EOB_RUN_TOO_LONG;
              return false;
            }
            *eobrun += br->ReadBits(r);
          }
          break;
        }
        in_zero_run = true;
      }
      do {
        coeff_t thiscoef = coeffs[kJPEGNaturalOrder[k]];
        if (thiscoef != 0) {
          if (br->ReadBits(1)) {
            if ((thiscoef & p1) == 0) {
              if (thiscoef >= 0) {
                thiscoef += p1;
              } else {
                thiscoef += m1;
              }
            }
          }
          coeffs[kJPEGNaturalOrder[k]] = thiscoef;
        } else {
          if (--r < 0) {
            break;
          }
        }
        k++;
      } while (k <= Se);
      if (s) {
        if (k > Se) {
          BRUNSLI_LOG_INFO() << "Out-of-band coefficient " << k << " band was "
                             << Ss << "-" << Se << BRUNSLI_ENDL();
          jpg->error = JPEGReadError::OUT_OF_BAND_COEFF;
          return false;
        }
        coeffs[kJPEGNaturalOrder[k]] = s;
      }
    }
  }
  if (in_zero_run) {
    BRUNSLI_LOG_INFO() << "Extra zero run before end-of-block."
                       << BRUNSLI_ENDL();
    jpg->error = JPEGReadError::EXTRA_ZERO_RUN;
    return false;
  }
  if (*eobrun > 0) {
    for (; k <= Se; k++) {
      coeff_t thiscoef = coeffs[kJPEGNaturalOrder[k]];
      if (thiscoef != 0) {
        if (br->ReadBits(1)) {
          if ((thiscoef & p1) == 0) {
            if (thiscoef >= 0) {
              thiscoef += p1;
            } else {
              thiscoef += m1;
            }
          }
        }
        coeffs[kJPEGNaturalOrder[k]] = thiscoef;
      }
    }
  }
  --(*eobrun);
  return true;
}

bool ProcessRestart(const uint8_t* data, const size_t len,
                    int* next_restart_marker, BitReaderState* br,
                    JPEGData* jpg) {
  size_t pos = 0;
  if (!br->FinishStream(jpg, &pos)) {
    jpg->error = JPEGReadError::INVALID_SCAN;
    return false;
  }
  int expected_marker = 0xd0 + *next_restart_marker;
  EXPECT_MARKER();
  int marker = data[pos + 1];
  if (marker != expected_marker) {
    BRUNSLI_LOG_INFO() << "Did not find expected restart"
                       << " marker " << expected_marker << " actual=" << marker
                       << BRUNSLI_ENDL();
    jpg->error = JPEGReadError::WRONG_RESTART_MARKER;
    return false;
  }
  br->Reset(pos + 2);
  *next_restart_marker += 1;
  *next_restart_marker &= 0x7;
  return true;
}

bool ProcessScan(const uint8_t* data, const size_t len,
                 const std::vector<HuffmanTableEntry>& dc_huff_lut,
                 const std::vector<HuffmanTableEntry>& ac_huff_lut,
                 uint16_t scan_progression[kMaxComponents][kDCTBlockSize],
                 bool is_progressive,
                 size_t* pos,
                 JPEGData* jpg,
                 bool* is_truncated) {
  jpg->read_scan_numbers++;
  BRUNSLI_LOG_DEBUG() << "------ProcessScan begin------ " << BRUNSLI_ENDL();
  
  if (!ProcessSOS(data, len, pos, jpg)) {
    return false;
  }
  JPEGScanInfo* scan_info = &jpg->scan_info.back();
  bool is_interleaved = (scan_info->components.size() > 1);
  int MCUs_per_row;
  int MCU_rows;
  if (is_interleaved) {
    MCUs_per_row = jpg->MCU_cols;
    MCU_rows = jpg->MCU_rows;
  } else {
    const JPEGComponent& c = jpg->components[scan_info->components[0].comp_idx];
    MCUs_per_row =
        DivCeil(jpg->width * c.h_samp_factor, 8 * jpg->max_h_samp_factor);
    MCU_rows =
        DivCeil(jpg->height * c.v_samp_factor, 8 * jpg->max_v_samp_factor);
  }
  coeff_t last_dc_coeff[kMaxComponents] = {0};
  BitReaderState br(data, len, *pos);
  int restarts_to_go = jpg->restart_interval;
  int next_restart_marker = 0;
  int eobrun = -1;
  int block_scan_index = 0;
  size_t last_mcu_row = 0;
  const int Al = is_progressive ? scan_info->Al : 0;
  const int Ah = is_progressive ? scan_info->Ah : 0;
  const int Ss = is_progressive ? scan_info->Ss : 0;
  const int Se = is_progressive ? scan_info->Se : 63;
  const uint16_t scan_bitmask = Ah == 0 ? (0xffff << Al) : (1u << Al);
  const uint16_t refinement_bitmask = (1 << Al) - 1;
  for (int i = 0; i < scan_info->components.size(); ++i) {
    int comp_idx = scan_info->components[i].comp_idx;
    for (int k = Ss; k <= Se; ++k) {
      if (scan_progression[comp_idx][k] & scan_bitmask) {
        BRUNSLI_LOG_INFO() << "Overlapping scans: component = " << comp_idx
                           << " k = " << k
                           << " prev_mask: " << scan_progression[i][k]
                           << " cur_mask: " << scan_bitmask << BRUNSLI_ENDL();
        jpg->error = JPEGReadError::OVERLAPPING_SCANS;
        return false;
      }
      if (scan_progression[comp_idx][k] & refinement_bitmask) {
        BRUNSLI_LOG_INFO() << "Invalid scan order,"
                           << " a more refined scan was already done:"
                           << " component = " << comp_idx << " k = " << k
                           << " prev_mask: " << scan_progression[i][k]
                           << " cur_mask: " << scan_bitmask << BRUNSLI_ENDL();
        jpg->error = JPEGReadError::INVALID_SCAN_ORDER;
        return false;
      }
      scan_progression[comp_idx][k] |= scan_bitmask;
    }
  }
  if (Al > 10) {
    BRUNSLI_LOG_INFO() << "Scan parameter Al = " << Al
                       << " is not supported in brunsli." << BRUNSLI_ENDL();
    jpg->error = JPEGReadError::NON_REPRESENTABLE_AC_COEFF;
    return false;
  }
  BRUNSLI_LOG_DEBUG() << "MCU_rows: " << MCU_rows << "  MCUs_per_row: " << MCUs_per_row << BRUNSLI_ENDL();
  for (int mcu_y = 0; mcu_y < MCU_rows; ++mcu_y) {
    //for truncated jpeg pos
    if (!br.get_pos(jpg, &jpg->last_mcu_row_pos)) {
      *is_truncated = true;
      BRUNSLI_LOG_DEBUG() << "is_truncated1, the current rows " << mcu_y + 1 << BRUNSLI_ENDL();
      BRUNSLI_LOG_DEBUG() << "is_truncated1, the last rows " << last_mcu_row + 1 << BRUNSLI_ENDL();
      goto READ_END;
    }
    last_mcu_row = mcu_y - 1;

    for (int mcu_x = 0; mcu_x < MCUs_per_row; ++mcu_x) {      
      //for truncated jpeg
      if (!br.get_pos(jpg, &jpg->last_mcu_pos)) {
        *is_truncated = true;
        BRUNSLI_LOG_DEBUG() << "is_truncated2, the current rows " << mcu_y + 1 << BRUNSLI_ENDL();
        BRUNSLI_LOG_DEBUG() << "is_truncated2, the last rows " << last_mcu_row + 1 << BRUNSLI_ENDL();
        goto READ_END;
      }
      // Handle the restart intervals.
      if (jpg->restart_interval > 0) {
        if (restarts_to_go == 0) {
          if (ProcessRestart(data, len,
                             &next_restart_marker, &br, jpg)) {
            restarts_to_go = jpg->restart_interval;
            memset(last_dc_coeff, 0, sizeof(last_dc_coeff));
            if (eobrun > 0) {
              BRUNSLI_LOG_INFO()
                  << "End-of-block run too long." << BRUNSLI_ENDL();
              jpg->error = JPEGReadError::EOB_RUN_TOO_LONG;
              return false;
            }
            eobrun = -1;  // fresh start
          } else {
            //for truncated jpeg
            *is_truncated = true;
            BRUNSLI_LOG_DEBUG() << "is_truncated3, the current rows " << mcu_y + 1 << BRUNSLI_ENDL();
            BRUNSLI_LOG_DEBUG() << "is_truncated3, the last rows " << last_mcu_row + 1 << BRUNSLI_ENDL();
            goto READ_END;
          }
        }
        --restarts_to_go;
      }
      // Decode one MCU.
      for (int i = 0; i < scan_info->components.size(); ++i) {
        JPEGComponentScanInfo* si = &scan_info->components[i];
        JPEGComponent* c = &jpg->components[si->comp_idx];
        const HuffmanTableEntry* dc_lut =
            &dc_huff_lut[si->dc_tbl_idx * kJpegHuffmanLutSize];
        const HuffmanTableEntry* ac_lut =
            &ac_huff_lut[si->ac_tbl_idx * kJpegHuffmanLutSize];
        int nblocks_y = is_interleaved ? c->v_samp_factor : 1;
        int nblocks_x = is_interleaved ? c->h_samp_factor : 1;
        for (int iy = 0; iy < nblocks_y; ++iy) {
          for (int ix = 0; ix < nblocks_x; ++ix) {            
            int block_y = mcu_y * nblocks_y + iy;
            int block_x = mcu_x * nblocks_x + ix;
            int block_idx = block_y * c->width_in_blocks + block_x;
            bool reset_state = false;
            int num_zero_runs = 0;
            coeff_t* coeffs = &c->coeffs[block_idx * kDCTBlockSize];
            if (Ah == 0) {
              if (!DecodeDCTBlock(dc_lut, ac_lut,
                                  Ss, Se, Al,
                                  &eobrun, &reset_state, &num_zero_runs,
                                  &br, jpg, &last_dc_coeff[si->comp_idx],
                                  coeffs)) {
                //for truncated jpeg
                *is_truncated = true;
                BRUNSLI_LOG_DEBUG() << "is_truncated4, the current rows " << mcu_y + 1 << BRUNSLI_ENDL();
                BRUNSLI_LOG_DEBUG() << "is_truncated4, the last rows " << last_mcu_row + 1 << BRUNSLI_ENDL();
                goto READ_END;
              }
            } else {
              if (!RefineDCTBlock(ac_lut, Ss, Se, Al,
                                  &eobrun, &reset_state,
                                  &br, jpg, coeffs)) {
                //for truncated jpeg
                *is_truncated = true;
                BRUNSLI_LOG_DEBUG() << "is_truncated5, the current rows " << mcu_y + 1 << BRUNSLI_ENDL();
                BRUNSLI_LOG_DEBUG() << "is_truncated5, the last rows " << last_mcu_row + 1 << BRUNSLI_ENDL();
                goto READ_END;
              }
            }
            if (reset_state) {
              scan_info->reset_points.insert(block_scan_index);
            }
            if (num_zero_runs > 0) {
              JPEGScanInfo::ExtraZeroRunInfo info;
              info.block_idx = block_scan_index;
              info.num_extra_zero_runs = num_zero_runs;
              scan_info->extra_zero_runs.push_back(info);
            }
            ++block_scan_index;
            // add the max block index for truncated jpeg
            c->max_block_index[iy] ++;
          }
        }
      }
    }
  }

  READ_END:

  for (int i = 0; i < jpg->components.size(); ++i) {
    for (int j = 0; j < jpg->components[i].v_samp_factor; ++j) {
        BRUNSLI_LOG_DEBUG() << "jpg->components[i].max_block_index[j]:" << jpg->components[i].max_block_index[j] << BRUNSLI_ENDL();
    }
  } 

  if (!(*is_truncated)) {
    BRUNSLI_LOG_DEBUG() << "finish mcu row " << last_mcu_row + 2 << " scan numbers " << jpg->read_scan_numbers << BRUNSLI_ENDL();
    if (eobrun > 0) {
      BRUNSLI_LOG_INFO() << "End-of-block run too long." << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::EOB_RUN_TOO_LONG;
      return false;
    }
    if (!br.FinishStream(jpg, pos)) {
      *is_truncated = true;
      *pos = jpg->last_mcu_row_pos;
      for (int i = 0; i < jpg->components.size(); ++i) {
        for (int j = 0; j < jpg->components[i].v_samp_factor; ++j) {
          jpg->components[i].max_block_index[j] = (last_mcu_row + 1) * jpg->components[i].h_samp_factor * MCUs_per_row;
        }
      }
      return true;
    }
    if (*pos > len) {
      BRUNSLI_LOG_INFO() << "Unexpected end of file during scan. pos=" << *pos
                        << " len=" << len << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::UNEXPECTED_EOF;
      return false;
    }
  } else {
    BRUNSLI_LOG_DEBUG() << "last_mcu_row " << last_mcu_row + 1 << " scan numbers " << jpg->read_scan_numbers << BRUNSLI_ENDL();
    *pos = jpg->last_mcu_row_pos;
    for (int i = 0; i < jpg->components.size(); ++i) {
      for (int j = 0; j < jpg->components[i].v_samp_factor; ++j) {
        jpg->components[i].max_block_index[j] = (last_mcu_row + 1) * jpg->components[i].h_samp_factor * MCUs_per_row;
        BRUNSLI_LOG_DEBUG() << "jpg->components[i].max_block_index[j]:" << jpg->components[i].max_block_index[j] << BRUNSLI_ENDL();
      }
    }
  }
  BRUNSLI_LOG_DEBUG() << "------ProcessScan end------ " << BRUNSLI_ENDL();
  return true;
}

// Changes the quant_idx field of the components to refer to the index of the
// quant table in the jpg->quant array.
bool FixupIndexes(JPEGData* jpg) {
  for (int i = 0; i < jpg->components.size(); ++i) {
    JPEGComponent* c = &jpg->components[i];
    bool found_index = false;
    for (int j = 0; j < jpg->quant.size(); ++j) {
      if (jpg->quant[j].index == c->quant_idx) {
        c->quant_idx = j;
        found_index = true;
        break;
      }
    }
    if (!found_index) {
      BRUNSLI_LOG_INFO() << "Quantization table with index " << c->quant_idx
                         << " not found." << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::QUANT_TABLE_NOT_FOUND;
      return false;
    }
  }
  return true;
}

size_t FindNextMarker(const uint8_t* data, const size_t len, size_t pos) {
  // kIsValidMarker[i] == 1 means (0xc0 + i) is a valid marker.
  static const uint8_t kIsValidMarker[] = {
    1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
  };
  size_t num_skipped = 0;
  while (pos + 1 < len &&
         (data[pos] != 0xff || data[pos + 1] < 0xc0 ||
          !kIsValidMarker[data[pos + 1] - 0xc0])) {
    ++pos;
    ++num_skipped;
  }
  return num_skipped;
}

}  // namespace

bool ReadJpeg(const uint8_t* data, const size_t len, JpegReadMode mode,
              JPEGData* jpg) {
  size_t pos = 0;
  // Check SOI marker.
  EXPECT_MARKER();
  int marker = data[pos + 1];
  pos += 2;
  if (marker != 0xd8) {
    BRUNSLI_LOG_INFO() << "Did not find expected SOI marker, actual=" << marker
                       << BRUNSLI_ENDL();
    jpg->error = JPEGReadError::SOI_NOT_FOUND;
    return false;
  }
  int lut_size = kMaxHuffmanTables * kJpegHuffmanLutSize;
  std::vector<HuffmanTableEntry> dc_huff_lut(lut_size);
  std::vector<HuffmanTableEntry> ac_huff_lut(lut_size);
  bool found_sof = false;
  bool found_dri = false;
  uint16_t scan_progression[kMaxComponents][kDCTBlockSize] = { { 0 } };

  jpg->padding_bits.resize(0);
  bool is_progressive = false;  // default
  bool is_truncated = false;    // add for truncated jpeg
  do {
    // Read next marker.
    size_t num_skipped = FindNextMarker(data, len, pos);
    if (num_skipped > 0) {
      // Add a fake marker to indicate arbitrary in-between-markers data.
      jpg->marker_order.push_back(0xff);
      jpg->inter_marker_data.push_back(
          std::string(reinterpret_cast<const char*>(&data[pos]), num_skipped));
      pos += num_skipped;
    }
    EXPECT_MARKER();
    marker = data[pos + 1];
    pos += 2;
    bool ok = true;
    switch (marker) {
      case 0xc0:
      case 0xc1:
      case 0xc2:
        is_progressive = (marker == 0xc2);
        ok = ProcessSOF(data, len, mode, &pos, jpg);
        found_sof = true;
        break;
      case 0xc4:
        ok = ProcessDHT(data, len, mode, &dc_huff_lut, &ac_huff_lut, &pos, jpg);
        break;
      case 0xd0:
      case 0xd1:
      case 0xd2:
      case 0xd3:
      case 0xd4:
      case 0xd5:
      case 0xd6:
      case 0xd7:
        // RST markers do not have any data.
        break;
      case 0xd9:
        // Found end marker.
        break;
      case 0xda:
        if (mode == JPEG_READ_ALL) {
          ok = ProcessScan(data, len, dc_huff_lut, ac_huff_lut,
                           scan_progression, is_progressive, &pos, jpg, &is_truncated);
        }
        break;
      case 0xdb:
        ok = ProcessDQT(data, len, &pos, jpg);
        break;
      case 0xdd:
        ok = ProcessDRI(data, len, &pos, &found_dri, jpg);
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
        if (mode != JPEG_READ_TABLES) {
          ok = ProcessAPP(data, len, &pos, jpg);
        }
        break;
      case 0xfe:
        if (mode != JPEG_READ_TABLES) {
          ok = ProcessCOM(data, len, &pos, jpg);
        }
        break;
      default:
        BRUNSLI_LOG_INFO() << "Unsupported marker: " << marker << " pos=" << pos
                           << " len=" << len << BRUNSLI_ENDL();
        jpg->error = JPEGReadError::UNSUPPORTED_MARKER;
        ok = false;
        break;
    }
    if (!ok) {
      return false;
    }
    jpg->marker_order.push_back(marker);
    if (mode == JPEG_READ_HEADER && found_sof) {
      break;
    }
    if (is_truncated) {
      // dummy EOI marker for decode
      jpg->marker_order.push_back(0xd9);              
      break;
    }
  } while (marker != 0xd9);

  if (!is_truncated && !found_sof) {
    BRUNSLI_LOG_INFO() << "Missing SOF marker." << BRUNSLI_ENDL();
    jpg->error = JPEGReadError::SOF_NOT_FOUND;
    return false;
  }

  if (is_truncated) {
    jpg->is_truncated = true;
  }

  // Supplemental checks.
  if (mode == JPEG_READ_ALL) {
    if (pos < len) {
      if (is_truncated) {
        BRUNSLI_LOG_DEBUG() << "pos: " << jpg->last_mcu_row_pos - kDCTBlockSize << BRUNSLI_ENDL();
        jpg->tail_data.assign(reinterpret_cast<const char*>(&data[jpg->last_mcu_row_pos - kDCTBlockSize]),
                            len - jpg->last_mcu_row_pos + kDCTBlockSize);
      } else {
        jpg->tail_data.assign(reinterpret_cast<const char*>(&data[pos]),
                            len - pos);
      }
    }
    if (!FixupIndexes(jpg)) {
      return false;
    }
    if (jpg->huffman_code.empty()) {
      // Section B.2.4.2: "If a table has never been defined for a particular
      // destination, then when this destination is specified in a scan header,
      // the results are unpredictable."
      BRUNSLI_LOG_INFO() << "Need at least one Huffman code table."
                         << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::HUFFMAN_TABLE_ERROR;
      return false;
    }
    if (jpg->huffman_code.size() >= kMaxDHTMarkers) {
      BRUNSLI_LOG_INFO() << "Too many Huffman tables." << BRUNSLI_ENDL();
      jpg->error = JPEGReadError::HUFFMAN_TABLE_ERROR;
      return false;
    }
  }
  return true;
}

}  // namespace brunsli
