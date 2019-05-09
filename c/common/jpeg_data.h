// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Data structures that represent the contents of a jpeg file.

#ifndef BRUNSLI_COMMON_JPEG_DATA_H_
#define BRUNSLI_COMMON_JPEG_DATA_H_

#include <set>
#include <string>
#include <vector>

#include "./types.h"

namespace brunsli {

static const int kDCTBlockSize = 64;
static const int kMaxComponents = 4;
static const int kMaxQuantTables = 4;
static const int kMaxHuffmanTables = 4;
static const int kJpegHuffmanMaxBitLength = 16;
static const int kJpegHuffmanAlphabetSize = 256;
static const int kJpegDCAlphabetSize = 12;
static const int kMaxDHTMarkers = 512;
static const int kMaxDimPixels = 65535;

static const uint8_t kDefaultQuantMatrix[2][64] = {
  { 16,  11,  10,  16,  24,  40,  51,  61,
    12,  12,  14,  19,  26,  58,  60,  55,
    14,  13,  16,  24,  40,  57,  69,  56,
    14,  17,  22,  29,  51,  87,  80,  62,
    18,  22,  37,  56,  68, 109, 103,  77,
    24,  35,  55,  64,  81, 104, 113,  92,
    49,  64,  78,  87, 103, 121, 120, 101,
    72,  92,  95,  98, 112, 100, 103,  99 },
  { 17,  18,  24,  47,  99,  99,  99,  99,
    18,  21,  26,  66,  99,  99,  99,  99,
    24,  26,  56,  99,  99,  99,  99,  99,
    47,  66,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99 }
};

const int kJPEGNaturalOrder[80] = {
  0,   1,  8, 16,  9,  2,  3, 10,
  17, 24, 32, 25, 18, 11,  4,  5,
  12, 19, 26, 33, 40, 48, 41, 34,
  27, 20, 13,  6,  7, 14, 21, 28,
  35, 42, 49, 56, 57, 50, 43, 36,
  29, 22, 15, 23, 30, 37, 44, 51,
  58, 59, 52, 45, 38, 31, 39, 46,
  53, 60, 61, 54, 47, 55, 62, 63,
  // extra entries for safety in decoder
  63, 63, 63, 63, 63, 63, 63, 63,
  63, 63, 63, 63, 63, 63, 63, 63
};

const int kJPEGZigZagOrder[64] = {
  0,   1,  5,  6, 14, 15, 27, 28,
  2,   4,  7, 13, 16, 26, 29, 42,
  3,   8, 12, 17, 25, 30, 41, 43,
  9,  11, 18, 24, 31, 40, 44, 53,
  10, 19, 23, 32, 39, 45, 52, 54,
  20, 22, 33, 38, 46, 51, 55, 60,
  21, 34, 37, 47, 50, 56, 59, 61,
  35, 36, 48, 49, 57, 58, 62, 63
};

enum struct JPEGReadError {
  OK = 0,
  SOI_NOT_FOUND,
  SOF_NOT_FOUND,
  UNEXPECTED_EOF,
  MARKER_BYTE_NOT_FOUND,
  UNSUPPORTED_MARKER,
  WRONG_MARKER_SIZE,
  INVALID_PRECISION,
  INVALID_WIDTH,
  INVALID_HEIGHT,
  INVALID_NUMCOMP,
  INVALID_SAMP_FACTOR,
  INVALID_START_OF_SCAN,
  INVALID_END_OF_SCAN,
  INVALID_SCAN_BIT_POSITION,
  INVALID_COMPS_IN_SCAN,
  INVALID_HUFFMAN_INDEX,
  INVALID_QUANT_TBL_INDEX,
  INVALID_QUANT_VAL,
  INVALID_MARKER_LEN,
  INVALID_SAMPLING_FACTORS,
  INVALID_HUFFMAN_CODE,
  INVALID_SYMBOL,
  NON_REPRESENTABLE_DC_COEFF,
  NON_REPRESENTABLE_AC_COEFF,
  INVALID_SCAN,
  OVERLAPPING_SCANS,
  INVALID_SCAN_ORDER,
  EXTRA_ZERO_RUN,
  DUPLICATE_DRI,
  DUPLICATE_SOF,
  WRONG_RESTART_MARKER,
  DUPLICATE_COMPONENT_ID,
  COMPONENT_NOT_FOUND,
  HUFFMAN_TABLE_NOT_FOUND,
  HUFFMAN_TABLE_ERROR,
  QUANT_TABLE_NOT_FOUND,
  EMPTY_DHT,
  EMPTY_DQT,
  OUT_OF_BAND_COEFF,
  EOB_RUN_TOO_LONG,
  IMAGE_TOO_LARGE,
};

// Quantization values for an 8x8 pixel block.
struct JPEGQuantTable {
  JPEGQuantTable() : values(kDCTBlockSize), precision(0),
                     index(0), is_last(true) {}

  std::vector<int> values;
  int precision;
  // The index of this quantization table as it was parsed from the input JPEG.
  // Each DQT marker segment contains an 'index' field, and we save this index
  // here. Valid values are 0 to 3.
  int index;
  // Set to true if this table is the last one within its marker segment.
  bool is_last;
};

// Huffman code and decoding lookup table used for DC and AC coefficients.
struct JPEGHuffmanCode {
  JPEGHuffmanCode() : counts(kJpegHuffmanMaxBitLength + 1),
                      values(kJpegHuffmanAlphabetSize + 1),
                      slot_id(0),
                      is_last(true) {}

  // Bit length histogram.
  std::vector<int> counts;
  // Symbol values sorted by increasing bit lengths.
  std::vector<int> values;
  // The index of the Huffman code in the current set of Huffman codes. For AC
  // component Huffman codes, 0x10 is added to the index.
  int slot_id;
  // Set to true if this Huffman code is the last one within its marker segment.
  bool is_last;
};

// Huffman table indexes used for one component of one scan.
struct JPEGComponentScanInfo {
  int comp_idx;
  int dc_tbl_idx;
  int ac_tbl_idx;
};

// Contains information that is used in one scan.
struct JPEGScanInfo {
  // Parameters used for progressive scans (named the same way as in the spec):
  //   Ss : Start of spectral band in zig-zag sequence.
  //   Se : End of spectral band in zig-zag sequence.
  //   Ah : Successive approximation bit position, high.
  //   Al : Successive approximation bit position, low.
  int Ss;
  int Se;
  int Ah;
  int Al;
  std::vector<JPEGComponentScanInfo> components;

  // Extra information required for bit-precise JPEG file reconstruction.

  // Set of block indexes where the jpeg encoder has to flush the end-of-block
  // runs and refinement bits.
  std::set<int> reset_points;
  // The number of extra zero runs (Huffman symbol 0xf0) before the end of
  // block (if nonzero), indexed by block index.
  // All of these symbols can be omitted without changing the pixel values, but
  // some jpeg encoders put these at the end of blocks.
  typedef struct {
    int block_idx;
    int num_extra_zero_runs;
  } ExtraZeroRunInfo;
  std::vector<ExtraZeroRunInfo> extra_zero_runs;
};

typedef int16_t coeff_t;

// Represents one component of a jpeg file.
struct JPEGComponent {
  JPEGComponent() : id(0),
                    h_samp_factor(1),
                    v_samp_factor(1),
                    quant_idx(0),
                    width_in_blocks(0),
                    height_in_blocks(0) {}

  // One-byte id of the component.
  int id;
  // Horizontal and vertical sampling factors.
  // In interleaved mode, each minimal coded unit (MCU) has
  // h_samp_factor x v_samp_factor DCT blocks from this component.
  int h_samp_factor;
  int v_samp_factor;
  // The index of the quantization table used for this component.
  int quant_idx;
  // The dimensions of the component measured in 8x8 blocks.
  int width_in_blocks;
  int height_in_blocks;
  int num_blocks;
  // The DCT coefficients of this component, laid out block-by-block, divided
  // through the quantization matrix values.
  std::vector<coeff_t> coeffs;
};

// Represents a parsed jpeg file.
struct JPEGData {
  JPEGData() : width(0),
               height(0),
               version(0),
               max_h_samp_factor(1),
               max_v_samp_factor(1),
               MCU_rows(0),
               MCU_cols(0),
               restart_interval(0),
               original_jpg(NULL),
               original_jpg_size(0),
               error(JPEGReadError::OK),
               has_zero_padding_bit(false) {}

  int width;
  int height;
  int version;
  int max_h_samp_factor;
  int max_v_samp_factor;
  int MCU_rows;
  int MCU_cols;
  int restart_interval;
  std::vector<std::string> app_data;
  std::vector<std::string> com_data;
  std::vector<JPEGQuantTable> quant;
  std::vector<JPEGHuffmanCode> huffman_code;
  std::vector<JPEGComponent> components;
  std::vector<JPEGScanInfo> scan_info;
  std::vector<uint8_t> marker_order;
  std::vector<std::string> inter_marker_data;
  std::string tail_data;
  const uint8_t* original_jpg;
  size_t original_jpg_size;
  JPEGReadError error;

  // Extra information required for bit-precise JPEG file reconstruction.

  bool has_zero_padding_bit;
  std::vector<int> padding_bits;
};

inline bool JPEGDataIs420(const JPEGData& jpg) {
  return (jpg.components.size() == 3 &&
          jpg.max_h_samp_factor == 2 &&
          jpg.max_v_samp_factor == 2 &&
          jpg.components[0].h_samp_factor == 2 &&
          jpg.components[0].v_samp_factor == 2 &&
          jpg.components[1].h_samp_factor == 1 &&
          jpg.components[1].v_samp_factor == 1 &&
          jpg.components[2].h_samp_factor == 1 &&
          jpg.components[2].v_samp_factor == 1);
}

inline bool JPEGDataIs444(const JPEGData& jpg) {
  return (jpg.components.size() == 3 &&
          jpg.max_h_samp_factor == 1 &&
          jpg.max_v_samp_factor == 1 &&
          jpg.components[0].h_samp_factor == 1 &&
          jpg.components[0].v_samp_factor == 1 &&
          jpg.components[1].h_samp_factor == 1 &&
          jpg.components[1].v_samp_factor == 1 &&
          jpg.components[2].h_samp_factor == 1 &&
          jpg.components[2].v_samp_factor == 1);
}

}  // namespace brunsli

#endif  // BRUNSLI_COMMON_JPEG_DATA_H_
