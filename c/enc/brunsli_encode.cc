// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <brunsli/brunsli_encode.h>

#include <algorithm>
#include <cstdlib>
#include <iterator>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <brotli/encode.h>
#include "../common/constants.h"
#include "../common/context.h"
#include "../common/distributions.h"
#include "../common/lehmer_code.h"
#include "../common/platform.h"
#include "../common/predict.h"
#include "../common/quant_matrix.h"
#include <brunsli/types.h>
#include "./ans_encode.h"
#include "./cluster.h"
#include "./context_map_encode.h"
#include "./histogram_encode.h"
#include <brunsli/jpeg_data_reader.h>
#include "./state.h"
#include "./write_bits.h"

namespace brunsli {

static const int kNumDirectCodes = 8;
static const int kBrotliQuality = 6;
static const int kBrotliWindowBits = 18;

using ::brunsli::internal::enc::BlockI32;
using ::brunsli::internal::enc::ComponentMeta;
using ::brunsli::internal::enc::DataStream;
using ::brunsli::internal::enc::EntropyCodes;
using ::brunsli::internal::enc::EntropySource;
using ::brunsli::internal::enc::Histogram;
using ::brunsli::internal::enc::State;

using ::brunsli::internal::enc::SelectContextBits;

// Returns an upper bound on the encoded size of the jpeg internals section.
size_t EstimateAuxDataSize(const JPEGData& jpg) {
  size_t size = (jpg.marker_order.size() + 272 * jpg.huffman_code.size() +
                 7 * jpg.scan_info.size() + 16);
  for (size_t i = 0; i < jpg.scan_info.size(); ++i) {
    size += 7 * jpg.scan_info[i].reset_points.size();
    size += 7 * jpg.scan_info[i].extra_zero_runs.size();
  }
  size_t nsize = jpg.has_zero_padding_bit ? jpg.padding_bits.size() : 0;
  // We have maximum 4 * 9 bits for describing nsize, plus nsize padding bits,
  // plus maximum 7 bits for going to the next byte boundary.
  size += (nsize + 43) >> 3;
  for (size_t i = 0; i < jpg.inter_marker_data.size(); ++i) {
    size += 5 + jpg.inter_marker_data[i].size();
  }
  return size;
}

size_t GetMaximumBrunsliEncodedSize(const JPEGData& jpg) {
  // Rough estimate is 1.2 * uncompressed size plus some more for the header.
  size_t hdr_size = 1 << 20;
  hdr_size += EstimateAuxDataSize(jpg);
  for (const auto& data : jpg.app_data) {
    hdr_size += data.size();
  }
  for (const auto& data : jpg.com_data) {
    hdr_size += data.size();
  }
  hdr_size += jpg.tail_data.size();
  size_t num_pixels = jpg.width * jpg.height * jpg.components.size();
  // TODO(eustas): are we certain about the multiplier?
  return static_cast<size_t>(num_pixels * 1.2) + hdr_size;
}

size_t Base128Size(size_t val) {
  size_t size = 1;
  for (; val >= 128; val >>= 7) ++size;
  return size;
}

size_t EncodeBase128(size_t val, uint8_t* data) {
  size_t len = 0;
  do {
    data[len++] = (val & 0x7f) | (val >= 128 ? 0x80 : 0);
    val >>= 7;
  } while (val > 0);
  return len;
}

void EncodeBase128Fix(size_t val, size_t len, uint8_t* data) {
  for (size_t i = 0; i < len; ++i) {
    *(data++) = (val & 0x7f) | (i + 1 < len ? 0x80 : 0);
    val >>= 7;
  }
}

bool TransformApp0Marker(const std::vector<uint8_t>& s,
                         std::vector<uint8_t>* out) {
  if (s.size() != 17) return false;
  if (memcmp(s.data(), AppData_0xe0, 9) != 0) return false;
  if ((s[9] == 1 || s[9] == 2) &&  // version / 1.1 or 1.2
      s[10] < 4 &&                 // density units
      s[15] == 0 && s[16] == 0) {  // thumbnail size / no thumbnail
    const uint8_t x_dens_hi = s[11];
    const uint8_t x_dens_lo = s[12];
    int x_dens = (x_dens_hi << 8) + x_dens_lo;
    const uint8_t y_dens_hi = s[13];
    const uint8_t y_dens_lo = s[14];
    int y_dens = (y_dens_hi << 8) + y_dens_lo;
    int density_ix = -1;
    for (int k = 0; k < kMaxApp0Densities; ++k) {
      if (x_dens == kApp0Densities[k] && y_dens == x_dens) {
        density_ix = k;
      }
    }
    if (density_ix >= 0) {
      uint8_t app0_status = (s[9] - 1) | s[10] << 1 | density_ix << 3;
      *out = std::vector<uint8_t>(1);
      out->at(0) = app0_status;
      return true;
    }
  }
  return false;
}

bool TransformApp2Marker(const std::vector<uint8_t>& s,
                         std::vector<uint8_t>* out) {
  if (s.size() == 3161 && !memcmp(s.data(), AppData_0xe2, 84) &&
      !memcmp(s.data() + 85, AppData_0xe2 + 85, 3161 - 85)) {
    std::vector<uint8_t> code(2);
    code[0] = 0x80;
    code[1] = s[84];
    *out = code;
    return true;
  }
  return false;
}

bool TransformApp12Marker(const std::vector<uint8_t>& s,
                          std::vector<uint8_t>* out) {
  if (s.size() == 18 && !memcmp(s.data(), AppData_0xec, 15) &&
      !memcmp(s.data() + 16, AppData_0xec + 16, 18 - 16)) {
    std::vector<uint8_t> code(2);
    code[0] = 0x81;
    code[1] = s[15];
    *out = code;
    return true;
  }
  return false;
}

bool TransformApp14Marker(const std::vector<uint8_t>& s,
                          std::vector<uint8_t>* out) {
  if (s.size() == 15 && !memcmp(&s[0], AppData_0xee, 10) &&
      !memcmp(&s[11], AppData_0xee + 11, 15 - 11)) {
    std::vector<uint8_t> code(2);
    code[0] = 0x82;
    code[1] = s[10];
    *out = code;
    return true;
  }
  return false;
}

std::vector<uint8_t> TransformAppMarker(const std::vector<uint8_t>& s,
                                        size_t* transformed_marker_count) {
  std::vector<uint8_t> out;
  if (TransformApp0Marker(s, &out)) {
    (*transformed_marker_count)++;
    return out;
  }
  if (TransformApp2Marker(s, &out)) {
    (*transformed_marker_count)++;
    return out;
  }
  if (TransformApp12Marker(s, &out)) {
    (*transformed_marker_count)++;
    return out;
  }
  if (TransformApp14Marker(s, &out)) {
    (*transformed_marker_count)++;
    return out;
  }
  return s;
}

int GetQuantTableId(const JPEGQuantTable& q, bool is_chroma,
                    uint8_t dst[kDCTBlockSize]) {
  for (int j = 0; j < kNumStockQuantTables; ++j) {
    bool match_found = true;
    for (int k = 0; match_found && k < kDCTBlockSize; ++k) {
      if (q.values[k] != kStockQuantizationTables[is_chroma][j][k]) {
        match_found = false;
      }
    }
    if (match_found) {
      return j;
    }
  }
  return kNumStockQuantTables + FindBestMatrix(&q.values[0], is_chroma, dst);
}

void EncodeVarint(int n, int max_bits, Storage* storage) {
  int b;
  BRUNSLI_DCHECK(n < (1 << max_bits));
  for (b = 0; n != 0 && b < max_bits; ++b) {
    if (b + 1 != max_bits) {
      WriteBits(1, 1, storage);
    }
    WriteBits(1, n & 1, storage);
    n >>= 1;
  }
  if (b < max_bits) {
    WriteBits(1, 0, storage);
  }
}

// encodes an integer with packets of 'nbits' bits, limited to 'max_symbols'
// emitted symbols.
void EncodeLimitedVarint(size_t bits, int nbits, int max_symbols,
                         Storage* storage) {
  const size_t mask = (static_cast<size_t>(1) << nbits) - 1;
  for (int b = 0; b < max_symbols; ++b) {
    WriteBits(1, bits != 0, storage);
    if (bits == 0) break;
    WriteBits(nbits, bits & mask, storage);
    bits >>= nbits;
  }
}

bool EncodeQuantTables(const JPEGData& jpg, Storage* storage) {
  if (jpg.quant.empty() || jpg.quant.size() > 4) {
    // If ReadJpeg() succeeded with JPEG_READ_ALL mode, this should not happen.
    return false;
  }
  WriteBits(2, jpg.quant.size() - 1, storage);
  for (size_t i = 0; i < jpg.quant.size(); ++i) {
    const JPEGQuantTable& q = jpg.quant[i];
    for (int k = 0; k < kDCTBlockSize; ++k) {
      const int j = kJPEGNaturalOrder[k];
      if (q.values[j] == 0) {
        // Note: ReadJpeg() checks this case and discards such JPEG files.
        return false;
      }
    }

    uint8_t quant_approx[kDCTBlockSize];
    const int code = GetQuantTableId(q, i > 0, quant_approx);
    WriteBits(1, (code >= kNumStockQuantTables), storage);
    if (code < kNumStockQuantTables) {
      WriteBits(3, code, storage);
    } else {
      int q_factor = code - kNumStockQuantTables;
      BRUNSLI_DCHECK(q_factor < kQFactorLimit);
      WriteBits(kQFactorBits, q_factor, storage);
      int last_diff = 0;  // difference predictor
      for (int k = 0; k < kDCTBlockSize; ++k) {
        const int j = kJPEGNaturalOrder[k];
        const int new_diff = q.values[j] - quant_approx[j];
        int diff = new_diff - last_diff;
        last_diff = new_diff;
        WriteBits(1, diff != 0, storage);
        if (diff) {
          WriteBits(1, diff < 0, storage);
          if (diff < 0) diff = -diff;
          diff -= 1;
          // This only happens on 16-bit precision with crazy values,
          // e.g. [..., 65535, 1, 65535,...]
          if (diff > 65535) return false;
          EncodeVarint(diff, 16, storage);
        }
      }
    }
  }
  for (size_t i = 0; i < jpg.components.size(); ++i) {
    WriteBits(2, jpg.components[i].quant_idx, storage);
  }
  return true;
}

bool EncodeHuffmanCode(const JPEGHuffmanCode& huff, bool is_known_last,
                       Storage* storage) {
  WriteBits(2, huff.slot_id & 0xf, storage);
  WriteBits(1, huff.slot_id >> 4, storage);
  if (!is_known_last) {
    WriteBits(1, huff.is_last, storage);
  } else if (!huff.is_last) {
    return false;
  }
  int is_dc_table = (huff.slot_id >> 4) == 0;
  int total_count = 0;
  int space = 1 << kJpegHuffmanMaxBitLength;
  int max_len = kJpegHuffmanMaxBitLength;
  int max_count = is_dc_table ? kJpegDCAlphabetSize : kJpegHuffmanAlphabetSize;
  int found_match = 0;
  int stock_table_idx = 0;
  if (is_dc_table) {
    for (int i = 0; i < kNumStockDCHuffmanCodes && !found_match; ++i) {
      if (memcmp(&huff.counts[1], kStockDCHuffmanCodeCounts[i],
                 sizeof(kStockDCHuffmanCodeCounts[i])) == 0 &&
          memcmp(&huff.values[0], kStockDCHuffmanCodeValues[i],
                 sizeof(kStockDCHuffmanCodeValues[i])) == 0) {
        found_match = 1;
        stock_table_idx = i;
      }
    }
  } else {
    for (int i = 0; i < kNumStockACHuffmanCodes && !found_match; ++i) {
      if (memcmp(&huff.counts[1], kStockACHuffmanCodeCounts[i],
                 sizeof(kStockACHuffmanCodeCounts[i])) == 0 &&
          memcmp(&huff.values[0], kStockACHuffmanCodeValues[i],
                 sizeof(kStockACHuffmanCodeValues[i])) == 0) {
        found_match = 1;
        stock_table_idx = i;
      }
    }
  }
  WriteBits(1, found_match, storage);
  if (found_match) {
    WriteBits(1, stock_table_idx, storage);
    return true;
  }
  while (max_len > 0 && huff.counts[max_len] == 0) --max_len;
  if (huff.counts[0] != 0 || max_len == 0) {
    return false;
  }
  WriteBits(4, max_len - 1, storage);
  space -= (1 << (kJpegHuffmanMaxBitLength - max_len));
  for (int i = 1; i <= max_len; ++i) {
    int count = huff.counts[i] - (i == max_len ? 1 : 0);
    int count_limit = std::min(max_count - total_count,
                               space >> (kJpegHuffmanMaxBitLength - i));
    if (count > count_limit) {
      BRUNSLI_LOG_DEBUG() << "len = " << i << " count = " << count
                          << " limit = " << count_limit << " space = " << space
                          << " total = " << total_count << BRUNSLI_ENDL();
      return false;
    }
    if (count_limit > 0) {
      int nbits = Log2FloorNonZero(count_limit) + 1;
      WriteBits(nbits, count, storage);
      total_count += count;
      space -= count * (1 << (kJpegHuffmanMaxBitLength - i));
    }
  }
  if (huff.values[total_count] != kJpegHuffmanAlphabetSize) {
    return false;
  }

  PermutationCoder p;
  p.Init(
      is_dc_table
          ? std::vector<uint8_t>(kDefaultDCValues, std::end(kDefaultDCValues))
          : std::vector<uint8_t>(kDefaultACValues, std::end(kDefaultACValues)));
  for (int i = 0; i < total_count; ++i) {
    const int val = huff.values[i];
    int code, nbits;
    if (!p.RemoveValue(val, &code, &nbits)) {
      return false;
    }
    EncodeLimitedVarint(code, 2, (nbits + 1) >> 1, storage);
  }
  return true;
}

bool EncodeScanInfo(const JPEGScanInfo& si, Storage* storage) {
  WriteBits(6, si.Ss, storage);
  WriteBits(6, si.Se, storage);
  WriteBits(4, si.Ah, storage);
  WriteBits(4, si.Al, storage);
  WriteBits(2, si.num_components - 1, storage);
  for (size_t i = 0; i < si.num_components; ++i) {
    const JPEGComponentScanInfo& csi = si.components[i];
    WriteBits(2, csi.comp_idx, storage);
    WriteBits(2, csi.dc_tbl_idx, storage);
    WriteBits(2, csi.ac_tbl_idx, storage);
  }
  int last_block_idx = -1;
  for (const auto& block_idx : si.reset_points) {
    WriteBits(1, 1, storage);
    BRUNSLI_DCHECK(block_idx >= last_block_idx + 1);
    EncodeVarint(block_idx - last_block_idx - 1, 28, storage);
    last_block_idx = block_idx;
  }
  WriteBits(1, 0, storage);

  last_block_idx = 0;
  for (size_t i = 0; i < si.extra_zero_runs.size(); ++i) {
    int block_idx = si.extra_zero_runs[i].block_idx;
    int num = si.extra_zero_runs[i].num_extra_zero_runs;
    BRUNSLI_DCHECK(block_idx >= last_block_idx);
    for (int j = 0; j < num; ++j) {
      WriteBits(1, 1, storage);
      EncodeVarint(block_idx - last_block_idx, 28, storage);
      last_block_idx = block_idx;
    }
  }
  WriteBits(1, 0, storage);

  return true;
}

int MatchComponentIds(const std::vector<JPEGComponent>& comps) {
  if (comps.size() == 1 && comps[0].id == 1) {
    return kComponentIdsGray;
  }
  if (comps.size() == 3) {
    if (comps[0].id == 1 && comps[1].id == 2 && comps[2].id == 3) {
      return kComponentIds123;
    } else if (comps[0].id == 'R' && comps[1].id == 'G' && comps[2].id == 'B') {
      return kComponentIdsRGB;
    }
  }
  return kComponentIdsCustom;
}

void JumpToByteBoundary(Storage* storage) {
  int nbits = storage->pos & 7;
  if (nbits > 0) {
    WriteBits(8 - nbits, 0, storage);
  }
}

bool EncodeAuxData(const JPEGData& jpg, Storage* storage) {
  if (jpg.marker_order.empty() || jpg.marker_order.back() != 0xd9) {
    return false;
  }
  bool have_dri = false;
  size_t num_scans = 0;
  for (size_t i = 0; i < jpg.marker_order.size(); ++i) {
    uint8_t marker = jpg.marker_order[i];
    if (marker < 0xc0) {
      return false;
    }
    WriteBits(6, marker - 0xc0, storage);
    if (marker == 0xdd) have_dri = true;
    if (marker == 0xda) ++num_scans;
  }
  if (have_dri) {
    WriteBits(16, jpg.restart_interval, storage);
  }

  BRUNSLI_DCHECK(jpg.huffman_code.size() < kMaxDHTMarkers);
  for (size_t i = 0; i < jpg.huffman_code.size(); ++i) {
    const bool is_known_last = ((i + 1) == jpg.huffman_code.size());
    WriteBits(1, is_known_last, storage);
    if (!EncodeHuffmanCode(jpg.huffman_code[i], is_known_last, storage)) {
      return false;
    }
  }

  if (num_scans != jpg.scan_info.size()) {
    return false;
  }
  for (size_t i = 0; i < jpg.scan_info.size(); ++i) {
    if (!EncodeScanInfo(jpg.scan_info[i], storage)) {
      return false;
    }
  }
  WriteBits(2, jpg.quant.size() - 1, storage);
  for (size_t i = 0; i < jpg.quant.size(); ++i) {
    WriteBits(2, jpg.quant[i].index, storage);
    if (i != jpg.quant.size() - 1) {
      WriteBits(1, jpg.quant[i].is_last, storage);
    } else if (!jpg.quant[i].is_last) {
      return false;
    }
    WriteBits(4, jpg.quant[i].precision, storage);
  }
  int comp_ids = MatchComponentIds(jpg.components);
  WriteBits(2, comp_ids, storage);
  if (comp_ids == kComponentIdsCustom) {
    for (size_t i = 0; i < jpg.components.size(); ++i) {
      WriteBits(8, jpg.components[i].id, storage);
    }
  }
  size_t nsize = jpg.has_zero_padding_bit ? jpg.padding_bits.size() : 0;
  if (nsize > PaddingBitsLimit(jpg)) return false;
  // we limit to 32b for nsize
  EncodeLimitedVarint(nsize, 8, 4, storage);
  if (nsize > 0) {
    for (size_t i = 0; i < nsize; ++i) {
      WriteBits(1, jpg.padding_bits[i], storage);
    }
  }
  JumpToByteBoundary(storage);
  for (size_t i = 0; i < jpg.inter_marker_data.size(); ++i) {
    const auto& s = jpg.inter_marker_data[i];
    uint8_t buffer[(sizeof(size_t) * 8 + 6) / 7];
    size_t len = EncodeBase128(s.size(), buffer);
    storage->AppendBytes(buffer, len);
    storage->AppendBytes(s.data(), s.size());
  }
  return true;
}

Histogram::Histogram() { Clear(); }

void Histogram::Clear() {
  memset(data_, 0, sizeof(data_));
  total_count_ = 0;
}

void Histogram::AddHistogram(const Histogram& other) {
  for (int i = 0; i < BRUNSLI_ANS_MAX_SYMBOLS; ++i) {
    data_[i] += other.data_[i];
  }
  total_count_ += other.total_count_;
}

void Histogram::Add(size_t val) {
  BRUNSLI_DCHECK(val < BRUNSLI_ANS_MAX_SYMBOLS);
  ++data_[val];
  ++total_count_;
}

void Histogram::Merge(const Histogram& other) {
  if (other.total_count_ == 0) return;
  total_count_ += other.total_count_;
  for (size_t i = 0; i < BRUNSLI_ANS_MAX_SYMBOLS; ++i)
    data_[i] += other.data_[i];
}

void ComputeCoeffOrder(const BlockI32& num_zeros, uint32_t* order) {
  std::vector<std::pair<int, int>> pos_and_val(kDCTBlockSize);
  for (int i = 0; i < kDCTBlockSize; ++i) {
    pos_and_val[i].first = i;
    pos_and_val[i].second = num_zeros[kJPEGNaturalOrder[i]];
  }
  std::stable_sort(
      pos_and_val.begin(), pos_and_val.end(),
      [](const std::pair<int, int>& a, const std::pair<int, int>& b) -> bool {
        return a.second < b.second;
      });
  for (size_t i = 0; i < kDCTBlockSize; ++i) {
    order[i] = kJPEGNaturalOrder[pos_and_val[i].first];
  }
}

void EntropySource::Resize(size_t num_bands) {
  num_bands_ = num_bands;
  histograms_.resize(num_bands * kNumAvrgContexts);
}

void EntropySource::AddCode(size_t code, size_t histo_ix) {
  histograms_[histo_ix].Add(code);
}

std::unique_ptr<EntropyCodes> EntropySource::Finish(
    const std::vector<size_t>& offsets) {
  std::vector<Histogram> histograms;
  histograms.swap(histograms_);
  return std::unique_ptr<EntropyCodes>(
      new EntropyCodes(histograms, num_bands_, offsets));
}

void EntropySource::Merge(const EntropySource& other) {
  BRUNSLI_DCHECK(histograms_.size() >= other.histograms_.size());
  for (size_t i = 0; i < other.histograms_.size(); ++i) {
    histograms_[i].Merge(other.histograms_[i]);
  }
}

EntropyCodes::EntropyCodes(const std::vector<Histogram>& histograms,
                           size_t num_bands,
                           const std::vector<size_t>& offsets) {
  brunsli::ClusterHistograms(histograms, kNumAvrgContexts, num_bands, offsets,
                             kMaxNumberOfHistograms, &clustered_,
                             &context_map_);
}

void EntropyCodes::EncodeContextMap(Storage* storage) const {
  brunsli::EncodeContextMap(context_map_, clustered_.size(), storage);
}

void EntropyCodes::BuildAndStoreEntropyCodes(Storage* storage) {
  ans_tables_.resize(clustered_.size());
  for (size_t i = 0; i < clustered_.size(); ++i) {
    BuildAndStoreANSEncodingData(&clustered_[i].data_[0], &ans_tables_[i],
                                 storage);
  }
}

const ANSTable* EntropyCodes::GetANSTable(int context) const {
  const int entropy_ix = context_map_[context];
  return &ans_tables_[entropy_ix];
}

DataStream::DataStream()
    : pos_(3),
      bw_pos_(0),
      ac_pos0_(1),
      ac_pos1_(2),
      low_(0),
      high_(~0),
      bw_val_(0),
      bw_bitpos_(0) {}

void DataStream::Resize(size_t max_num_code_words) {
  code_words_.resize(max_num_code_words);
}

void DataStream::ResizeForBlock() {
  if (pos_ + kSlackForOneBlock > code_words_.size()) {
    static const double kGrowMult = 1.2;
    const size_t new_size =
        static_cast<size_t>(kGrowMult * code_words_.capacity()) +
        kSlackForOneBlock;
    code_words_.resize(new_size);
  }
}

void DataStream::AddCode(size_t code, size_t band, size_t context,
                         EntropySource* s) {
  size_t histo_ix = band * kNumAvrgContexts + context;
  CodeWord word;
  word.context = static_cast<uint32_t>(histo_ix);
  word.code = static_cast<uint32_t>(code);
  word.nbits = 0;
  word.value = 0;
  BRUNSLI_DCHECK(pos_ < code_words_.size());
  code_words_[pos_++] = word;
  s->AddCode(code, histo_ix);
}

void DataStream::AddBits(int nbits, int bits) {
  bw_val_ |= (bits << bw_bitpos_);
  bw_bitpos_ += nbits;
  if (bw_bitpos_ > 16) {
    CodeWord word;
    word.context = 0;
    word.code = 0;
    word.nbits = 16;
    word.value = bw_val_ & 0xffff;
    code_words_[bw_pos_] = word;
    bw_pos_ = pos_;
    ++pos_;
    bw_val_ >>= 16;
    bw_bitpos_ -= 16;
  }
}

void DataStream::FlushArithmeticCoder() {
  code_words_[ac_pos0_].value = high_ >> 16;
  code_words_[ac_pos1_].value = high_ & 0xffff;
  code_words_[ac_pos0_].nbits = 16;
  code_words_[ac_pos1_].nbits = 16;
  low_ = 0;
  high_ = ~0;
}

void DataStream::FlushBitWriter() {
  code_words_[bw_pos_].nbits = 16;
  code_words_[bw_pos_].value = bw_val_ & 0xffff;
}

// Encodes the next bit to the bit stream, based on the 8-bit precision
// probability, i.e. P(bit = 0) = prob / 256. Statistics are updated in 'p'.
void DataStream::AddBit(Prob* const p, int bit) {
  const uint8_t prob = p->get_proba();
  p->Add(bit);
  const uint32_t diff = high_ - low_;
  const uint32_t split = low_ + (((uint64_t)diff * prob) >> 8);
  if (bit) {
    low_ = split + 1;
  } else {
    high_ = split;
  }
  if (((low_ ^ high_) >> 16) == 0) {
    code_words_[ac_pos0_].value = high_ >> 16;
    code_words_[ac_pos0_].nbits = 16;
    ac_pos0_ = ac_pos1_;
    ac_pos1_ = pos_;
    ++pos_;
    low_ <<= 16;
    high_ <<= 16;
    high_ |= 0xffff;
  }
}

void DataStream::EncodeCodeWords(EntropyCodes* s, Storage* storage) {
  FlushBitWriter();
  FlushArithmeticCoder();
  ANSCoder ans;
  for (int i = pos_ - 1; i >= 0; --i) {
    CodeWord* const word = &code_words_[i];
    if (word->nbits == 0) {
      const ANSEncSymbolInfo info =
          s->GetANSTable(word->context)->info_[word->code];
      word->value = ans.PutSymbol(info, &word->nbits);
    }
  }
  const uint32_t state = ans.GetState();
  uint16_t* out = reinterpret_cast<uint16_t*>(storage->data);
  const uint16_t* out_start = out;
  // Mixed-endian for historical reasons.
  BRUNSLI_UNALIGNED_STORE16LE(out++, state >> 16);
  BRUNSLI_UNALIGNED_STORE16LE(out++, state);
  for (int i = 0; i < pos_; ++i) {
    const CodeWord& word = code_words_[i];
    if (word.nbits) {
      BRUNSLI_UNALIGNED_STORE16LE(out++, word.value);
    }
  }
  storage->pos += (out - out_start) * 16;
}

void EncodeNumNonzeros(size_t val, Prob* p, DataStream* data_stream) {
  BRUNSLI_DCHECK(val < (1u << kNumNonZeroBits));

  // To simplity BST navigation, we use 1-based indexing.
  Prob* bst = p - 1;
  size_t ctx = 1;

  for (size_t mask = 1 << (kNumNonZeroBits - 1); mask != 0; mask >>= 1) {
    const int bit = ((val & mask) != 0);
    data_stream->AddBit(bst + ctx, bit);
    ctx = 2 * ctx + bit;
  }
}

// Or'ing of all coeffs [1..63], for quick zero-test:
coeff_t CollectAllCoeffs(const coeff_t coeffs[kDCTBlockSize]) {
  coeff_t all_coeffs = 0;
  for (int k = 1; all_coeffs == 0 && k < kDCTBlockSize; ++k) {
    all_coeffs |= coeffs[k];
  }
  return all_coeffs;
}

void EncodeCoeffOrder(const uint32_t* order, DataStream* data_stream) {
  uint32_t order_zigzag[kDCTBlockSize];
  for (size_t i = 0; i < kDCTBlockSize; ++i) {
    order_zigzag[i] = kJPEGZigZagOrder[order[i]];
  }
  uint32_t lehmer[kDCTBlockSize];
  ComputeLehmerCode(order_zigzag, kDCTBlockSize, lehmer);
  int tail = kDCTBlockSize - 1;
  while (tail >= 1 && lehmer[tail] == 0) {
    --tail;
  }
  for (int i = 1; i <= tail; ++i) {
    ++lehmer[i];
  }
  static const int kSpan = 16;
  for (int i = 0; i < kDCTBlockSize; i += kSpan) {
    const int start = (i > 0) ? i : 1;
    const int end = i + kSpan;
    int has_non_zero = 0;
    for (int j = start; j < end; ++j) has_non_zero |= lehmer[j];
    if (!has_non_zero) {  // all zero in the span -> escape
      data_stream->AddBits(1, 0);
      continue;
    } else {
      data_stream->AddBits(1, 1);
    }
    for (int j = start; j < end; ++j) {
      int v;
      BRUNSLI_DCHECK(lehmer[j] <= kDCTBlockSize);
      for (v = lehmer[j]; v >= 7; v -= 7) {
        data_stream->AddBits(3, 7);
      }
      data_stream->AddBits(3, v);
    }
  }
}

uint32_t FrameTypeCode(const JPEGData& jpg) {
  uint32_t code = 0;
  int shift = 0;
  for (size_t i = 0; i < jpg.components.size() && i < 4; ++i) {
    uint32_t h_samp = jpg.components[i].h_samp_factor - 1;
    uint32_t v_samp = jpg.components[i].v_samp_factor - 1;
    code |= (h_samp << (shift + 4)) | (v_samp << shift);
    shift += 8;
  }
  return code;
}

bool EncodeSignature(size_t len, uint8_t* data, size_t* pos) {
  if (len < kBrunsliSignatureSize || *pos > len - kBrunsliSignatureSize) {
    return false;
  }
  memcpy(&data[*pos], kBrunsliSignature, kBrunsliSignatureSize);
  *pos += kBrunsliSignatureSize;
  return true;
}

static void EncodeValue(uint8_t tag, size_t value, uint8_t* data, size_t* pos) {
  data[(*pos)++] = ValueMarker(tag);
  *pos += EncodeBase128(value, data + *pos);
}

bool EncodeHeader(const JPEGData& jpg, State* state, uint8_t* data,
                  size_t* len) {
  BRUNSLI_UNUSED(state);
  size_t version = jpg.version;
  bool is_fallback = (version & 1);
  // Fallback can not be combined with anything else.
  if (is_fallback && (version != 1)) return false;
  // Non-fallback image can not be empty.
  if ((!is_fallback && (jpg.width == 0 || jpg.height == 0)) ||
      jpg.components.empty() || jpg.components.size() > kMaxComponents) {
    return false;
  }
  // Only 3 bits are defined.
  if (version & ~7u) return false;

  size_t version_comp = (jpg.components.size() - 1) | (version << 2);
  size_t subsampling = FrameTypeCode(jpg);

  size_t pos = 0;
  EncodeValue(kBrunsliHeaderWidthTag, jpg.width, data, &pos);
  EncodeValue(kBrunsliHeaderHeightTag, jpg.height, data, &pos);
  EncodeValue(kBrunsliHeaderVersionCompTag, version_comp, data, &pos);
  EncodeValue(kBrunsliHeaderSubsamplingTag, subsampling, data, &pos);

  *len = pos;
  return true;
}

bool EncodeMetaData(const JPEGData& jpg, State* state, uint8_t* data,
                    size_t* len) {
  BRUNSLI_UNUSED(state);
  // Concatenate all the (possibly transformed) metadata pieces into one string.
  std::vector<uint8_t> metadata;
  size_t transformed_marker_count = 0;
  for (size_t i = 0; i < jpg.app_data.size(); ++i) {
    const auto& s = jpg.app_data[i];
    Append(&metadata, TransformAppMarker(s, &transformed_marker_count));
  }
  if (transformed_marker_count > kBrunsliShortMarkerLimit) {
    BRUNSLI_LOG_ERROR() << "Too many short markers: "
                        << transformed_marker_count << BRUNSLI_ENDL();
    return false;
  }
  for (const auto& s : jpg.com_data) {
    Append(&metadata, s);
  }
  if (!jpg.tail_data.empty()) {
    const uint8_t marker[] = {0xD9};
    Append(&metadata, marker, 1);
    Append(&metadata, jpg.tail_data);
  }
  if (metadata.empty()) {
    *len = 0;
    return true;
  } else if (metadata.size() == 1) {
    *len = 1;
    data[0] = metadata[0];
    return true;
  }

  // Write base-128 encoding of the original metadata size.
  size_t pos = EncodeBase128(metadata.size(), data);

  // Write the compressed metadata directly to the output.
  size_t compressed_size = *len - pos;
  if (!BrotliEncoderCompress(kBrotliQuality, kBrotliWindowBits,
                             BROTLI_DEFAULT_MODE, metadata.size(),
                             metadata.data(), &compressed_size, &data[pos])) {
    BRUNSLI_LOG_ERROR() << "Brotli compression failed:"
                        << " input size = " << metadata.size()
                        << " pos = " << pos << " len = " << *len
                        << BRUNSLI_ENDL();
    return false;
  }
  pos += compressed_size;
  *len = pos;
  return true;
}

bool EncodeJPEGInternals(const JPEGData& jpg, State* state, uint8_t* data,
                         size_t* len) {
  BRUNSLI_UNUSED(state);
  Storage storage(data, *len);

  if (!EncodeAuxData(jpg, &storage)) {
    return false;
  }

  *len = storage.GetBytesUsed();
  return true;
}

bool EncodeQuantData(const JPEGData& jpg, State* state, uint8_t* data,
                     size_t* len) {
  BRUNSLI_UNUSED(state);
  Storage storage(data, *len);

  if (!EncodeQuantTables(jpg, &storage)) {
    return false;
  }

  *len = storage.GetBytesUsed();
  return true;
}

bool EncodeHistogramData(const JPEGData& jpg, State* state, uint8_t* data,
                         size_t* len) {
  Storage storage(data, *len);

  for (size_t i = 0; i < jpg.components.size(); ++i) {
    WriteBits(3, state->meta[i].context_bits, &storage);
  }

  state->entropy_codes->EncodeContextMap(&storage);

  state->entropy_codes->BuildAndStoreEntropyCodes(&storage);

  *len = storage.GetBytesUsed();
  return true;
}

bool EncodeDCData(const JPEGData& jpg, State* state, uint8_t* data,
                  size_t* len) {
  BRUNSLI_UNUSED(jpg);
  Storage storage(data, *len);

  state->data_stream_dc.EncodeCodeWords(state->entropy_codes, &storage);

  *len = storage.GetBytesUsed();
  return true;
}

bool EncodeACData(const JPEGData& jpg, State* state, uint8_t* data,
                  size_t* len) {
  BRUNSLI_UNUSED(jpg);
  Storage storage(data, *len);

  state->data_stream_ac.EncodeCodeWords(state->entropy_codes, &storage);

  *len = storage.GetBytesUsed();
  return true;
}

typedef bool (*EncodeSectionDataFn)(const JPEGData& jpg, State* state,
                                    uint8_t* data, size_t* len);

bool EncodeSection(const JPEGData& jpg, State* s, uint8_t tag,
                   EncodeSectionDataFn write_section, size_t section_size_bytes,
                   size_t len, uint8_t* data, size_t* pos) {
  // Write the marker byte for the section.
  const size_t pos_start = *pos;
  const uint8_t marker = SectionMarker(tag);
  data[(*pos)++] = marker;

  // Skip enough bytes for a valid (though not necessarily optimal) base-128
  // encoding of the size of the section.
  *pos += section_size_bytes;

  size_t section_size = len - *pos;
  if (!write_section(jpg, s, &data[*pos], &section_size)) {
    return false;
  }
  *pos += section_size;

  if ((section_size >> (7 * section_size_bytes)) > 0) {
    BRUNSLI_LOG_ERROR() << "Section 0x" << std::hex << marker << " size "
                        << std::dec << section_size << " too large for "
                        << section_size_bytes << " bytes base128 number."
                        << BRUNSLI_ENDL();
    return false;
  }

  // Write the final size of the section after the marker byte.
  EncodeBase128Fix(section_size, section_size_bytes, &data[pos_start + 1]);
  return true;
}

namespace internal {
namespace enc {

size_t SampleNumNonZeros(ComponentMeta* m) {
  size_t num_blocks = m->width_in_blocks * m->height_in_blocks;
  if (num_blocks < 32 * 32) return kDCTBlockSize * num_blocks;

  const coeff_t* coeffs = m->ac_coeffs;
  size_t stride = m->ac_stride;
  size_t width_in_blocks = m->width_in_blocks;
  BlockI32& num_zeros = m->num_zeros;

  // For faster compression we only go over a sample of the blocks here.
  static const int kStride = 5;
  size_t total_nonzeros = 0;
  for (size_t i = 0; i < num_blocks; i += kStride) {
    size_t x = i % width_in_blocks;
    size_t y = i / width_in_blocks;
    const coeff_t* block = coeffs + x * kDCTBlockSize + y * stride;
    for (size_t k = 0; k < kDCTBlockSize; ++k) {
      if (block[k] == 0) ++num_zeros[k];
    }
    total_nonzeros += kDCTBlockSize;
  }
  for (size_t i = 0; i < kDCTBlockSize; ++i) total_nonzeros -= num_zeros[i];
  num_zeros[0] = 0;  // DC coefficient is always the first one.
  return total_nonzeros * kStride;
}

int SelectContextBits(size_t num_symbols) {
  static const int kContextBits[33] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2, 3,
      3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6,
  };
  size_t log2_size = Log2FloorNonZero(static_cast<uint32_t>(num_symbols));
  int scheme = kContextBits[log2_size];
  BRUNSLI_DCHECK(scheme < kNumSchemes);
  return scheme;
}

bool PredictDCCoeffs(State* state) {
  std::vector<ComponentMeta>& meta = state->meta;
  for (size_t i = 0; i < meta.size(); ++i) {
    ComponentMeta& m = meta[i];
    const int width = m.width_in_blocks;
    const int height = m.height_in_blocks;
    const int ac_stride = m.ac_stride;
    const int dc_stride = m.dc_stride;
    for (int y = 0; y < height; ++y) {
      const coeff_t* coeffs = m.ac_coeffs + ac_stride * y;
      coeff_t* pred_errors = m.dc_prediction_errors + dc_stride * y;
      for (int x = 0; x < width; ++x) {
        int err =
            coeffs[0] - PredictWithAdaptiveMedian(coeffs, x, y, ac_stride);
        if (std::abs(err) > kBrunsliMaxDCAbsVal) {
          BRUNSLI_LOG_INFO() << "Invalid DC coefficient: " << coeffs[0]
                             << " after prediction: " << err << BRUNSLI_ENDL();
          return false;
        }
        coeffs += kDCTBlockSize;
        *(pred_errors++) = err;
      }
    }
  }
  return true;
}

bool CalculateMeta(const JPEGData& jpg, State* state) {
  const size_t num_components = jpg.components.size();
  std::vector<ComponentMeta>& meta = state->meta;
  meta.resize(num_components);
  for (size_t i = 0; i < num_components; ++i) {
    const JPEGComponent& c = jpg.components[i];
    ComponentMeta& m = meta[i];
    if (c.quant_idx >= jpg.quant.size()) return false;
    const JPEGQuantTable& q = jpg.quant[c.quant_idx];
    m.h_samp = c.h_samp_factor;
    m.v_samp = c.v_samp_factor;
    m.width_in_blocks = jpg.MCU_cols * m.h_samp;
    m.height_in_blocks = jpg.MCU_rows * m.v_samp;
    m.ac_coeffs = &c.coeffs[0];
    m.ac_stride = m.width_in_blocks * kDCTBlockSize;
    m.dc_stride = m.width_in_blocks;
    m.b_stride = m.width_in_blocks;
    memcpy(m.quant.data(), &q.values[0], kDCTBlockSize * sizeof(m.quant[0]));
  }
  return true;
}

void EncodeDC(State* state) {
  const std::vector<ComponentMeta>& meta = state->meta;
  const size_t num_components = meta.size();
  const int mcu_rows = meta[0].height_in_blocks / meta[0].v_samp;
  EntropySource& entropy_source = state->entropy_source;
  DataStream& data_stream = state->data_stream_dc;

  std::vector<ComponentStateDC> comps(num_components);
  size_t total_num_blocks = 0;
  for (int i = 0; i < num_components; ++i) {
    const ComponentMeta& m = meta[i];
    comps[i].SetWidth(m.width_in_blocks);
    total_num_blocks += m.width_in_blocks * m.height_in_blocks;
  }
  entropy_source.Resize(num_components);
  data_stream.Resize(3u * total_num_blocks + 128u);

  // We encode image components in the following interleaved manner:
  //   v_samp[0] rows of 8x8 blocks from component 0
  //   v_samp[1] rows of 8x8 blocks from component 1
  //   v_samp[2] rows of 8x8 blocks from component 2
  //   v_samp[3] rows of 8x8 blocks from component 3 (if present)
  //
  // E.g. in a YUV420 image, we encode 2 rows of 8x8 blocks from Y and then
  // 1 row of 8x8 blocks from U and 1 row of 8x8 blocks from V.
  //
  // In the terminology of the JPEG standard, we encode one row of MCUs at a
  // time, but within this MCU row, we encode the components non-interleaved.
  for (int mcu_y = 0; mcu_y < mcu_rows; ++mcu_y) {
    for (size_t i = 0; i < num_components; ++i) {
      ComponentStateDC* c = &comps[i];
      const ComponentMeta& m = meta[i];
      const int width = c->width;
      const int ac_stride = m.ac_stride;
      const int dc_stride = m.dc_stride;
      const int b_stride = m.b_stride;
      int y = mcu_y * m.v_samp;
      int* prev_sgn = &c->prev_sign[1];
      int* prev_abs = &c->prev_abs_coeff[2];
      for (int iy = 0; iy < m.v_samp; ++iy, ++y) {
        const coeff_t* dc_coeffs_in = m.dc_prediction_errors + y * dc_stride;
        const coeff_t* ac_coeffs_in = m.ac_coeffs + y * ac_stride;
        uint8_t* block_state = m.block_state + y * b_stride;
        for (int x = 0; x < width; ++x) {
          data_stream.ResizeForBlock();
          const coeff_t coeff = dc_coeffs_in[0];
          const int sign = (coeff > 0) ? 1 : (coeff < 0) ? 2 : 0;
          const int absval = (sign == 2) ? -coeff : coeff;
          const coeff_t all_coeffs = coeff | CollectAllCoeffs(ac_coeffs_in);
          const bool is_empty_block = (all_coeffs == 0);
          const int is_empty_ctx =
              IsEmptyBlockContext(&c->prev_is_nonempty[1], x);
          data_stream.AddBit(&c->is_empty_block_prob[is_empty_ctx],
                             !is_empty_block);
          c->prev_is_nonempty[x + 1] = !is_empty_block;
          *block_state = is_empty_block;
          if (!is_empty_block) {
            const int is_zero = (coeff == 0);
            data_stream.AddBit(&c->is_zero_prob, is_zero);
            if (!is_zero) {
              const int avrg_ctx = WeightedAverageContextDC(prev_abs, x);
              const int sign_ctx = prev_sgn[x] * 3 + prev_sgn[x - 1];
              data_stream.AddBit(&c->sign_prob[sign_ctx], sign - 1);
              const size_t zdens_ctx = i;
              if (absval <= kNumDirectCodes) {
                data_stream.AddCode(absval - 1, zdens_ctx,
                                    static_cast<uint32_t>(avrg_ctx),
                                    &entropy_source);
              } else {
                int nbits = Log2FloorNonZero(absval - kNumDirectCodes + 1) - 1;
                data_stream.AddCode(kNumDirectCodes + nbits, zdens_ctx,
                                    avrg_ctx, &entropy_source);
                int extra_bits = absval - (kNumDirectCodes - 1 + (2 << nbits));
                int first_extra_bit = (extra_bits >> nbits) & 1;
                data_stream.AddBit(&c->first_extra_bit_prob[nbits],
                                   first_extra_bit);
                if (nbits > 0) {
                  extra_bits &= (1 << nbits) - 1;
                  data_stream.AddBits(nbits, extra_bits);
                }
              }
            }
          }
          prev_sgn[x] = sign;
          prev_abs[x] = absval;
          ++block_state;
          ++dc_coeffs_in;
          ac_coeffs_in += kDCTBlockSize;
        }
      }
    }
  }
}

void EncodeAC(State* state) {
  const std::vector<ComponentMeta>& meta = state->meta;
  const size_t num_components = meta.size();
  const int mcu_rows = meta[0].height_in_blocks / meta[0].v_samp;
  EntropySource& entropy_source = state->entropy_source;
  DataStream& data_stream = state->data_stream_ac;
  const uint8_t* context_modes =
      kContextAlgorithm + (state->use_legacy_context_model ? 64 : 0);

  size_t num_code_words = 0;
  std::vector<ComponentState> comps(num_components);
  for (size_t i = 0; i < num_components; ++i) {
    const ComponentMeta& m = meta[i];
    const size_t num_blocks = m.width_in_blocks * m.height_in_blocks;
    num_code_words += 2u * m.approx_total_nonzeros + 1024u + 3u * num_blocks;

    // TODO(eustas): what is better - use shared order or "group" order?
    ComputeCoeffOrder(m.num_zeros, &comps[i].order[0]);
    // TODO(eustas): this computation could be shared between "groups".
    ComputeACPredictMultipliers(m.quant.data(), &comps[i].mult_row[0],
                                &comps[i].mult_col[0]);
    comps[i].SetWidth(m.width_in_blocks);
  }

  entropy_source.Resize(state->num_contexts);
  data_stream.Resize(num_code_words);

  for (int i = 0; i < num_components; ++i) {
    EncodeCoeffOrder(&comps[i].order[0], &data_stream);
  }

  // We encode image components in the following interleaved manner:
  //   v_samp[0] rows of 8x8 blocks from component 0
  //   v_samp[1] rows of 8x8 blocks from component 1
  //   v_samp[2] rows of 8x8 blocks from component 2
  //   v_samp[3] rows of 8x8 blocks from component 3 (if present)
  //
  // E.g. in a YUV420 image, we encode 2 rows of 8x8 blocks from Y and then
  // 1 row of 8x8 blocks from U and 1 row of 8x8 blocks from V.
  //
  // In the terminology of the JPEG standard, we encode one row of MCUs at a
  // time, but within this MCU row, we encode the components non-interleaved.
  for (int mcu_y = 0; mcu_y < mcu_rows; ++mcu_y) {
    for (size_t i = 0; i < num_components; ++i) {
      ComponentState* const c = &comps[i];
      const ComponentMeta& m = meta[i];
      const int cur_ctx_bits = m.context_bits;
      const uint32_t* cur_order = c->order;
      const int width = c->width;
      int y = mcu_y * m.v_samp;
      const int ac_stride = m.ac_stride;
      const int b_stride = m.b_stride;
      int prev_row_delta = (1 - 2 * (y & 1)) * (width + 3) * kDCTBlockSize;
      for (int iy = 0; iy < m.v_samp; ++iy, ++y) {
        const coeff_t* coeffs_in = m.ac_coeffs + y * ac_stride;
        const uint8_t* block_state = m.block_state + y * b_stride;
        const coeff_t* prev_row_coeffs = coeffs_in - ac_stride;
        const coeff_t* prev_col_coeffs = coeffs_in - kDCTBlockSize;
        int* prev_sgn = &c->prev_sign[kDCTBlockSize];
        int* prev_abs =
            &c->prev_abs_coeff[((y & 1) * (width + 3) + 2) * kDCTBlockSize];
        for (int x = 0; x < width; ++x) {
          data_stream.ResizeForBlock();
          coeff_t coeffs[kDCTBlockSize] = {0};
          int last_nz = 0;
          const bool is_empty_block = *block_state;
          if (!is_empty_block) {
            for (int k = 1; k < kDCTBlockSize; ++k) {
              const int k_nat = cur_order[k];
              coeffs[k] = coeffs_in[k_nat];
              if (coeffs[k]) last_nz = k;
            }
            const uint8_t nzero_context =
                NumNonzerosContext(c->prev_num_nonzeros.data(), x, y);
            EncodeNumNonzeros(
                last_nz,
                c->num_nonzero_prob + kNumNonZeroTreeSize * nzero_context,
                &data_stream);
          }
          for (int k = kDCTBlockSize - 1; k > last_nz; --k) {
            prev_sgn[k] = 0;
            prev_abs[k] = 0;
          }
          size_t num_nzeros = 0;
          coeff_t encoded_coeffs[kDCTBlockSize] = {0};
          for (int k = last_nz; k >= 1; --k) {
            coeff_t coeff = coeffs[k];
            const int is_zero = (coeff == 0);
            if (k < last_nz) {
              const int bucket = kNonzeroBuckets[num_nzeros - 1];
              const int is_zero_ctx = bucket * kDCTBlockSize + k;
              Prob* const p = &c->is_zero_prob[is_zero_ctx];
              data_stream.AddBit(p, is_zero);
            }
            if (!is_zero) {
              const int sign = (coeff > 0 ? 0 : 1);
              const int absval = sign ? -coeff : coeff;

              const int k_nat = cur_order[k];
              size_t context_type = context_modes[k_nat];
              size_t avg_ctx = 0;
              size_t sign_ctx = kMaxAverageContext;
              if ((context_type & 1) && (y > 0)) {
                if (y > 0) {
                  size_t offset = k_nat & 7;
                  ACPredictContextRow(
                      prev_row_coeffs + offset, encoded_coeffs + offset,
                      &c->mult_col[offset * 8], &avg_ctx, &sign_ctx);
                }
              } else if ((context_type & 2) && (x > 0)) {
                if (x > 0) {
                  size_t offset = k_nat & ~7;
                  ACPredictContextCol(
                      prev_col_coeffs + offset, encoded_coeffs + offset,
                      &c->mult_row[offset], &avg_ctx, &sign_ctx);
                }
              } else if (!context_type) {
                avg_ctx = WeightedAverageContext(prev_abs + k, prev_row_delta);
                sign_ctx = prev_sgn[k] * 3 + prev_sgn[k - kDCTBlockSize];
              }
              sign_ctx = sign_ctx * kDCTBlockSize + k;
              Prob* const sign_p = &c->sign_prob[sign_ctx];
              data_stream.AddBit(sign_p, sign);
              prev_sgn[k] = sign + 1;
              const size_t zdens_ctx =
                  m.context_offset +
                  ZeroDensityContext(num_nzeros, k, cur_ctx_bits);
              if (absval <= kNumDirectCodes) {
                data_stream.AddCode(absval - 1, zdens_ctx, avg_ctx,
                                    &entropy_source);
              } else {
                const int base_code = absval - kNumDirectCodes + 1;
                const int nbits = Log2FloorNonZero(base_code) - 1;
                data_stream.AddCode(kNumDirectCodes + nbits, zdens_ctx,
                                    static_cast<uint32_t>(avg_ctx),
                                    &entropy_source);
                const int extra_bits = base_code - (2 << nbits);
                const int first_extra_bit = (extra_bits >> nbits) & 1;
                Prob* const p = &c->first_extra_bit_prob[k * 10 + nbits];
                data_stream.AddBit(p, first_extra_bit);
                if (nbits > 0) {
                  const int left_over_bits = extra_bits & ((1 << nbits) - 1);
                  data_stream.AddBits(nbits, left_over_bits);
                }
              }
              ++num_nzeros;
              encoded_coeffs[k_nat] = coeff;
              prev_abs[k] = absval;
            } else {
              prev_sgn[k] = 0;
              prev_abs[k] = 0;
            }
          }
          BRUNSLI_DCHECK(num_nzeros <= kNumNonZeroTreeSize);
          c->prev_num_nonzeros[x] = static_cast<uint8_t>(num_nzeros);
          ++block_state;
          coeffs_in += kDCTBlockSize;
          prev_sgn += kDCTBlockSize;
          prev_abs += kDCTBlockSize;
          prev_row_coeffs += kDCTBlockSize;
          prev_col_coeffs += kDCTBlockSize;
        }
        prev_row_delta *= -1;
      }
    }
  }
}

std::unique_ptr<EntropyCodes> PrepareEntropyCodes(State* state) {
  std::vector<ComponentMeta>& meta = state->meta;
  const size_t num_components = meta.size();
  // Prepend DC context group (starts at 0).
  std::vector<size_t> group_context_offsets(1 + num_components);
  for (size_t i = 0; i < num_components; ++i) {
    group_context_offsets[i + 1] = meta[i].context_offset;
  }
  return state->entropy_source.Finish(group_context_offsets);
}

bool BrunsliSerialize(State* state, const JPEGData& jpg, uint32_t skip_sections,
                      uint8_t* data, size_t* len) {
  size_t pos = 0;

  // TODO(eustas): refactor to remove repetitive params.
  bool ok = true;

  const auto encode_section = [&](uint8_t tag, EncodeSectionDataFn fn,
                                  size_t size) {
    return EncodeSection(jpg, state, tag, fn, size, *len, data, &pos);
  };

  if (!(skip_sections & (1u << kBrunsliSignatureTag))) {
    ok = EncodeSignature(*len, data, &pos);
    if (!ok) return false;
  }

  if (!(skip_sections & (1u << kBrunsliHeaderTag))) {
    ok = encode_section(kBrunsliHeaderTag, EncodeHeader, 1);
    if (!ok) return false;
  }

  if (!(skip_sections & (1u << kBrunsliJPEGInternalsTag))) {
    ok = encode_section(kBrunsliJPEGInternalsTag, EncodeJPEGInternals,
                        Base128Size(EstimateAuxDataSize(jpg)));
    if (!ok) return false;
  }

  if (!(skip_sections & (1u << kBrunsliMetaDataTag))) {
    ok = encode_section(kBrunsliMetaDataTag, EncodeMetaData,
                        Base128Size(*len - pos));
    if (!ok) return false;
  }

  if (!(skip_sections & (1u << kBrunsliQuantDataTag))) {
    ok = encode_section(kBrunsliQuantDataTag, EncodeQuantData, 2);
    if (!ok) return false;
  }

  if (!(skip_sections & (1u << kBrunsliHistogramDataTag))) {
    ok = encode_section(kBrunsliHistogramDataTag, EncodeHistogramData,
                        Base128Size(*len - pos));
    if (!ok) return false;
  }

  if (!(skip_sections & (1u << kBrunsliDCDataTag))) {
    ok = encode_section(kBrunsliDCDataTag, EncodeDCData,
                        Base128Size(*len - pos));
    if (!ok) return false;
  }

  if (!(skip_sections & (1u << kBrunsliACDataTag))) {
    ok = encode_section(kBrunsliACDataTag, EncodeACData,
                        Base128Size(*len - pos));
    if (!ok) return false;
  }

  *len = pos;
  return true;
}

}  // namespace enc
}  // namespace internal

/* Regular Brunsli workflow.
 *
 * For "groups" workflow, few more stages are required, see comments.
 */
bool BrunsliEncodeJpeg(const JPEGData& jpg, uint8_t* data, size_t* len) {
  State state;
  std::vector<ComponentMeta>& meta = state.meta;
  size_t num_components = jpg.components.size();
  state.use_legacy_context_model = !(jpg.version & 2);

  if (!CalculateMeta(jpg, &state)) return false;
  // Groups workflow: update width_in_blocks, height_in_blocks, ac_coeffs.

  for (size_t i = 0; i < num_components; ++i) {
    meta[i].approx_total_nonzeros = SampleNumNonZeros(&state.meta[i]);
  }
  // Groups workflow: reduce approx_total_nonzeros.
  for (size_t i = 0; i < num_components; ++i) {
    meta[i].context_bits = SelectContextBits(meta[i].approx_total_nonzeros + 1);
  }
  // Groups workflow: distribute context_bits.

  // First `num_components` contexts are used for DC.
  size_t num_contexts = num_components;
  for (size_t i = 0; i < num_components; ++i) {
    meta[i].context_offset = num_contexts;
    num_contexts += kNumNonzeroContextSkip[meta[i].context_bits];
  }
  state.num_contexts = num_contexts;

  std::vector<std::vector<coeff_t>> dc_prediction_errors(num_components);
  for (size_t i = 0; i < num_components; ++i) {
    dc_prediction_errors[i].resize(meta[i].width_in_blocks *
                                   meta[i].height_in_blocks);
    meta[i].dc_prediction_errors = dc_prediction_errors[i].data();
  }

  if (!PredictDCCoeffs(&state)) return false;

  std::vector<std::vector<uint8_t>> block_state(num_components);
  for (size_t i = 0; i < num_components; ++i) {
    block_state[i].resize(meta[i].width_in_blocks * meta[i].height_in_blocks);
    meta[i].block_state = block_state[i].data();
  }

  EncodeDC(&state);

  EncodeAC(&state);

  // Groups workflow: merge histograms.
  std::unique_ptr<EntropyCodes> entropy_codes = PrepareEntropyCodes(&state);
  state.entropy_codes = entropy_codes.get();
  // Groups workflow: distribute codes.

  // Groups workflow: apply corresponding skip masks.
  return BrunsliSerialize(&state, jpg, 0, data, len);
}

#if defined(BRUNSLI_EXTRA_API)
// The memory usage of BrunsliEncodeJpeg() looks roughly like this:
//   +-----------------+------------------------+
//   | BrotliCompress  | State::entropy_source  |
//   | (brotli_peak)   | (entropy_source_size)  |
//   +-----------------+------------------------+
//                     | State::data_stream_dc  |
//                     | State::data_stream_ac  |
//                     | (data_stream_size)     |
//                     +------------------------+
//                     | vector<ComponentState> |
//                     | (component_state_size) |
//                     +------------------------+
size_t EstimateBrunsliEncodePeakMemoryUsage(size_t jpg_size,
                                            const JPEGData& jpg) {
  std::vector<uint8_t> tmp;
  size_t metadata_size = 0;
  for (const auto& s : jpg.app_data) {
    metadata_size += TransformApp0Marker(s, &tmp) ? tmp.size() : s.size();
  }
  for (const auto& s : jpg.com_data) {
    metadata_size += s.size();
  }
  size_t brotli_peak = 0;
  if (metadata_size > 1) {
    brotli_peak = BrotliEncoderEstimatePeakMemoryUsage(
        kBrotliQuality, kBrotliWindowBits, metadata_size);
  }
  size_t ncomp = jpg.components.size();
  size_t total_num_blocks = 0;
  size_t component_state_size = 0;
  for (size_t i = 0; i < ncomp; ++i) {
    const JPEGComponent& c = jpg.components[i];
    total_num_blocks += c.num_blocks;
    component_state_size += ComponentState::SizeInBytes(c.width_in_blocks);
  }
  // Since we do not have access to the coefficient data at this point, we
  // estimate the number of nonzero coefficients by assuming that each of them
  // takes about 5 bits compressed in the jpeg format.
  size_t nonzeros_per_block = std::max(
      size_t{1}, std::min(size_t{64}, (jpg_size * 8) / (total_num_blocks * 5)));
  size_t nonzeros = nonzeros_per_block * total_num_blocks;
  int context_bits = SelectContextBits(nonzeros);
  size_t ncontexts = ncomp * (kNumNonzeroContextSkip[context_bits] + 1);
  // We have ncontexts * kNumAvrgContext histograms for both the raw and
  // clustered histogram set.
  size_t entropy_source_size =
      2 * ncontexts * kNumAvrgContexts * sizeof(Histogram);
  size_t ncodewords = std::max(
      size_t{1} << 18, 2 * nonzeros + 6 * total_num_blocks + ncomp * 1024);
  size_t data_stream_size = ncodewords * 8;  // 8 = sizeof(DataStream::CodeWord)
  size_t brunsli_peak =
      entropy_source_size + data_stream_size + component_state_size;
  return std::max(brotli_peak, brunsli_peak);
}
#endif  // defined(BRUNSLI_EXTRA_API)

// bypass mode
const size_t kMaxBypassHeaderSize = 5 * 6;  // = 5x tag + EncodeBase128() call.
size_t GetBrunsliBypassSize(size_t jpg_size) {
  return jpg_size + kBrunsliSignatureSize + kMaxBypassHeaderSize;
}

bool EncodeOriginalJpg(const JPEGData& jpg, State* state, uint8_t* data,
                       size_t* len) {
  BRUNSLI_UNUSED(state);
  if (jpg.original_jpg == NULL || jpg.original_jpg_size > *len) {
    return false;
  }
  memcpy(data, jpg.original_jpg, jpg.original_jpg_size);
  *len = jpg.original_jpg_size;
  return true;
}

bool BrunsliEncodeJpegBypass(const uint8_t* jpg_data, size_t jpg_data_len,
                             uint8_t* data, size_t* len) {
  size_t pos = 0;
  if (!EncodeSignature(*len, data, &pos)) {
    return false;
  }
  JPEGData jpg;
  if (!ReadJpeg(jpg_data, jpg_data_len, JPEG_READ_HEADER, &jpg)) {
    // mark JPEG data as invalid
    jpg.width = 0;
    jpg.height = 0;
    jpg.components.resize(1);
    jpg.components[0].h_samp_factor = 1;
    jpg.components[0].v_samp_factor = 1;
  }
  jpg.version = 1;
  jpg.original_jpg = jpg_data;
  jpg.original_jpg_size = jpg_data_len;

  // TODO(eustas): could use simpler serialization here
  State state;

  if (!EncodeSection(jpg, &state, kBrunsliHeaderTag, EncodeHeader,
                     1, *len, data, &pos)) {
    return false;
  }
  if (!EncodeSection(jpg, &state, kBrunsliOriginalJpgTag, EncodeOriginalJpg,
                     Base128Size(jpg_data_len), *len, data, &pos)) {
    return false;
  }
  *len = pos;
  return true;
}

}  // namespace brunsli
