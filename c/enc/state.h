// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

/* Internal encoder state and workflow methods declaration. */

#ifndef BRUNSLI_ENC_STATE_H_
#define BRUNSLI_ENC_STATE_H_

#include <array>
#include <memory>
#include <vector>

#include "../common/distributions.h"
#include <brunsli/jpeg_data.h>
#include "../common/platform.h"
#include <brunsli/types.h>
#include "./ans_encode.h"
#include "./write_bits.h"

namespace brunsli {
namespace internal {
namespace enc {

typedef std::array<int32_t, kDCTBlockSize> BlockI32;

struct ComponentMeta {
  size_t context_offset;
  size_t approx_total_nonzeros;
  int32_t h_samp;
  int32_t v_samp;
  int32_t context_bits;
  int32_t ac_stride;
  int32_t dc_stride;
  int32_t b_stride;
  int32_t width_in_blocks;
  int32_t height_in_blocks;
  const coeff_t* ac_coeffs;
  coeff_t* dc_prediction_errors;
  // TODO(eustas): investigate bit fields.
  uint8_t* block_state;
  BlockI32 num_zeros;
  BlockI32 quant;
};

struct Histogram {
  Histogram();
  void Clear();
  void AddHistogram(const Histogram& other);
  void Add(size_t val);
  void Merge(const Histogram& other);

  int data_[BRUNSLI_ANS_MAX_SYMBOLS];
  int total_count_;
  double bit_cost_;
};

class EntropyCodes {
 public:
  EntropyCodes(const std::vector<Histogram>& histograms, size_t num_bands,
               const std::vector<size_t>& offsets);
  // GCC declares it won't apply RVO, even if it actually does.
  // EntropyCodes(const EntropyCodes&) = delete;
  void EncodeContextMap(Storage* storage) const;
  void BuildAndStoreEntropyCodes(Storage* storage);
  const ANSTable* GetANSTable(int context) const;

 private:
  static const size_t kMaxNumberOfHistograms = 256;

  std::vector<Histogram> clustered_;
  std::vector<uint32_t> context_map_;
  std::vector<ANSTable> ans_tables_;
};

// Manages building of the histograms of an entropy source.
class EntropySource {
 public:
  EntropySource() : num_bands_(0) {}
  void Resize(size_t num_bands);
  void AddCode(size_t code, size_t histo_ix);
  void Merge(const EntropySource& other);
  std::unique_ptr<EntropyCodes> Finish(const std::vector<size_t>& offsets);

 private:
  size_t num_bands_;
  std::vector<Histogram> histograms_;
};

// Manages the multiplexing of the ANS-coded and arithmetic coded bits.
class DataStream {
 public:
  DataStream();
  void Resize(size_t max_num_code_words);
  void ResizeForBlock();
  void AddCode(size_t code, size_t band, size_t context, EntropySource* s);
  void AddBits(int nbits, int bits);
  void FlushArithmeticCoder();
  void FlushBitWriter();
  // Encodes the next bit to the bit stream, based on the 8-bit precision
  // probability, i.e. P(bit = 0) = prob / 256. Statistics are updated in 'p'.
  void AddBit(Prob* const p, int bit);
  void EncodeCodeWords(EntropyCodes* s, Storage* storage);

 private:
  struct CodeWord {
    // Add a constructor that does nothing (unlike the default one) to avoid
    // initializing to values that are unused anyway.
    CodeWord() {}
    uint32_t context;
    uint16_t value;
    uint8_t code;
    uint8_t nbits;
  };

  static const size_t kSlackForOneBlock = 1024;

  int pos_;
  int bw_pos_;
  int ac_pos0_;
  int ac_pos1_;
  uint32_t low_;
  uint32_t high_;
  uint32_t bw_val_;
  int bw_bitpos_;
  std::vector<CodeWord> code_words_;
};

struct State {
  EntropySource entropy_source;
  EntropyCodes* entropy_codes;
  DataStream data_stream_dc;
  DataStream data_stream_ac;

  std::vector<ComponentMeta> meta;
  size_t num_contexts;
  bool use_legacy_context_model = false;
};

// Encoder workflow:
bool CalculateMeta(const JPEGData& jpg, State* state);
size_t SampleNumNonZeros(ComponentMeta* m);
int SelectContextBits(size_t num_symbols);
bool PredictDCCoeffs(State* state);
void EncodeDC(State* state);
void EncodeAC(State* state);
std::unique_ptr<EntropyCodes> PrepareEntropyCodes(State* state);
bool BrunsliSerialize(State* state, const JPEGData& jpg, uint32_t skip_sections,
                      uint8_t* data, size_t* len);

}  // namespace enc
}  // namespace internal
}  // namespace brunsli

#endif  // BRUNSLI_ENC_STATE_H_
