// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Functions for writing encoding / decoding Brunsli in "groups" mode.

#include "./groups.h"

#include <algorithm>
#include <atomic>
#include <vector>

#include "../common/constants.h"
#include "../common/context.h"
#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>
#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>
#include "../dec/state.h"
#include <brunsli/brunsli_encode.h>
#include <brunsli/jpeg_data_reader.h>
#include "../enc/state.h"

namespace brunsli {

namespace {

bool SkipSection(const uint8_t** data, size_t len) {
  size_t section_len = 0;
  uint64_t b = 0x80;
  size_t off = 1;
  for (size_t i = 0; (i < 9) && (b & 0x80u); ++i) {
    if (off >= len) return false;
    b = (*data)[off++];
    section_len |= (b & 0x7Fu) << (i * 7);
  }
  if ((b & 0x80u) != 0) return false;
  off += section_len;
  if (off > len) return false;
  *data += off;
  return true;
}

}  // namespace

void SequentialExecutor(const Runnable& runnable, size_t num_tasks) {
  for (size_t i = 0; i < num_tasks; ++i) runnable(i);
}

ParallelExecutor::ParallelExecutor(size_t num_threads)
    : num_threads(num_threads) {
  const auto worker = [this]() {
    while (true) {
      {
        std::unique_lock<std::mutex> lock(this->lock);
        start_latch.wait(lock, [this] { return next_task.load() < num_tasks; });
        busy_count++;
        if (terminate) {
          finish_latch.notify_one();
          return;
        }
      }
      while (true) {
        size_t my_task = next_task++;
        if (my_task >= num_tasks) break;
        (*runnable)(my_task);
      }
      {
        std::lock_guard<std::mutex> lock(this->lock);
        busy_count--;
        finish_latch.notify_one();
      }
    }
  };
  futures.reserve(num_threads);
  for (size_t i = 0; i < num_threads; ++i) {
    futures.push_back(std::async(std::launch::async, worker));
  }
}

ParallelExecutor::~ParallelExecutor() {
  std::unique_lock<std::mutex> lock(this->lock);
  terminate = true;
  next_task.store(1);
  this->num_tasks = 1;
  this->runnable = nullptr;
  start_latch.notify_all();
  finish_latch.wait(lock, [this] { return busy_count.load() == num_threads; });
}

Executor ParallelExecutor::getExecutor() {
  return [this](const Runnable& runnable, size_t num_tasks) {
    return execute(runnable, num_tasks);
  };
}

void ParallelExecutor::execute(const Runnable& runnable, size_t num_tasks) {
  std::unique_lock<std::mutex> lock(this->lock);
  next_task.store(0);
  this->num_tasks = num_tasks;
  this->runnable = &runnable;
  start_latch.notify_all();
  finish_latch.wait(lock, [this, num_tasks] {
    return (next_task.load() >= num_tasks) && (busy_count.load() == 0);
  });
}

bool EncodeGroups(const brunsli::JPEGData& jpg, uint8_t* data, size_t* len,
                  size_t ac_group_dim, size_t dc_group_dim,
                  Executor* executor) {
  using ::brunsli::internal::enc::BlockI32;
  using ::brunsli::internal::enc::ComponentMeta;
  using ::brunsli::internal::enc::DataStream;
  using ::brunsli::internal::enc::EntropyCodes;
  using ::brunsli::internal::enc::EntropySource;
  using ::brunsli::internal::enc::Histogram;
  using ::brunsli::internal::enc::SelectContextBits;
  using ::brunsli::internal::enc::State;

  if ((ac_group_dim & (ac_group_dim - 1)) != 0) return false;
  if ((dc_group_dim & (dc_group_dim - 1)) != 0) return false;
  if ((dc_group_dim % ac_group_dim) != 0) return false;
  if ((ac_group_dim % jpg.max_h_samp_factor) != 0) return false;
  if ((ac_group_dim % jpg.max_v_samp_factor) != 0) return false;

  size_t num_components = jpg.components.size();

  std::vector<size_t> approx_total_nonzeros(num_components);

  size_t width_in_blocks = jpg.MCU_cols * jpg.max_h_samp_factor;
  size_t height_in_blocks = jpg.MCU_rows * jpg.max_v_samp_factor;

  size_t w_ac = (width_in_blocks + ac_group_dim - 1) / ac_group_dim;
  size_t h_ac = (height_in_blocks + ac_group_dim - 1) / ac_group_dim;

  size_t w_dc = (width_in_blocks + dc_group_dim - 1) / dc_group_dim;
  size_t h_dc = (height_in_blocks + dc_group_dim - 1) / dc_group_dim;

  std::vector<std::vector<brunsli::coeff_t>> dc_prediction_errors(
      num_components);
  std::vector<std::vector<uint8_t>> block_state(num_components);
  for (size_t i = 0; i < num_components; ++i) {
    const JPEGComponent& c = jpg.components[i];
    dc_prediction_errors[i].resize(c.width_in_blocks * c.height_in_blocks);
    block_state[i].resize(c.width_in_blocks * c.height_in_blocks);
  }

  State state;
  std::vector<State> dc_state(w_dc * h_dc);
  std::vector<State> ac_state(w_ac * h_ac);

  if (!CalculateMeta(jpg, &state)) return false;
  for (size_t c = 0; c < num_components; ++c) {
    ComponentMeta& m = state.meta[c];
    m.dc_prediction_errors = dc_prediction_errors[c].data();
    m.block_state = block_state[c].data();
  }

  for (size_t y = 0; y < h_dc; ++y) {
    for (size_t x = 0; x < w_dc; ++x) {
      State& s = dc_state[x + y * w_dc];
      std::vector<ComponentMeta>& meta = s.meta;
      if (!CalculateMeta(jpg, &s)) return false;
      for (size_t c = 0; c < num_components; ++c) {
        ComponentMeta& m = meta[c];
        size_t h_group_dim = m.h_samp * dc_group_dim / jpg.max_h_samp_factor;
        size_t first_x = x * h_group_dim;
        size_t last_x =
            std::min<size_t>(first_x + h_group_dim, m.width_in_blocks);
        size_t v_group_dim = m.v_samp * dc_group_dim / jpg.max_v_samp_factor;
        size_t first_y = y * v_group_dim;
        size_t last_y =
            std::min<size_t>(first_y + v_group_dim, m.height_in_blocks);
        m.ac_coeffs += first_x * brunsli::kDCTBlockSize + first_y * m.ac_stride;
        m.width_in_blocks = last_x - first_x;
        m.height_in_blocks = last_y - first_y;
        m.dc_prediction_errors =
            dc_prediction_errors[c].data() + first_x + first_y * m.dc_stride;
        m.block_state = block_state[c].data() + first_x + first_y * m.b_stride;
      }
    }
  }

  for (size_t y = 0; y < h_ac; ++y) {
    for (size_t x = 0; x < w_ac; ++x) {
      State& s = ac_state[x + y * w_ac];
      std::vector<ComponentMeta>& meta = s.meta;
      if (!CalculateMeta(jpg, &s)) return false;
      for (size_t c = 0; c < num_components; ++c) {
        ComponentMeta& m = meta[c];
        size_t h_group_dim = m.h_samp * ac_group_dim / jpg.max_h_samp_factor;
        size_t first_x = x * h_group_dim;
        size_t last_x =
            std::min<size_t>(first_x + h_group_dim, m.width_in_blocks);
        size_t v_group_dim = m.v_samp * ac_group_dim / jpg.max_v_samp_factor;
        size_t first_y = y * v_group_dim;
        size_t last_y =
            std::min<size_t>(first_y + v_group_dim, m.height_in_blocks);
        m.ac_coeffs += first_x * brunsli::kDCTBlockSize + first_y * m.ac_stride;
        m.width_in_blocks = last_x - first_x;
        m.height_in_blocks = last_y - first_y;
        m.dc_prediction_errors =
            dc_prediction_errors[c].data() + first_x + first_y * m.dc_stride;
        m.block_state = block_state[c].data() + first_x + first_y * m.b_stride;
      }
    }
  }

  const auto sample_nonzeros = [num_components, &ac_state](size_t idx) {
    for (size_t c = 0; c < num_components; ++c) {
      ComponentMeta& m = ac_state[idx].meta[c];
      m.approx_total_nonzeros = SampleNumNonZeros(&m);
    }
  };
  (*executor)(sample_nonzeros, ac_state.size());

  // Groups workflow: reduce approx_total_nonzeros.
  for (size_t y = 0; y < h_ac; ++y) {
    for (size_t x = 0; x < w_ac; ++x) {
      for (size_t c = 0; c < num_components; ++c) {
        approx_total_nonzeros[c] +=
            ac_state[x + y * w_ac].meta[c].approx_total_nonzeros;
      }
    }
  }

  int32_t num_contexts = num_components;
  for (size_t c = 0; c < num_components; ++c) {
    ComponentMeta& m = state.meta[c];
    m.context_bits = SelectContextBits(approx_total_nonzeros[c] + 1);
    m.context_offset = num_contexts;
    num_contexts += brunsli::kNumNonzeroContextSkip[m.context_bits];
  }
  state.num_contexts = num_contexts;

  // Groups workflow: distribute context_bits.
  for (size_t y = 0; y < h_ac; ++y) {
    for (size_t x = 0; x < w_ac; ++x) {
      State& s = ac_state[x + y * w_ac];
      for (size_t c = 0; c < num_components; ++c) {
        ComponentMeta& m = state.meta[c];
        s.meta[c].context_bits = m.context_bits;
        s.meta[c].context_offset = m.context_offset;
      }
      s.num_contexts = state.num_contexts;
    }
  }

  for (size_t y = 0; y < h_dc; ++y) {
    for (size_t x = 0; x < w_dc; ++x) {
      State& s = dc_state[x + y * w_dc];
      for (size_t c = 0; c < num_components; ++c) {
        ComponentMeta& m = state.meta[c];
        s.meta[c].context_bits = m.context_bits;
        s.meta[c].context_offset = m.context_offset;
      }
      s.num_contexts = state.num_contexts;
    }
  }

  std::atomic<bool> failed{false};
  const auto encode_dc = [&failed, &dc_state](size_t idx) {
    if (failed.load()) return;
    if (!PredictDCCoeffs(&dc_state[idx])) failed.store(true);
    if (failed.load()) return;
    EncodeDC(&dc_state[idx]);
  };
  (*executor)(encode_dc, dc_state.size());
  if (failed.load()) return false;

  const auto encode_ac = [&ac_state](size_t idx) {
    EncodeAC(&ac_state[idx]);
  };
  (*executor)(encode_ac, ac_state.size());

  // Groups workflow: merge histograms.
  // TODO(eustas): SIMDify.
  state.entropy_source.Resize(num_contexts);
  for (size_t y = 0; y < h_dc; ++y) {
    for (size_t x = 0; x < w_dc; ++x) {
      state.entropy_source.Merge(dc_state[x + y * w_dc].entropy_source);
    }
  }
  for (size_t y = 0; y < h_ac; ++y) {
    for (size_t x = 0; x < w_ac; ++x) {
      state.entropy_source.Merge(ac_state[x + y * w_ac].entropy_source);
    }
  }

  std::unique_ptr<EntropyCodes> entropy_codes = PrepareEntropyCodes(&state);

  std::vector<std::vector<uint8_t>> output;
  output.resize(1 + dc_state.size() + ac_state.size());

  // TODO(eustas): pull entropy codes serialization "side effect".
  {
    std::vector<uint8_t>& part = output[0];
    state.entropy_codes = entropy_codes.get();
    size_t part_size = 20480;
    for (size_t i = 0; i < jpg.inter_marker_data.size(); ++i) {
      part_size += 5 + jpg.inter_marker_data[i].size();
    }
    for (const std::vector<uint8_t>& chunk : jpg.app_data) {
      part_size += chunk.size();
    }
    for (const std::vector<uint8_t>& chunk : jpg.com_data) {
      part_size += chunk.size();
    }
    part_size += jpg.tail_data.size();
    // TODO(eustas): take into account histograms.
    part.resize(part_size);
    uint32_t skip_flags =
        (1u << brunsli::kBrunsliDCDataTag) | (1u << brunsli::kBrunsliACDataTag);
    if (!BrunsliSerialize(&state, jpg, skip_flags, part.data(), &part_size)) {
      return false;
    }
    part.resize(part_size);
  }

  const auto serialize = [&](size_t idx) {
    if (failed.load()) return;
    std::vector<uint8_t>& part = output[idx];
    if (idx == 0) return;
    idx--;
    if (idx < dc_state.size()) {
      State& s = dc_state[idx];
      // TODO(eustas): reduce for subsampled
      size_t part_size = 128 * (128 + 16) * jpg.components.size();
      part.resize(part_size);
      s.entropy_codes = entropy_codes.get();
      uint32_t skip_flags = ~(1u << brunsli::kBrunsliDCDataTag);
      bool ok = BrunsliSerialize(&s, jpg, skip_flags, part.data(), &part_size);
      if (ok) {
        part.resize(part_size);
      } else {
        failed.store(true);
      }
      return;
    }
    idx -= dc_state.size();
    if (idx < ac_state.size()) {
      State& s = ac_state[idx];
      // TODO(eustas): reduce for subsampled
      size_t part_size = 32 * 32 * 63 * jpg.components.size();
      part.resize(part_size);
      s.entropy_codes = entropy_codes.get();
      uint32_t skip_flags = ~(1u << brunsli::kBrunsliACDataTag);
      bool ok = BrunsliSerialize(&s, jpg, skip_flags, part.data(), &part_size);
      if (ok) {
        part.resize(part_size);
      } else{
        failed.store(true);
      }
      return;
    }
    failed.store(true);
  };
  (*executor)(serialize, output.size());
  if (failed.load()) return false;

  size_t capacity = *len;
  size_t size = 0;
  for (const std::vector<uint8_t>& part : output) {
    if (size + part.size() > capacity) return false;
    memcpy(data, part.data(), part.size());
    size += part.size();
    data += part.size();
  }
  *len = size;

  return true;
}

bool DecodeGroups(const uint8_t* data, size_t len, brunsli::JPEGData* jpg,
                  size_t ac_group_dim, size_t dc_group_dim,
                  Executor* executor) {
  using ::brunsli::BrunsliStatus;
  using ::brunsli::internal::dec::BlockI32;
  using ::brunsli::internal::dec::ComponentMeta;
  using ::brunsli::internal::dec::PrepareMeta;
  using ::brunsli::internal::dec::ProcessJpeg;
  using ::brunsli::internal::dec::Stage;
  using ::brunsli::internal::dec::State;
  using ::brunsli::internal::dec::WarmupMeta;

  if ((ac_group_dim & (ac_group_dim - 1)) != 0) return false;
  if ((dc_group_dim & (dc_group_dim - 1)) != 0) return false;
  if ((dc_group_dim % ac_group_dim) != 0) return false;

  const uint8_t* data_end = data + len;
  const uint8_t* chunk_end = data;
  const uint8_t* chunk_start = chunk_end;
  // Signature / Header / Meta / Internals / Quant / Histo.
  for (size_t i = 0; i < 6; ++i) {
    if (!SkipSection(&chunk_end, data_end - chunk_end)) return false;
  }

  // Common sections.
  State state;
  state.data = chunk_start;
  state.len = chunk_end - chunk_start;
  chunk_start = chunk_end;

  BrunsliStatus status = ProcessJpeg(&state, jpg);
  if (status != BrunsliStatus::BRUNSLI_NOT_ENOUGH_DATA) return false;
  WarmupMeta(jpg, &state);

  if ((ac_group_dim % jpg->max_h_samp_factor) != 0) return false;
  if ((ac_group_dim % jpg->max_v_samp_factor) != 0) return false;

  size_t num_components = jpg->components.size();

  size_t width_in_blocks = jpg->MCU_cols * jpg->max_h_samp_factor;
  size_t height_in_blocks = jpg->MCU_rows * jpg->max_v_samp_factor;

  size_t w_ac = (width_in_blocks + ac_group_dim - 1) / ac_group_dim;
  size_t h_ac = (height_in_blocks + ac_group_dim - 1) / ac_group_dim;

  size_t w_dc = (width_in_blocks + dc_group_dim - 1) / dc_group_dim;
  size_t h_dc = (height_in_blocks + dc_group_dim - 1) / dc_group_dim;

  std::vector<const uint8_t*> dc_section_start(h_dc * w_dc);
  std::vector<size_t> dc_section_length(h_dc * w_dc);
  for (size_t y = 0; y < h_dc; ++y) {
    for (size_t x = 0; x < w_dc; ++x) {
      if (!SkipSection(&chunk_end, data_end - chunk_end)) return false;
      size_t idx = x + w_dc * y;
      dc_section_start[idx] = chunk_start;
      dc_section_length[idx] = chunk_end - chunk_start;
      chunk_start = chunk_end;
    }
  }

  std::vector<const uint8_t*> ac_section_start(h_ac * w_ac);
  std::vector<size_t> ac_section_length(h_ac * w_ac);
  for (size_t y = 0; y < h_ac; ++y) {
    for (size_t x = 0; x < w_ac; ++x) {
      if (!SkipSection(&chunk_end, data_end - chunk_end)) return false;
      size_t idx = x + w_ac * y;
      ac_section_start[idx] = chunk_start;
      ac_section_length[idx] = chunk_end - chunk_start;
      chunk_start = chunk_end;
    }
  }
  if (chunk_end != data_end) return false;

  std::atomic<bool> failed{false};
  const auto decode_dc = [&](size_t idx) {
    if (failed.load()) return;
    size_t y = idx / w_dc;
    size_t x = idx % w_dc;
    State dc_state;
    dc_state.stage = Stage::SECTION;
    dc_state.tags_met = ~(1u << brunsli::kBrunsliDCDataTag);
    dc_state.data = dc_section_start[idx];
    dc_state.len = dc_section_length[idx];

    dc_state.context_map = state.context_map;
    dc_state.entropy_codes = state.entropy_codes;

    std::vector<ComponentMeta>& meta = dc_state.meta;

    PrepareMeta(jpg, &dc_state);
    dc_state.is_storage_allocated = true;
    WarmupMeta(jpg, &dc_state);
    for (size_t c = 0; c < num_components; ++c) {
      ComponentMeta& m = meta[c];
      size_t h_group_dim = m.h_samp * dc_group_dim / jpg->max_h_samp_factor;
      size_t first_x = x * h_group_dim;
      size_t last_x =
          std::min<size_t>(first_x + h_group_dim, m.width_in_blocks);
      size_t v_group_dim = m.v_samp * dc_group_dim / jpg->max_v_samp_factor;
      size_t first_y = y * v_group_dim;
      size_t last_y =
          std::min<size_t>(first_y + v_group_dim, m.height_in_blocks);
      m.ac_coeffs += first_x * brunsli::kDCTBlockSize + first_y * m.ac_stride;
      m.block_state =
          state.meta[c].block_state + first_x + first_y * m.b_stride;
      m.width_in_blocks = last_x - first_x;
      m.height_in_blocks = last_y - first_y;
    }

    status = ProcessJpeg(&dc_state, jpg);
    if (status != BrunsliStatus::BRUNSLI_OK) failed.store(true);
  };
  (*executor)(decode_dc, dc_section_start.size());
  if (failed.load()) return false;

  const auto decode_ac = [&](size_t idx) {
    if (failed.load()) return;
    size_t y = idx / w_ac;
    size_t x = idx % w_ac;
    State ac_state;
    ac_state.stage = Stage::SECTION;
    ac_state.tags_met = ~(1u << brunsli::kBrunsliACDataTag);
    ac_state.data = ac_section_start[idx];
    ac_state.len = ac_section_length[idx];

    ac_state.context_map = state.context_map;
    ac_state.entropy_codes = state.entropy_codes;

    std::vector<ComponentMeta>& meta = ac_state.meta;

    PrepareMeta(jpg, &ac_state);
    ac_state.is_storage_allocated = true;
    WarmupMeta(jpg, &ac_state);
    for (size_t c = 0; c < num_components; ++c) {
      ComponentMeta& m = meta[c];
      size_t h_group_dim = m.h_samp * ac_group_dim / jpg->max_h_samp_factor;
      size_t first_x = x * h_group_dim;
      size_t last_x =
          std::min<size_t>(first_x + h_group_dim, m.width_in_blocks);
      size_t v_group_dim = m.v_samp * ac_group_dim / jpg->max_v_samp_factor;
      size_t first_y = y * v_group_dim;
      size_t last_y =
          std::min<size_t>(first_y + v_group_dim, m.height_in_blocks);
      m.context_bits = state.meta[c].context_bits;
      m.context_offset = state.meta[c].context_offset;
      m.ac_coeffs += first_x * brunsli::kDCTBlockSize + first_y * m.ac_stride;
      m.block_state =
          state.meta[c].block_state + first_x + first_y * m.b_stride;
      m.width_in_blocks = last_x - first_x;
      m.height_in_blocks = last_y - first_y;
    }

    status = ProcessJpeg(&ac_state, jpg);
    if (status != BrunsliStatus::BRUNSLI_OK) failed.store(true);
  };
  (*executor)(decode_ac, ac_section_start.size());
  if (failed.load()) return false;

  return true;
}

}  // namespace brunsli
