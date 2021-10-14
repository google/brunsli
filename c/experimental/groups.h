// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Functions for writing encoding / decoding Brunsli in "groups" mode.

#ifndef BRUNSLI_EXPERIMENTAL_GROUPS_H_
#define BRUNSLI_EXPERIMENTAL_GROUPS_H_

#include <functional>
#include <future>  // NOLINT(build/c++11)
#include <mutex>  // NOLINT(build/c++11)

#include <brunsli/jpeg_data.h>
#include <brunsli/types.h>

namespace brunsli {

typedef std::function<void(size_t)> Runnable;
typedef std::function<void(const Runnable&, size_t)> Executor;

void SequentialExecutor(const Runnable& runnable, size_t num_tasks);

// Poor man's thread pool.
class ParallelExecutor {
 public:
  explicit ParallelExecutor(size_t num_threads);
  ~ParallelExecutor();
  Executor getExecutor();

 private:
  void execute(const Runnable& runnable, size_t num_tasks);
  size_t num_threads;
  std::mutex lock;
  std::condition_variable start_latch;
  std::condition_variable finish_latch;
  size_t num_tasks = 0;
  const Runnable* runnable = nullptr;
  std::atomic<size_t> next_task{0};
  std::atomic<size_t> busy_count{0};
  bool terminate = false;
  std::vector<std::future<void>> futures;
};

bool DecodeGroups(const uint8_t* data, size_t len, brunsli::JPEGData* jpg,
                  size_t ac_group_dim, size_t dc_group_dim, Executor* executor);

bool EncodeGroups(const brunsli::JPEGData& jpg, uint8_t* data, size_t* len,
                  size_t ac_group_dim, size_t dc_group_dim, Executor* executor);

}  // namespace brunsli

#endif  // BRUNSLI_EXPERIMENTAL_GROUPS_H_
