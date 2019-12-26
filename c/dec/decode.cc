// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>

#include <brunsli/decode.h>
#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>
#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>

#include "groups.h"
#include "highwayhash/data_parallel.h"

/* C API for brunsli encoder */

extern "C" {

struct OutputStruct {
  size_t (*fun)(void* out_data, const uint8_t* buf, size_t size);
  void* data;
};

int StringWriter(void* data, const uint8_t* buf, size_t count) {
  std::string* output = reinterpret_cast<std::string*>(data);
  output->append(reinterpret_cast<const char*>(buf), count);
  return count;
}

int DecodeBrunsli(size_t in_size, const uint8_t* in, void* out_data,
                  DecodeBrunsliSink out_fun) {
  OutputStruct out = {out_fun, out_data};
  brunsli::JPEGData jpg;
  brunsli::BrunsliStatus status = brunsli::BrunsliDecodeJpeg(in, in_size, &jpg);
  if (status != brunsli::BRUNSLI_OK) {
    return 0;
  }
  brunsli::JPEGOutput writer(
      [](void* data, const uint8_t* buf, size_t count) {
        OutputStruct* sink = (OutputStruct*)data;
        int result = sink->fun(sink->data, buf, count) == count ? count : -1;
        return result;
      },
      &out);
  if (!brunsli::WriteJpeg(jpg, writer)) {
    return 0;
  }
  return 1;  // ok
}

size_t brunsli_decompress(const char* input, size_t input_size, std::string& output, size_t thread_nums) {
  brunsli::JPEGData jpg;
  const uint8_t* input_data = reinterpret_cast<const uint8_t*>(input);

  // 1. decode
  highwayhash::ThreadPool thread_pool(thread_nums);
  brunsli::Executor executor = [&](const brunsli::Runnable& runnable, size_t num_tasks) {
    thread_pool.Run(0, num_tasks, runnable);
  };
  bool ok = brunsli::DecodeGroups(input_data, input_size, &jpg, 32, 128, &executor);
  if (!ok) {
    fprintf(stderr, "Failed to parse Brunsli input.\n");
    return false;
  }

  //2. write jpeg
  brunsli::JPEGOutput writer(StringWriter, &output);
  ok = brunsli::WriteJpeg(jpg, writer);
  if (!ok) {
    fprintf(stderr, "Failed to serialize JPEG data.\n");
    return false;
  }

  return output.size();
}

} /* extern "C" */
