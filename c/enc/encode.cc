// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>

#include <brunsli/encode.h>
#include <brunsli/jpeg_data.h>
#include <brunsli/types.h>
#include <brunsli/brunsli_encode.h>
#include <brunsli/jpeg_data_reader.h>
#include "groups.h"
#include "highwayhash/data_parallel.h"

/* C API for brunsli encoder */

extern "C" {

int EncodeBrunsli(size_t insize, const unsigned char* in, void* outdata,
    size_t (*outfun)(void* outdata, const unsigned char* buf, size_t size)) {

  std::vector<uint8_t> output;
  brunsli::JPEGData jpg;
  if (!brunsli::ReadJpeg(in, insize, brunsli::JPEG_READ_ALL, &jpg)) {
    return 0;
  }
  size_t output_size = brunsli::GetMaximumBrunsliEncodedSize(jpg);
  output.resize(output_size);
  if (!brunsli::BrunsliEncodeJpeg(jpg, output.data(), &output_size)) {
    return 0;
  }
  output.resize(output_size);

  if (!outfun(outdata, reinterpret_cast<const unsigned char*>(output.data()),
      output.size())) {
    return 0;
  }
  return 1;  /* ok */
}

size_t brunsli_compress(const char* input, size_t input_size, std::string& output, size_t thread_nums) {
  brunsli::JPEGData jpg;
  const uint8_t* input_data = reinterpret_cast<const uint8_t*>(input);

  //1. read jpeg
  bool ok  = brunsli::ReadJpeg(input_data, input_size, brunsli::JPEG_READ_ALL, &jpg);
  if (!ok) {
    fprintf(stderr, "Failed to parse JPEG input.\n");
    return false;
  }

  //2. brunsli encode
  size_t output_size = brunsli::GetMaximumBrunsliEncodedSize(jpg);
  output.resize(output_size);
  uint8_t* output_data = reinterpret_cast<uint8_t*>(&output[0]);

  highwayhash::ThreadPool thread_pool(thread_nums);
  brunsli::Executor executor = [&](const brunsli::Runnable& runnable, size_t num_tasks) {
      thread_pool.Run(0, num_tasks, runnable);
  };
  ok = brunsli::EncodeGroups(jpg, output_data, &output_size, 32, 128, &executor);
  if (!ok) {
    fprintf(stderr, "Failed to transform JPEG to Brunsli\n");
    return false;
  }
  output.resize(output_size);
  
  return output_size;
}

}  /* extern "C" */