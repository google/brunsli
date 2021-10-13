// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <cstdio>
#include <cstdlib>
#include <string>

#include <brunsli/jpeg_data.h>
#include <brunsli/status.h>
#include <brunsli/types.h>
#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>

#if defined(BRUNSLI_EXPERIMENTAL_GROUPS)
#include "../experimental/groups.h"
#endif

#if defined(_WIN32)
#define fopen ms_fopen
static FILE* ms_fopen(const char* filename, const char* mode) {
  FILE* result = 0;
  fopen_s(&result, filename, mode);
  return result;
}
#endif  /* WIN32 */

size_t StringWriter(void* data, const uint8_t* buf, size_t count) {
  std::string* output = reinterpret_cast<std::string*>(data);
  output->append(reinterpret_cast<const char*>(buf), count);
  return count;
}

bool ReadFileInternal(FILE* file, std::string* content) {
  if (fseek(file, 0, SEEK_END) != 0) {
    fprintf(stderr, "Failed to seek end of input file.\n");
    return false;
  }
  int input_size = ftell(file);
  if (input_size == 0) {
    fprintf(stderr, "Input file is empty.\n");
    return false;
  }
  if (fseek(file, 0, SEEK_SET) != 0) {
    fprintf(stderr, "Failed to rewind input file to the beginning.\n");
    return false;
  }
  content->resize(input_size);
  size_t read_pos = 0;
  while (read_pos < content->size()) {
    const size_t bytes_read =
        fread(&content->at(read_pos), 1, content->size() - read_pos, file);
    if (bytes_read == 0) {
      fprintf(stderr, "Failed to read input file\n");
      return false;
    }
    read_pos += bytes_read;
  }
  return true;
}

bool ReadFile(const std::string& file_name, std::string* content) {
  FILE* file = fopen(file_name.c_str(), "rb");
  if (file == nullptr) {
    fprintf(stderr, "Failed to open input file.\n");
    return false;
  }
  bool ok = ReadFileInternal(file, content);
  if (fclose(file) != 0) {
    if (ok) {
      fprintf(stderr, "Failed to close input file.\n");
    }
    return false;
  }
  return ok;
}

bool WriteFileInternal(FILE* file, const std::string& content) {
  size_t write_pos = 0;
  while (write_pos < content.size()) {
    const size_t bytes_written =
        fwrite(&content[write_pos], 1, content.size() - write_pos, file);
    if (bytes_written == 0) {
      fprintf(stderr, "Failed to write output.\n");
      return false;
    }
    write_pos += bytes_written;
  }
  return true;
}

bool WriteFile(const std::string& file_name, const std::string& content) {
  FILE* file = fopen(file_name.c_str(), "wb");
  if (file == nullptr) {
    fprintf(stderr, "Failed to open file for writing.\n");
    return false;
  }
  bool ok = WriteFileInternal(file, content);
  if (fclose(file) != 0) {
    if (ok) {
      fprintf(stderr, "Failed to close output file.\n");
    }
    return false;
  }
  return ok;
}

bool ProcessFile(const std::string& file_name,
                 const std::string& outfile_name) {
  std::string input;
  bool ok = ReadFile(file_name, &input);
  if (!ok) return false;

  std::string output;
  {
    brunsli::JPEGData jpg;
    const uint8_t* input_data = reinterpret_cast<const uint8_t*>(input.data());

#if defined(BRUNSLI_EXPERIMENTAL_GROUPS)
    {
      brunsli::ParallelExecutor pool(4);
      brunsli::Executor executor = pool.getExecutor();
      ok = brunsli::DecodeGroups(input_data, input.size(), &jpg, 32, 128,
                                 &executor);
    }
#else
    brunsli::BrunsliStatus status =
        brunsli::BrunsliDecodeJpeg(input_data, input.size(), &jpg);
    ok = (status == brunsli::BRUNSLI_OK);
#endif

    input.clear();
    input.shrink_to_fit();
    if (!ok) {
      fprintf(stderr, "Failed to parse Brunsli input.\n");
      return false;
    }

    brunsli::JPEGOutput writer(StringWriter, &output);
    ok = brunsli::WriteJpeg(jpg, writer);
    if (!ok) {
      fprintf(stderr, "Failed to serialize JPEG data.\n");
      return false;
    }
  }

  ok = WriteFile(outfile_name, output);

  return ok;
}

int main(int argc, char** argv) {
  if (argc != 2 && argc != 3) {
    fprintf(stderr, "Usage: dbrunsli FILE [OUTPUT_FILE, default=FILE.jpg]\n");
    return EXIT_FAILURE;
  }
  const std::string file_name = std::string(argv[1]);
  if (file_name.empty()) {
    fprintf(stderr, "Empty input file name.\n");
    return EXIT_FAILURE;
  }
  const std::string outfile_name =
      argc == 2 ? file_name + ".jpg" : std::string(argv[2]);

  bool ok = ProcessFile(file_name, outfile_name);
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
