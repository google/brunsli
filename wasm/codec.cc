// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include <string>

#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>
#include <brunsli/brunsli_encode.h>
#include <brunsli/jpeg_data_reader.h>

extern "C" {

std::string* BrunsliToJpeg(const uint8_t* data, size_t length) {
  std::string* result = nullptr;
  brunsli::JPEGData jpg;
  brunsli::BrunsliStatus status =
      brunsli::BrunsliDecodeJpeg(data, length, &jpg);
  if (status != brunsli::BRUNSLI_OK) {
    printf("Decoding Brunsli failed with status %d\n", status);
    return result;
  }

  result = new std::string();
  auto write = [](void* data, const uint8_t* buf, size_t count) -> size_t {
    std::string* result = reinterpret_cast<std::string*>(data);
    result->append(reinterpret_cast<const char*>(buf), count);
    return count;
  };
  brunsli::JPEGOutput writer(write, result);
  if (!brunsli::WriteJpeg(jpg, writer)) {
    printf("Serializing JPEG failed\n");
    delete result;
    result = nullptr;
  }
  return result;
}

/*
 * Instance layout (uint32 array):
 *  0: decoder
 *  1: max_in_len
 *  2: in
 *  3: in_len - the only writable field
 *  4: out
 *  5: out_len
 */

static const size_t kBufferSize = 65536;

uint32_t* BrunsliDecoderInit() {
  uint32_t* instance = new uint32_t[6];
  // TODO(eustas): check for OOMs
  brunsli::BrunsliDecoder* decoder = new brunsli::BrunsliDecoder();
  uint8_t* in = reinterpret_cast<uint8_t*>(malloc(kBufferSize));
  uint8_t* out = reinterpret_cast<uint8_t*>(malloc(kBufferSize));
  instance[0] = reinterpret_cast<uint32_t>(decoder);
  instance[1] = kBufferSize;
  instance[2] = reinterpret_cast<uint32_t>(in);
  instance[3] = 0;
  instance[4] = reinterpret_cast<uint32_t>(out);
  instance[5] = 0;
  return instance;
}

/*
 * Parse input and produce output.
 *
 * Input should be passed via `in_len`; data have to be stored starting at `in`.
 * It is considered that all input is consumed.
 * After invocation `out_len` output bytes are placed at `out`.
 * It is considered that client consumes all the output right after inovcation.
 *
 * Result:
 *  0: ok - keep feeding input / consuming output
 *  1: done - no more input / output
 *  2: error
 */
uint32_t BrunsliDecoderProcess(uint32_t* instance) {
  brunsli::BrunsliDecoder* decoder =
      reinterpret_cast<brunsli::BrunsliDecoder*>(instance[0]);
  const uint8_t* next_in = reinterpret_cast<uint8_t*>(instance[2]);
  size_t available_in = instance[3];
  uint8_t* next_out = reinterpret_cast<uint8_t*>(instance[4]);
  size_t available_out = kBufferSize;
  brunsli::BrunsliDecoder::Status result = decoder->Decode(
      &available_in, &next_in, &available_out, &next_out);
  instance[5] = kBufferSize - available_out;
  if ((result == brunsli::BrunsliDecoder::NEEDS_MORE_INPUT) ||
      (result == brunsli::BrunsliDecoder::NEEDS_MORE_OUTPUT)) {
    return 0;
  }
  if (result == brunsli::BrunsliDecoder::DONE) return 1;
  return 2;
}

void BrunsliDecoderCleanup(uint32_t* instance) {
  if (instance == nullptr) return;
  brunsli::BrunsliDecoder* decoder =
      reinterpret_cast<brunsli::BrunsliDecoder*>(instance[0]);
  delete decoder;
  uint8_t* in = reinterpret_cast<uint8_t*>(instance[2]);
  delete[] in;
  uint8_t* out = reinterpret_cast<uint8_t*>(instance[4]);
  delete[] out;
  delete[] instance;
}

const char* GetJpegData(std::string* jpeg) {
  if (!jpeg) return 0;
  return jpeg->data();
}

size_t GetJpegLength(std::string* jpeg) {
  if (!jpeg) return 0;
  return jpeg->size();
}

void FreeJpeg(std::string* jpeg) { delete jpeg; }

std::string* JpegToBrunsli(const uint8_t* data, size_t length) {
  std::string* result = nullptr;
  brunsli::JPEGData jpg;
  if (!brunsli::ReadJpeg(data, length, brunsli::JPEG_READ_ALL, &jpg)) {
    printf("Parsing JPEG failed\n");
    return result;
  }
  size_t result_length = brunsli::GetMaximumBrunsliEncodedSize(jpg);
  result = new std::string();
  result->resize(result_length);
  bool ok = brunsli::BrunsliEncodeJpeg(
      jpg, reinterpret_cast<uint8_t*>(&result->at(0)), &result_length);
  if (ok) {
    result->resize(result_length);
  } else {
    printf("Encoding Brunsli failed\n");
    delete result;
    result = nullptr;
  }
  return result;
}

const char* GetBrunsliData(std::string* brunsli) {
  if (!brunsli) return 0;
  return brunsli->data();
}

size_t GetBrunsliLength(std::string* brunsli) {
  if (!brunsli) return 0;
  return brunsli->size();
}

void FreeBrunsli(std::string* brunsli) { delete brunsli; }

}  // extern "C"
