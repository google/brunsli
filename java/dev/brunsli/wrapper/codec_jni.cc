// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "codec_jni.h"

#include <cstring>
#include <string>

#include <brunsli/brunsli_decode.h>
#include <brunsli/jpeg_data_writer.h>
#include <brunsli/brunsli_encode.h>
#include <brunsli/jpeg_data_reader.h>

namespace {

jbyteArray ToByteArray(JNIEnv* env, const std::string& data) {
  jbyteArray result = env->NewByteArray(data.size());
  if (!result) return nullptr;
  env->SetByteArrayRegion(result, 0, data.size(),
                          reinterpret_cast<const jbyte*>(data.data()));
  return result;
}

}  // namespace

#ifdef __cplusplus
extern "C" {
#endif

JNIEXPORT jbyteArray JNICALL Java_dev_brunsli_wrapper_CodecJNI_nativeDecode(
    JNIEnv* env, jobject /*jobj*/, jobject input) {
  uint8_t* data =
      reinterpret_cast<uint8_t*>(env->GetDirectBufferAddress(input));
  if (!data) return nullptr;
  jlong length = env->GetDirectBufferCapacity(input);
  if (length == -1) return nullptr;

  ::brunsli::JPEGData jpg;
  ::brunsli::BrunsliStatus status = ::brunsli::BrunsliDecodeJpeg(
      data, length, &jpg);
  if (status != ::brunsli::BRUNSLI_OK) {
    return nullptr;
  }

  std::string output;
  auto write = [](void* data, const uint8_t* buf, size_t count) -> size_t {
    std::string* output = reinterpret_cast<std::string*>(data);
    output->append(reinterpret_cast<const char*>(buf), count);
    return count;
  };
  ::brunsli::JPEGOutput writer(write, &output);
  if (!::brunsli::WriteJpeg(jpg, writer)) {
    return nullptr;
  }

  return ToByteArray(env, output);
}

JNIEXPORT jbyteArray JNICALL Java_dev_brunsli_wrapper_CodecJNI_nativeEncode(
    JNIEnv* env, jobject /*jobj*/, jobject input) {
  uint8_t* data =
      reinterpret_cast<uint8_t*>(env->GetDirectBufferAddress(input));
  if (!data) return nullptr;
  jlong length = env->GetDirectBufferCapacity(input);
  if (length == -1) return nullptr;

  ::brunsli::JPEGData jpg;
  if (!::brunsli::ReadJpeg(data, length, ::brunsli::JPEG_READ_ALL, &jpg)) {
    return nullptr;
  }
  size_t output_length = ::brunsli::GetMaximumBrunsliEncodedSize(jpg);
  std::string output;
  output.resize(output_length);
  bool ok = ::brunsli::BrunsliEncodeJpeg(
      jpg, reinterpret_cast<uint8_t*>(&output[0]), &output_length);
  if (!ok) return nullptr;
  output.resize(output_length);

  return ToByteArray(env, output);
}

#ifdef __cplusplus
}
#endif
