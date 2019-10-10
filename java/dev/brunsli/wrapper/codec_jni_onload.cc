// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.


#include <jni.h>

#include "codec_jni.h"

#ifdef __cplusplus
extern "C" {
#endif

static const JNINativeMethod kCodecMethods[] = {
    {"nativeDecode", "(Ljava/nio/ByteBuffer;)[B",
     reinterpret_cast<void*>(Java_dev_brunsli_wrapper_CodecJNI_nativeDecode)},
    {"nativeEncode", "(Ljava/nio/ByteBuffer;)[B",
     reinterpret_cast<void*>(Java_dev_brunsli_wrapper_CodecJNI_nativeEncode)}};

JNIEXPORT jint JNI_OnLoad(JavaVM* vm, void* reserved) {
  JNIEnv* env;
  if (vm->GetEnv(reinterpret_cast<void**>(&env), JNI_VERSION_1_6) != JNI_OK) {
    return -1;
  }

  jclass clazz = env->FindClass("dev/brunsli/wrapper/CodecJNI");
  if (clazz == nullptr) {
    return -1;
  }

  int num_methods = sizeof(kCodecMethods) / sizeof(kCodecMethods[0]);
  jint result = env->RegisterNatives(clazz, kCodecMethods, num_methods);
  if (result < 0) {
    return -1;
  }

  return JNI_VERSION_1_6;
}

#ifdef __cplusplus
}
#endif
