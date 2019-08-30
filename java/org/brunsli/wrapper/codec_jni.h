// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#ifndef ORG_BRUNSLI_WRAPPER_CODEC_JNI_H_
#define ORG_BRUNSLI_WRAPPER_CODEC_JNI_H_

#include <jni.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Decodes Brunsli stream to JPEG stream.
 *
 * @param data direct ByteBuffer containing Brunsli stream
 * @returns byte array JPEG stream on success; otherwise null
 */
JNIEXPORT jbyteArray JNICALL
Java_org_brunsli_wrapper_CodecJNI_nativeDecode(
    JNIEnv* env, jobject /*jobj*/, jobject input);

/**
 * Encodes JPEG stream to Brunsli stream.
 *
 * @param data direct ByteBuffer containing JPEG stream
 * @returns byte array Brunsli stream on success; otherwise null
 */
JNIEXPORT jbyteArray JNICALL
Java_org_brunsli_wrapper_CodecJNI_nativeEncode(
    JNIEnv* env, jobject /*jobj*/, jobject input);

#ifdef __cplusplus
}
#endif

#endif  // ORG_BRUNSLI_WRAPPER_CODEC_JNI_H_
