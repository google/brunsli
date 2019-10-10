// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

package dev.brunsli.wrapper;

import java.nio.ByteBuffer;

/**
 * Native functions declaration.
 */
public class CodecJNI {
  static native byte[] nativeDecode(ByteBuffer data);
  static native byte[] nativeEncode(ByteBuffer data);

  private CodecJNI() {}
}
