// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

package dev.brunsli.wrapper;

import java.io.IOException;
import java.nio.ByteBuffer;

/**
 * Brunsli codec wrapper.
 */
public class Codec {
  /**
   * Decodes Brunsli stream to JPEG stream.
   *
   * Input data is accessed without copying.
   *
   * @param data direct ByteBuffer containing Brunsli stream
   * @return JPEG stream
   */
  public static byte[] decode(ByteBuffer data) throws IOException {
    if (!data.isDirect()) {
      throw new IllegalArgumentException("only direct buffers allowed");
    }
    byte[] result = CodecJNI.nativeDecode(data);
    if (result == null) {
      throw new IOException("Brunsli decoder failed");
    }
    return result;
  }

  /**
   * Decodes Brunsli stream to JPEG stream.
   *
   * Input is wrapped (copied) into a direct ByteBuffer
   *
   * @param data Brunsli stream
   * @return JPEG stream
   */
  public static byte[] decode(byte[] data) throws IOException {
    ByteBuffer input = ByteBuffer.allocateDirect(data.length);
    input.put(data);
    return decode(input);
  }

  /**
   * Encodes JPEG stream to Brunsli stream.
   *
   * Input data is accessed without copying.
   *
   * @param data direct ByteBuffer containing JPEG stream
   * @return Brunsli stream
   */
  public static byte[] encode(ByteBuffer data) throws IOException {
    if (!data.isDirect()) {
      throw new IllegalArgumentException("only direct buffers allowed");
    }
    byte[] result = CodecJNI.nativeEncode(data);
    if (result == null) {
      throw new IOException("Brunsli encoder failed");
    }
    return result;
  }

  /**
   * Encodes JPEG stream to Brunsli stream.
   *
   * Input is wrapped (copied) into a direct ByteBuffer
   *
   * @param data JPEG stream
   * @return Brunsli stream
   */
  public static byte[] encode(byte[] data) throws IOException {
    ByteBuffer input = ByteBuffer.allocateDirect(data.length);
    input.put(data);
    return encode(input);
  }

  private Codec() {}
}
