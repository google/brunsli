// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

package org.brunsli.wrapper;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link org.brunsli.wrapper.Codec}. */
@RunWith(JUnit4.class)
public class CodecTest {

  static {
    String jniLibrary = System.getProperty("BRUNSLI_JNI_LIBRARY");
    if (jniLibrary != null) {
      System.load(new java.io.File(jniLibrary).getAbsolutePath());
    }
  }

  @Test
  public void testRoundtrip() throws IOException {
    byte[] original = {(byte) 0xff, (byte) 0xd8, (byte) 0xff, (byte) 0xe0, 0x00, 0x10, 0x4a, 0x46,
        0x49, 0x46, 0x00, 0x01, 0x01, 0x01, 0x01, 0x2c, 0x01, 0x2c, 0x00, 0x00, (byte) 0xff,
        (byte) 0xdb, 0x00, 0x43, 0x00, 0x02, 0x01, 0x01, 0x01, 0x01, 0x01, 0x02, 0x01, 0x01, 0x01,
        0x02, 0x02, 0x02, 0x02, 0x02, 0x04, 0x03, 0x02, 0x02, 0x02, 0x02, 0x05, 0x04, 0x04, 0x03,
        0x04, 0x06, 0x05, 0x06, 0x06, 0x06, 0x05, 0x06, 0x06, 0x06, 0x07, 0x09, 0x08, 0x06, 0x07,
        0x09, 0x07, 0x06, 0x06, 0x08, 0x0b, 0x08, 0x09, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x06, 0x08,
        0x0b, 0x0c, 0x0b, 0x0a, 0x0c, 0x09, 0x0a, 0x0a, 0x0a, (byte) 0xff, (byte) 0xdb, 0x00, 0x43,
        0x01, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x05, 0x03, 0x03, 0x05, 0x0a, 0x07, 0x06, 0x07,
        0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a,
        0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a,
        0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a,
        0x0a, 0x0a, 0x0a, 0x0a, 0x0a, (byte) 0xff, (byte) 0xc0, 0x00, 0x11, 0x08, 0x00, 0x10, 0x00,
        0x10, 0x03, 0x01, 0x11, 0x00, 0x02, 0x11, 0x01, 0x03, 0x11, 0x01, (byte) 0xff, (byte) 0xc4,
        0x00, 0x15, 0x00, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x07, (byte) 0xff, (byte) 0xc4, 0x00, 0x1f, 0x10, 0x00, 0x02,
        0x01, 0x04, 0x02, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01,
        0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x09, 0x21, 0x00, 0x08, 0x11, (byte) 0xff, (byte) 0xc4,
        0x00, 0x16, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x08, 0x09, (byte) 0xff, (byte) 0xc4, 0x00, 0x21, 0x11, 0x00,
        0x02, 0x00, 0x05, 0x05, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x01, 0x02, 0x03, 0x04, 0x05, 0x11, 0x21, 0x00, 0x06, 0x07, 0x12, 0x13, 0x41, 0x51,
        (byte) 0xff, (byte) 0xda, 0x00, 0x0c, 0x03, 0x01, 0x00, 0x02, 0x11, 0x03, 0x11, 0x00, 0x3f,
        0x00, (byte) 0x9f, 0x71, (byte) 0xe9, (byte) 0xc7, (byte) 0x9e, (byte) 0xd0, (byte) 0xf7,
        0x47, 0x68, 0x63, (byte) 0xf5, 0x35, 0x3a, (byte) 0xfb, 0x24, (byte) 0x97, 0x5a,
        (byte) 0xcb, (byte) 0x92, 0x1b, 0x5e, 0x55, (byte) 0x95, 0x58, 0x5a, 0x24, 0x36,
        (byte) 0xf2, (byte) 0xb1, 0x2c, (byte) 0xae, 0x03, (byte) 0xca, (byte) 0xae, 0x15,
        (byte) 0x82, (byte) 0xc9, 0x11, (byte) 0xed, 0x18, 0x7c, 0x71, (byte) 0xd7, (byte) 0x93,
        (byte) 0xae, (byte) 0xdf, (byte) 0xdb, (byte) 0xf3, 0x55, (byte) 0x99, (byte) 0xa4, 0x25,
        0x1b, (byte) 0xc7, (byte) 0xb5, (byte) 0x99, (byte) 0x85, (byte) 0xb1, (byte) 0x8b,
        (byte) 0xfd, (byte) 0xbf, (byte) 0xe8, (byte) 0xf9, (byte) 0xad, (byte) 0x97, (byte) 0xe5,
        (byte) 0xfe, 0x5f, (byte) 0xa1, (byte) 0xf1, (byte) 0xad, 0x0e, 0x65, 0x56, 0x66, 0x10,
        (byte) 0xa9, 0x08, 0x5e, (byte) 0x90, 0x60, (byte) 0xc4, (byte) 0xec, 0x7b, (byte) 0xdd,
        (byte) 0x8a, (byte) 0x8c, 0x29, 0x52, 0x41, 0x2a, (byte) 0xc3, 0x0c, 0x0d, (byte) 0xc6,
        (byte) 0x9c, (byte) 0x85, (byte) 0xf1, (byte) 0xe7, (byte) 0xb4, 0x3d, 0x2e, (byte) 0xda,
        0x19, 0x05, 0x4d, 0x36, (byte) 0xbe, (byte) 0xc9, 0x22, (byte) 0xd6, (byte) 0xb1, 0x64,
        (byte) 0x82, (byte) 0xd7, (byte) 0x8a, (byte) 0xe5, 0x57, (byte) 0xe6, (byte) 0x89,
        (byte) 0xcd, (byte) 0xc0, (byte) 0xb4, 0x4d, 0x2a, 0x02, (byte) 0xf1, 0x2a, 0x06, 0x62,
        (byte) 0xb1, (byte) 0xca, 0x7a, 0x45, 0x1f, 0x10, (byte) 0xf5, (byte) 0xe3, 0x70, 0x6d,
        (byte) 0xf9, (byte) 0xaa, 0x34, (byte) 0xd3, (byte) 0x90, (byte) 0x8d, (byte) 0xe3,
        (byte) 0xda, (byte) 0xca, (byte) 0xc6, (byte) 0xd9, (byte) 0xc5, (byte) 0xfe, 0x5b,
        (byte) 0xf0, (byte) 0xfc, (byte) 0xd3, (byte) 0x88, 0x39, 0x7e, (byte) 0x87, (byte) 0xc9,
        0x54, 0x39, 0x65, 0x69, (byte) 0x98, 0x46, (byte) 0xa4, 0x61, 0x7a, 0x46, (byte) 0x83, 0x0f,
        (byte) 0xb0, (byte) 0xe9, 0x66, 0x0a, 0x70, (byte) 0xc5, (byte) 0x88, 0x00, (byte) 0xb2,
        (byte) 0x8c, (byte) 0xb1, 0x37, 0x3a, 0x71, (byte) 0xe9, (byte) 0xc8, 0x66, (byte) 0xd0,
        (byte) 0xf4, (byte) 0xbb, 0x68, 0x63, (byte) 0xf4, (byte) 0xd5, 0x3b, 0x07, 0x24,
        (byte) 0x8b, 0x5a, (byte) 0xc5, (byte) 0x92, 0x1b, (byte) 0xa6, 0x55, (byte) 0x8a,
        (byte) 0xd8, 0x56, 0x27, 0x37, 0x02, (byte) 0xd1, 0x2c, 0x4e, 0x42, 0x4a, (byte) 0xc8, 0x19,
        (byte) 0x8a, (byte) 0xc7, 0x10, (byte) 0xed, (byte) 0xd4, 0x7c, 0x41, (byte) 0xdf,
        (byte) 0x8d, (byte) 0xbf, (byte) 0xb8, 0x26, (byte) 0xa8, (byte) 0xd3, 0x48, 0x0b,
        (byte) 0xb7, (byte) 0x8f, 0x6b, (byte) 0xb2, (byte) 0x8b, 0x67, 0x16, (byte) 0xfb, 0x6f,
        (byte) 0xc1, (byte) 0xf7, 0x4e, 0x5f, (byte) 0xe2, 0x0a, 0x1f, 0x25, 0x50, (byte) 0xe6,
        0x59, 0x65, (byte) 0xa1, 0x1a, (byte) 0x91, (byte) 0x85, (byte) 0xe7, 0x06, 0x34, 0x4e,
        (byte) 0xc3, (byte) 0xa5, (byte) 0x98, (byte) 0xb0, (byte) 0xca, (byte) 0x86, 0x20, 0x02,
        (byte) 0xcc, 0x70, (byte) 0xa4, (byte) 0xdc, (byte) 0xe9, (byte) 0xc8, 0x5f, 0x21,
        (byte) 0x9b, 0x43, (byte) 0xdd, 0x1d, (byte) 0xa1, (byte) 0x90, 0x53, 0x53, 0x6c, 0x1c,
        (byte) 0x92, 0x5d, 0x6b, 0x2e, 0x48, 0x2e, (byte) 0x98, (byte) 0xae, 0x2b, 0x7e, 0x58,
        (byte) 0x90, (byte) 0xdb, (byte) 0xca, (byte) 0xc4, (byte) 0xd1, 0x21, 0x29, 0x13, 0x38,
        0x56, 0x0b, 0x24, (byte) 0xa3, (byte) 0xa7, 0x61, (byte) 0xf1, (byte) 0xcf, 0x7e, 0x37,
        0x06, (byte) 0xe0, (byte) 0x9a, (byte) 0xac, (byte) 0xcd, 0x38, 0x0e, (byte) 0xde, 0x3d,
        (byte) 0xae, (byte) 0xaa, 0x6d, (byte) 0x8c, 0x5b, (byte) 0xe5, (byte) 0xff, 0x00, 0x4f,
        (byte) 0xdd, 0x38, (byte) 0x83, (byte) 0x88, 0x28, 0x7c, 0x6b, 0x43, (byte) 0x96, 0x66,
        (byte) 0x96, (byte) 0x84, 0x2a, 0x42, 0x17, (byte) 0x9c, 0x68, (byte) 0xd0, (byte) 0xfb,
        0x1e, (byte) 0xf7, 0x60, (byte) 0xc7, 0x2c, 0x14, (byte) 0x90, 0x4a, (byte) 0xa9,
        (byte) 0xca, (byte) 0x83, 0x71, (byte) 0xaf, (byte) 0xff, (byte) 0xd9};
    byte[] compressed = Codec.encode(original);
    assertTrue(compressed.length < original.length);
    byte[] decompressed = Codec.decode(compressed);
    assertArrayEquals(original, decompressed);
  }
}
