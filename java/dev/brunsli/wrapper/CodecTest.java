// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

package dev.brunsli.wrapper;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.Base64;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link dev.brunsli.wrapper.Codec}. */
@RunWith(JUnit4.class)
public class CodecTest {
  static {
    String jniLibrary = System.getProperty("BRUNSLI_JNI_LIBRARY");
    if (jniLibrary != null) {
      System.load(new java.io.File(jniLibrary).getAbsolutePath());
    }
  }

  // clang-format off
  private static final String ORIG_BASE64 =
      "/9j/4AAQSkZJRgABAQEBLAEsAAD/2wBDAAIBAQEBAQIBAQECAgICAgQDAgICAgUEBAMEBgUG"
          + "BgYFBgYGBwkIBgcJBwYGCAsICQoKCgoKBggLDAsKDAkKCgr/2wBDAQICAgICAgUDAw"
          + "UKBwYHCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoK"
          + "CgoKCgr/wAARCAAQABADAREAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAf/xA"
          + "AfEAACAQQCAwAAAAAAAAAAAAABAgMEBQYHCSEACBH/xAAWAQEBAQAAAAAAAAAAAAAA"
          + "AAAACAn/xAAhEQACAAUFAQEAAAAAAAAAAAABAgMEBREhAAYHEhNBUf/aAAwDAQACEQ"
          + "MRAD8An3Hpx57Q90doY/U1Ovskl1rLkhteVZVYWiQ28rEsrgPKrhWCyRHtGHxx15Ou"
          + "39vzVZmkJRvHtZmFsYv9v+j5rZfl/l+h8a0OZVZmEKkIXpBgxOx73YqMKVJBKsMMDc"
          + "achfHntD0u2hkFTTa+ySLWsWSC14rlV+aJzcC0TSoC8SoGYrHKekUfEPXjcG35qjTT"
          + "kI3j2srG2cX+W/D804g5fofJVDllaZhGpGF6RoMPsOlmCnDFiACyjLE3OnHpyGbQ9L"
          + "toY/TVOwcki1rFkhumVYrYVic3AtEsTkJKyBmKxxDt1HxB342/uCao00gLt49rsotn"
          + "FvtvwfdOX+IKHyVQ5llloRqRhecGNE7DpZiwyoYgAsxwpNzpyF8hm0PdHaGQU1NsHJ"
          + "Jday5ILpiuK35YkNvKxNEhKRM4Vgsko6dh8c9+NwbgmqzNOA7ePa6qbYxb5f8AT904"
          + "g4gofGtDlmaWhCpCF5xo0Pse92DHLBSQSqnKg3Gv/9k=";

  // clang-format on

  @Test
  public void testRoundtrip() throws IOException {
    byte[] original = Base64.getDecoder().decode(ORIG_BASE64);
    byte[] compressed = Codec.encode(original);
    assertTrue(compressed.length < original.length);
    byte[] decompressed = Codec.decode(compressed);
    assertArrayEquals(original, decompressed);
  }
}
