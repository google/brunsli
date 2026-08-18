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

  @Test
  public void testRoundtrip() throws IOException {
    // clang-format off
    String orig64 =
        "/9j/4AAQSkZJRgABAQEBLAEsAAD/2wBDAAIBAQEBAQIBAQECAgICAgQD"
            + "AgICAgUEBAMEBgUGBgYFBgYGBwkIBgcJBwYGCAsICQoKCgoKBggLDAsKDAkKCgr/2wBDAQICAgIC"
            + "AgUDAwUKBwYHCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoKCgoK"
            + "Cgr/wAARCAAQABADAREAAhEBAxEB/8QAFQABAQAAAAAAAAAAAAAAAAAAAAf/xAAfEAACAQQCAwAA"
            + "AAAAAAAAAAABAgMEBQYHCSEACBH/xAAWAQEBAQAAAAAAAAAAAAAAAAAACAn/xAAhEQACAAUFAQEA"
            + "AAAAAAAAAAABAgMEBREhAAYHEhNBUf/aAAwDAQACEQMRAD8An3Hpx57Q90doY/U1Ovskl1rLkhte"
            + "VZVYWiQ28rEsrgPKrhWCyRHtGHxx15Ou39vzVZmkJRvHtZmFsYv9v+j5rZfl/l+h8a0OZVZmEKkI"
            + "XpBgxOx73YqMKVJBKsMMDcachfHntD0u2hkFTTa+ySLWsWSC14rlV+aJzcC0TSoC8SoGYrHKekUf"
            + "EPXjcG35qjTTkI3j2srG2cX+W/D804g5fofJVDllaZhGpGF6RoMPsOlmCnDFiACyjLE3OnHpyGbQ"
            + "9LtoY/TVOwcki1rFkhumVYrYVic3AtEsTkJKyBmKxxDt1HxB342/uCao00gLt49rsotnFvtvwfdO"
            + "X+IKHyVQ5llloRqRhecGNE7DpZiwyoYgAsxwpNzpyF8hm0PdHaGQU1NsHJJday5ILpiuK35YkNvK"
            + "xNEhKRM4Vgsko6dh8c9+NwbgmqzNOA7ePa6qbYxb5f8AT904g4gofGtDlmaWhCpCF5xo0Pse92DH"
            + "LBSQSqnKg3Gv/9k=";
    // clang-format on
    byte[] original = Base64.getDecoder().decode(orig64);
    byte[] compressed = Codec.encode(original);
    assertTrue(compressed.length < original.length);
    byte[] decompressed = Codec.decode(compressed);
    assertArrayEquals(original, decompressed);
  }
}
