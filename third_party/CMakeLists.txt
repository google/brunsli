# Copyright (c) Google LLC 2019
#
# Use of this source code is governed by an MIT-style
# license that can be found in the LICENSE file or at
# https://opensource.org/licenses/MIT.

# brotli
if (NOT TARGET brotlicommon-static)
  if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/brotli/c/include/brotli/decode.h")
    message(SEND_ERROR "Brotli is required to build brunsli, run \"git submodule update --init\" to get it, then run cmake again")
  endif ()
  add_subdirectory(brotli EXCLUDE_FROM_ALL)
endif ()

if (NOT TARGET gtest_main)
  add_subdirectory(googletest EXCLUDE_FROM_ALL)
endif ()
