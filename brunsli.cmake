# Copyright (c) Google LLC 2019
#
# Use of this source code is governed by an MIT-style
# license that can be found in the LICENSE file or at
# https://opensource.org/licenses/MIT.

file(GLOB BRUNSLI_COMMON_SOURCES
  c/common/*.cc
)

# TODO(eustas): split public/private headers.
file(GLOB BRUNSLI_COMMON_HEADERS
  c/common/*.h
)

set(BRUNSLI_DEC_SOURCES
  c/dec/ans_decode.cc
  c/dec/bit_reader.cc
  c/dec/brunsli_decode.cc
  c/dec/context_map_decode.cc
  c/dec/histogram_decode.cc
  c/dec/huffman_decode.cc
  c/dec/huffman_table.cc
  c/dec/jpeg_data_writer.cc
  c/dec/state.cc
)

# TODO(eustas): split public/private headers.
file(GLOB BRUNSLI_DEC_HEADERS
  c/dec/*.h
)

set(BRUNSLI_ENC_SOURCES
  c/enc/ans_encode.cc
  c/enc/brunsli_encode.cc
  c/enc/context_map_encode.cc
  c/enc/histogram_encode.cc
  c/enc/huffman_encode.cc
  c/enc/huffman_tree.cc
  c/enc/jpeg_data_reader.cc
  c/enc/jpeg_huffman_decode.cc
  c/enc/write_bits.cc
)

# TODO(eustas): split public/private headers.
file(GLOB BRUNSLI_ENC_HEADERS
  c/enc/*.h
)

set(BRUNSLI_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/c/include")
mark_as_advanced(BRUNSLI_INCLUDE_DIRS)

add_library(brunslicommon-static STATIC
  ${BRUNSLI_COMMON_SOURCES}
  ${BRUNSLI_COMMON_HEADERS}
)

add_library(brunslidec-static STATIC
  ${BRUNSLI_DEC_SOURCES}
  ${BRUNSLI_DEC_HEADERS}
)
target_link_libraries(brunslidec-static PRIVATE
  brotlidec-static
  brunslicommon-static
)

add_library(brunslienc-static STATIC
  ${BRUNSLI_ENC_SOURCES}
  ${BRUNSLI_ENC_HEADERS}
)
target_link_libraries(brunslienc-static PRIVATE
  brotlienc-static
  brunslicommon-static
)

set(BRUNSLI_LIBRARIES brunslicommon-static brunslidec-static brunslienc-static)

if(NOT BRUNSLI_EMSCRIPTEN)
add_library(brunslidec-c SHARED
  c/dec/decode.cc
)
target_link_libraries(brunslidec-c PRIVATE brunslidec-static)
add_library(brunslienc-c SHARED
  c/enc/encode.cc
)
target_link_libraries(brunslienc-c PRIVATE brunslienc-static)
list(APPEND BRUNSLI_LIBRARIES brunslidec-c brunslienc-c)
endif()  # BRUNSLI_EMSCRIPTEN

foreach(lib IN LISTS BRUNSLI_LIBRARIES)
  target_include_directories(${lib} PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/c/include"
    "${CMAKE_CURRENT_SOURCE_DIR}"
  )
  set_property(TARGET ${lib} PROPERTY POSITION_INDEPENDENT_CODE ON)
endforeach()

add_executable(cbrunsli c/tools/cbrunsli.cc)
target_link_libraries(cbrunsli PRIVATE
  brunslienc-static
)
add_executable(dbrunsli c/tools/dbrunsli.cc)
target_link_libraries(dbrunsli PRIVATE
  brunslidec-static
)
if(BRUNSLI_EMSCRIPTEN)
set(WASM_MODULES brunslicodec-wasm brunslidec-wasm brunslienc-wasm)
foreach(module IN LISTS WASM_MODULES)
add_executable(${module} wasm/codec.cc)
target_link_libraries(${module} PRIVATE brunslidec-static brunslienc-static)
endforeach()
set(WASM_BASE_FLAGS "\
  -O3 \
  --closure 1 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -flto \
  --llvm-lto 1 \
  -s DISABLE_EXCEPTION_CATCHING=1 \
")
set(WASM_LINK_FLAGS "\
  ${WASM_BASE_FLAGS} \
  -s MODULARIZE=1 \
  -s FILESYSTEM=0 \
")
set_target_properties(cbrunsli PROPERTIES LINK_FLAGS "\
  ${WASM_BASE_FLAGS} \
  -s NODERAWFS=1 \
")
set_target_properties(dbrunsli PROPERTIES LINK_FLAGS "\
  ${WASM_BASE_FLAGS} \
  -s NODERAWFS=1 \
")
set(WASM_COMMON_EXPORT "\"_malloc\",\"_free\"")
set(WASM_DEC_EXPORT "\"_BrunsliToJpeg\",\"_GetJpegData\",\"_GetJpegLength\",\"_FreeJpeg\",\"_BrunsliDecoderInit\",\"_BrunsliDecoderProcess\",\"_BrunsliDecoderCleanup\"")
set(WASM_ENC_EXPORT "\"_JpegToBrunsli\",\"_GetBrunsliData\",\"_GetBrunsliLength\",\"_FreeBrunsli\"")
set_target_properties(brunslicodec-wasm PROPERTIES LINK_FLAGS "\
  ${WASM_LINK_FLAGS}\
  -s EXPORT_NAME=\"BrunsliCodecModule\"\
  -s EXPORTED_FUNCTIONS='[${WASM_COMMON_EXPORT},${WASM_DEC_EXPORT},${WASM_ENC_EXPORT}]'\
")
set_target_properties(brunslidec-wasm PROPERTIES LINK_FLAGS "\
  ${WASM_LINK_FLAGS}\
  -s EXPORT_NAME=\"BrunsliDecModule\"\
  -s EXPORTED_FUNCTIONS='[${WASM_COMMON_EXPORT},${WASM_DEC_EXPORT}]'\
")
set_target_properties(brunslienc-wasm PROPERTIES LINK_FLAGS "\
  ${WASM_LINK_FLAGS}\
  -s EXPORT_NAME=\"BrunsliEncModule\"\
  -s EXPORTED_FUNCTIONS='[${WASM_COMMON_EXPORT},${WASM_ENC_EXPORT}]'\
")
endif()  # BRUNSLI_EMSCRIPTEN

# Installation
if(NOT BRUNSLI_EMSCRIPTEN)
  install(
    TARGETS brunslidec-c brunslienc-c
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  )

  install(
    DIRECTORY ${BRUNSLI_INCLUDE_DIRS}/brunsli
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
endif() # BRUNSLI_EMSCRIPTEN

# Gather artifacts in a single directory for easier uploading.
set_target_properties(cbrunsli dbrunsli ${BRUNSLI_LIBRARIES} PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/artifacts"
  LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/artifacts"
  RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/artifacts"
)

if (${BUILD_TESTING})

include(GoogleTest)

set(BRUNSLI_TEST_ITEMS
    bit_reader
    build_huffman_table
    c_api
    context
    distributions
    fallback
    headerless
    huffman_tree
    lehmer_code
    quant_matrix
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests)
foreach (TEST_ITEM IN LISTS BRUNSLI_TEST_ITEMS)
  set(TEST_NAME ${TEST_ITEM}_test)
  add_executable(${TEST_NAME}
    c/tests/${TEST_NAME}.cc
    c/dec/decode.cc  # "static" brunslidec-c
    c/enc/encode.cc  # "static" brunslienc-c
    c/tests/test_utils.cc  # test utils
  )
  target_link_libraries(${TEST_NAME}
    brunslicommon-static
    brunslidec-static
    brunslienc-static
    gtest_main
  )
  gtest_discover_tests(${TEST_NAME})
endforeach()

endif()  # BUILD_TESTING
