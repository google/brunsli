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

file(GLOB BRUNSLI_DEC_SOURCES
  c/dec/*.cc
)

# TODO(eustas): split public/private headers.
file(GLOB BRUNSLI_DEC_HEADERS
  c/dec/*.h
)

file(GLOB BRUNSLI_ENC_SOURCES
  c/enc/*.cc
)

# TODO(eustas): split public/private headers.
file(GLOB BRUNSLI_ENC_HEADERS
  c/enc/*.h
)

add_library(brunslicommon-static STATIC
  ${BRUNSLI_COMMON_SOURCES}
  ${BRUNSLI_COMMON_HEADERS}
)

add_library(brunslidec-static STATIC
  ${BRUNSLI_DEC_SOURCES}
  ${BRUNSLI_DEC_HEADERS}
)

add_library(brunslienc-static STATIC
  ${BRUNSLI_ENC_SOURCES}
  ${BRUNSLI_ENC_HEADERS}
)

target_link_libraries(brunslidec-static PRIVATE
  brotlidec-static
  brunslicommon-static
)

target_link_libraries(brunslienc-static PRIVATE
  brotlienc-static
  brunslicommon-static
)

foreach(lib brunslicommon-static brunslidec-static brunslienc-static)
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang" OR
     "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    target_compile_options(${lib} PUBLIC
      # Debug flags
      -dwarf-column-info
      -debug-info-kind=line-tables-only
      -dwarf-version=4
      -debugger-tuning=gdb

      # F_FLAGS
      -fmerge-all-constants
      -fno-builtin-fwrite
      -fno-builtin-fread
      -fno-signed-char
      -fsized-deallocation
      -fnew-alignment=8
      -fno-cxx-exceptions
      -fno-exceptions
      -fno-slp-vectorize
      -fno-vectorize

      # WARN_FLAGS
      -Wformat-security
      -Wno-char-subscripts
      -Wno-error=deprecated-declarations
      -Wno-sign-compare
      -Wno-strict-overflow
      -Wno-unused-function
      -Wthread-safety-analysis
      -Wno-unknown-warning-option
      -Wno-unused-command-line-argument
      -Wno-ignored-optimization-argument
      -Wno-ambiguous-member-template
      -Wno-pointer-sign
      -Wno-address-of-packed-member
      -Wno-enum-compare-switch
      -Wno-expansion-to-defined
      -Wno-extern-c-compat
      -Wno-gnu-alignof-expression
      -Wno-gnu-designator
      -Wno-gnu-variable-sized-type-not-at-end
      -Wno-ignored-attributes
      -Wno-ignored-qualifiers
      -Wno-inconsistent-missing-override
      -Wno-invalid-source-encoding
      -Wno-mismatched-tags
      -Wno-potentially-evaluated-expression
      -Wno-return-std-move
      -Wno-self-assign-overloaded
      -Wno-tautological-constant-compare
      -Wno-tautological-constant-in-range-compare
      -Wno-tautological-type-limit-compare
      -Wno-tautological-undefined-compare
      -Wno-tautological-unsigned-zero-compare
      -Wno-tautological-unsigned-enum-zero-compare
      -Wno-undefined-func-template
      -Wno-unknown-pragmas
      -Wno-unused-const-variable
      -Wno-unused-lambda-capture
      -Wno-unused-local-typedef
      -Wno-unused-private-field
      -Wno-private-header
      -Wfloat-overflow-conversion
      -Wfloat-zero-conversion
      -Wfor-loop-analysis
      -Wgnu-redeclared-enum
      -Winfinite-recursion
      -Wliteral-conversion
      -Wself-assign
      -Wstring-conversion
      -Wtautological-overlap-compare
      -Wunused-comparison
      -Wvla
      -Wno-reserved-user-defined-literal
      -Wno-return-type-c-linkage
      -Wno-deprecated
      -Wno-invalid-offsetof
      -Wno-literal-suffix
      -Woverloaded-virtual
      -Wnon-virtual-dtor
      -Wdeprecated-increment-bool
      -Wc++11-compat
      -Wno-c++11-compat-binary-literal
      -Wc++2a-extensions
      -Wno-register
      -Wno-dynamic-exception-spec
      -Wprivate-header
      -Wno-builtin-macro-redefined

      # Language flags
      -disable-free
      -disable-llvm-verifier
      -discard-value-names
      # Note: this works only because this is the only -Xclang passed.
      -Xclang -relaxed-aliasing
      -fmath-errno
    )
  endif()

  target_include_directories(${lib}
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
endforeach()

add_executable(cbrunsli c/tools/cbrunsli.cc)
target_link_libraries(cbrunsli PRIVATE
  brunslienc-static
)

add_executable(dbrunsli c/tools/dbrunsli.cc)
target_link_libraries(dbrunsli PRIVATE
  brunslidec-static
)
