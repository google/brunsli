cmake_minimum_required(VERSION 3.1)

project(BRUNSLI C CXX)

# Recommend clang for building.
if(NOT "${CMAKE_C_COMPILER_ID}" MATCHES "Clang" OR
   NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  message(WARNING "Using ${CMAKE_CXX_COMPILER_ID} compiler.\n"
    "For best results, use clang instead:\n  CC=clang CXX=clang++ cmake ..")
  if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    if ("${CMAKE_C_COMPILER_VERSION}" VERSION_LESS 6)
      message(FATAL_ERROR
        "Minimum Clang version required is Clang 6, please update.")
    endif()
  endif()
  if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS 6)
      message(FATAL_ERROR
        "Minimum Clang version required is Clang 6, please update.")
    endif()
  endif()
endif()

# CMAKE_EXPORT_COMPILE_COMMANDS is used to generate the compilation database
# used by clang-tidy.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang" OR
   "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

  # Global compiler flags for all targets here and in subdirectories.
  add_definitions(
    # Define "register" as empty since it was deprecated already.
    -Dregister=

    # Avoid changing the binary based on the current time and date.
    -D__DATE__="redacted"
    -D__TIMESTAMP__="redacted"
    -D__TIME__="redacted"
  )

  # In CMake before 3.12 it is problematic to pass repeated flags like -Xclang.
  # For this reason we place them in CMAKE_CXX_FLAGS instead.
  # See https://gitlab.kitware.com/cmake/cmake/issues/15826

  # Pretty colorful messages within reasonable limits.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -Xclang -ferror-limit -Xclang 19 -Xclang -fmessage-length -Xclang 0 \
    -fdiagnostics-show-option -fcolor-diagnostics")

  # Machine flags.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -mrelax-all \
    -Xclang -mrelocation-model -Xclang pic \
    -Xclang -pic-level -Xclang 2 \
    -Xclang -mdisable-fp-elim \
    -Xclang -mconstructor-aliases \
    -mpie-copy-relocations \
    -Xclang -munwind-tables")

  # CPU flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -mavx2 \
    -mfma \
    -Xclang -mprefer-vector-width=128 \
    -Xclang -target-cpu -Xclang haswell \
    -Xclang -target-feature -Xclang +avx2")

  add_compile_options(
    # Ignore this to allow redefining __DATE__ and others.
    -Wno-builtin-macro-redefined
    
    # Global warning settings.
    -Wall
    -Werror
  )  
endif()

# Force build with optimizations in release mode.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_subdirectory(third_party)

# The Brunsli library definition.
include(brunsli.cmake)
