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

# Force build with optimizations in release mode.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_subdirectory(third_party)

# The Brunsli library definition.
include(brunsli.cmake)
