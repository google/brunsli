cmake_minimum_required(VERSION 3.1)

project(BRUNSLI C CXX)

# Force build with optimizations in release mode.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

include(CheckCSourceCompiles)
check_c_source_compiles(
  "#if defined(__EMSCRIPTEN__)
   int main() {return 0;}
   #endif"
  BRUNSLI_EMSCRIPTEN
)
if (BRUNSLI_EMSCRIPTEN)
  message("-- Compiler is EMSCRIPTEN")
else()
  message("-- Compiler is not EMSCRIPTEN")
endif()

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_subdirectory(third_party)

# The Brunsli library definition.
include(brunsli.cmake)
