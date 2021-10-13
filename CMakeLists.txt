cmake_minimum_required(VERSION 3.11)

project(BRUNSLI C CXX)

include(CheckCSourceCompiles)
check_c_source_compiles(
  "#if defined(__EMSCRIPTEN__)
   int main() {return 0;}
   #endif"
  BRUNSLI_EMSCRIPTEN
)
if (BRUNSLI_EMSCRIPTEN)
  message("-- Compiler is EMSCRIPTEN")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -flto -fno-exceptions -fno-rtti")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -flto -fno-exceptions -fno-rtti")
  set(CMAKE_STATIC_LIBRARY_SUFFIX ".bc")
  set(CMAKE_EXECUTABLE_SUFFIX ".js")
else()
  message("-- Compiler is not EMSCRIPTEN")
endif()

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

include(FetchContent)

# Add GTest
FetchContent_Declare(
  GTest
  GIT_REPOSITORY https://github.com/google/googletest
  GIT_TAG e2239ee6043f73722e7aa812a459f54a28552929 # v1.11.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(GTest)

# Add Brotli
FetchContent_Declare(
  Brotli
  GIT_REPOSITORY https://github.com/google/brotli
  GIT_TAG e61745a6b7add50d380cfd7d3883dd6c62fc2c71 # v1.0.9
)
set(BROTLI_DISABLE_TESTS ON CACHE INTERNAL "")
FetchContent_MakeAvailable(Brotli)

include(CTest)

# The Brunsli library definition.
include(brunsli.cmake)
