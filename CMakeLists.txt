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
  message(STATUS "Switching to static build for EMSCRIPTEN")
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -flto -fno-exceptions -fno-rtti")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -flto -fno-exceptions -fno-rtti")
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
  GIT_TAG 52eb8108c5bdec04579160ae17225d66034bd723 # v1.17.0
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(GTest)

# Add Brotli
FetchContent_Declare(
  Brotli
  GIT_REPOSITORY https://github.com/google/brotli
  GIT_TAG 028fb5a23661f123017c060daa546b55cf4bde29 # v1.2.0
)
if (NOT BRUNSLI_EMSCRIPTEN)
  # We have already set "BUILD_SHARED_LIBS"
  set(BROTLI_BUILD_FOR_PACKAGE ON CACHE INTERNAL "")
endif()
set(BROTLI_DISABLE_TESTS ON CACHE INTERNAL "")
set(BROTLI_BUILD_TOOLS OFF CACHE INTERNAL "")
FetchContent_MakeAvailable(Brotli)

include(CTest)

# The Brunsli library definition.
include(brunsli.cmake)
