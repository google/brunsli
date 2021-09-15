### Introduction
![CMake build](https://github.com/google/brunsli/workflows/Build-Test-Upload/badge.svg?branch=master)


Brunsli is a lossless JPEG repacking library.

Brunsli allows for a 22% decrease in file size while allowing the original
JPEG to be recovered byte-by-byte.

It is possible to try how much Brunsli will save on your images on the site [brunsli.dev](https://brunsli.dev). Images are transcoded in browser, no data is transmitted or stored. Codec is powered by [WASM](https://webassembly.org/) technology. _Safari uses "interpretation" mode for WASM at first few runs after the page load. It is much slower. To avoid long wait, please feed codec with small images first._

## Build instructions

Run the following commands to clone, configure and build Brunsli:

    git clone --depth=1 https://github.com/google/brunsli.git
    cd brunsli
    git submodule update --init --recursive
    cmake -DCMAKE_BUILD_TYPE=Release -B out
    cmake --build out --config Release

## Prebuilt binaries

For some platforms (e.g. Windows) libraries and executables are uploaded as "artifacts" as a part of continous integration process.
Unfortunately, there is no static link to access those. Please follow the [GitHub manual](https://docs.github.com/en/actions/configuring-and-managing-workflows/persisting-workflow-data-using-artifacts#downloading-and-deleting-artifacts-after-a-workflow-run-is-complete).
