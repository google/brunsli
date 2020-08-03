### Introduction
![CMake build](https://github.com/google/brunsli/workflows/CMake%20build/badge.svg?branch=master)


Brunsli is a lossless JPEG repacking library.

Brunsli allows for a 22% decrease in file size while allowing the original
JPEG to be recovered byte-by-byte.

It is possible to try how much Brunsli will save on your images on the site [brunsli.dev](https://brunsli.dev). Images are transcoded in browser, no data is transmitted or stored. Codec is powered by [WASM](https://webassembly.org/) technology. _Safari uses "interpretation" mode for WASM at first few runs after the page load. It is much slower. To avoid long wait, please feed codec with small images first._

<p align="center"><img alt="JPEG XL Logo" src="https://jpeg.org/images/jpegxl-logo.png" width="70px"></p>

VERY GOOD NEWS: Brunsli is on its way to become standardized. Brunsli has been specified as the lossless JPEG transport layer in the [Committee Draft of JPEG XL Image Coding System](https://arxiv.org/abs/1908.03565) and is ready to power faster and more economical transfer and storage of photographs.

We are committed to making JPEG XL a first-class citizen of the open-source and closed-source worlds, and we will integrate it into image and networking related tools.

The currently planned/on-going integration work includes:

- [x] one-shot C API / dynamic library
- [x] WASM module
- [ ] Node.js module
- [x] Java bindings
- [x] Python libraries support (OpenCV, imageio, PythonMagic, PIL, etc.)
- [ ] Python bindings
- [ ] Nginx transcoding module
- [x] Nginx serving module
- [ ] Apache transcoding module
- [x] Apache serving module

Stay tuned!

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
