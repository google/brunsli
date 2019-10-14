### Introduction
[![Travis Build Status](https://travis-ci.org/google/brunsli.svg?branch=master)](https://travis-ci.org/google/brunsli) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/google/brunsli?branch=master&svg=true)](https://ci.appveyor.com/project/google/brunsli)

Brunsli is a lossless JPEG repacking library.

Brunsli allows for a 22% decrease in file size while allowing the original
JPEG to be recovered byte-by-byte.

<p align="center"><img alt="JPEG XL Logo" src="https://jpeg.org/images/jpegxl-logo.png" width="70px"></p>

VERY GOOD NEWS: Brunsli is on its way to become standardized. Brunsli has been specified as the lossless JPEG transport layer in the [Committee Draft of JPEG XL Image Coding System](https://arxiv.org/abs/1908.03565) and is ready to power faster and more economical transfer and storage of photographs.

We are committed making JPEG XL a first-class citizen of the open-source and closed-source worlds, and we will integrate it into image and networking related tools.

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

### Build instructions

The basic commands to clone, build, and install brunsli are:

    $ git clone https://github.com/google/brunsli.git
    $ cd brunsli
    $ git submodule update --init --recursive
    $ cmake ./
    $ make -j
    $ make -j install
