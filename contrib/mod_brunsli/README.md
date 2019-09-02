[Brunsli][] is a fast lossless JPEG recompressor that is included in the
[committee draft of the JPEG XL standard][CD]. This Apache module allows serving
precompressed Brunsli files.

[Brunsli]: https://github.com/google/brunsli
[CD]: https://arxiv.org/abs/1908.03565

# Installation

Compile the module with:

```console
$ apxs -c mod_brunsli.c
```

And install it with:

```console
$ apxs -i -a mod_brunsli.la
```

You can also do both steps at once with:

```console
$ apxs -i -a -c mod_brunsli.c
```

# Usage

Make sure that your Apache configuration contains a `LoadModule` directive for
`mod_brunsli` (`apxs -a` should have done that).

Then, for any file `foo.jpg`, if `foo.jpg.j` exists, then the latter will
automatically be served to clients that support Brunsli.
