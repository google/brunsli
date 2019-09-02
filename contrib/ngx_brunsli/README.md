# ngx_brunsli

[Brunsli][] is a lossless JPEG compressor. It is very fast, reduces JPEG file size
by approximately 22%, and is included in the [committee draft of the JPEG XL
standard][CD].

[Brunsli]: https://github.com/google/brunsli
[CD]: https://arxiv.org/abs/1908.03565

ngx_brunsli is an nginx module that serves pre-compressed Brunsli files.

## Status

Both Brunsli library and nginx module are under active development.

## Installation

### CentOS / RHEL 7

    yum install https://extras.getpagespeed.com/release-el7-latest.rpm
    yum install nginx nginx-module-nbr

Follow installation package prompts to enable the dynamic Brunsli module in `nginx.conf`:

    load_module modules/ngx_http_brunsli_static_module.so;

### Other Platforms - Dynamically loaded

    $ cd nginx-1.x.x
    $ ./configure --with-compat --add-dynamic-module=/path/to/ngx_brunsli
    $ make modules

You will need to use **exactly** the same `./configure` arguments as your Nginx configuration and append `--with-compat --add-dynamic-module=/path/to/ngx_brunsli` to the end, otherwise you will get a "module is not binary compatible" error on startup. You can run `nginx -V` to get the configuration arguments for your Nginx installation.

`make modules` will result in `ngx_http_brunsli_static_module.so` in the `objs` directory. Copy that module to `/usr/lib/nginx/modules/` then add the `load_module` lines above to `nginx.conf`.


### Other Platforms - Statically compiled

    $ cd nginx-1.x.x
    $ ./configure --add-module=/path/to/ngx_brunsli
    $ make && make install

This will compile the module directly into Nginx.

## Configuration directives

### `brunsli_static`

- **syntax**: `brunsli_static on|off|always`
- **default**: `off`
- **context**: `http`, `server`, `location`

Enables or disables checking of the existence of pre-compressed files with`.j`
extension. With the `always` value, pre-compressed file is used in all cases,
without checking if the client supports it.

## Contributing

See [Contributing](CONTRIBUTING.md).

## License

    Copyright (C) 2002-2015 Igor Sysoev
    Copyright (C) 2011-2015 Nginx, Inc.
    Copyright (C) 2015 Google Inc.
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:
    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in the
       documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
    OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
    OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
    SUCH DAMAGE.
