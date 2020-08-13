# jxl Content-Encoding explained

## What’s all this then?

HTTP Content-Encoding allows transfer of resources in a compressed form.
Most popular Content-Encodings – `deflate` (`gzip`) and `br` – are LZ77-based general-purpose data compressors; those fit well for compression of HTML, JS and CSS.

Video and image resources take 1-st and 2-nd place in the Internet traffic volume chart. Those are considered succint.
However, one of the most popular image formats – JPEG – could be "repacked" to a more dense form.

`jxl` Content-Encoding allows ~22% traffic reduction for JPEG images.

## Getting started

First, install appropriate server plugin (e.g. `ngx_brunsli` for NGINX or `mod_brunsli` for Apache from [https://github.com/google/brunsli/tree/master/contrib](contrib)).

Now server is capable to respond with compressed resources when it sees appropriate encoding in 'Accept-Encoding' header.

Restriction: it is likely, that browser will send `jxl` in 'Accept-Encoding' header only over encrypted connections (HTTPS) to avoid problems with faulty proxies.

## Links

Codec demo site: https://brunsli.dev/

Specification: https://arxiv.org/pdf/1908.03565.pdf
