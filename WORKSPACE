workspace(name = "dev_brunsli")

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

local_repository(
    name = "ignore_dev_brunsli_java",
    path = "java",
)

http_archive(
    name = "com_google_googletest",
    urls = ["https://github.com/google/googletest/archive/e2239ee6043f73722e7aa812a459f54a28552929.zip"],
    sha256 = "8daa1a71395892f7c1ec5f7cb5b099a02e606be720d62f1a6a98f8f8898ec826",
    strip_prefix = "googletest-e2239ee6043f73722e7aa812a459f54a28552929",
)

http_archive(
    name = "org_brotli",
    urls = ["https://github.com/google/brotli/archive/ed738e842d2fbdf2d6459e39267a633c4a9b2f5d.zip"],
    sha256 = "a68ec12a898abc9cf248f21362620562041b7aab4d623ecd736f39bedf5002a0",
    strip_prefix = "brotli-ed738e842d2fbdf2d6459e39267a633c4a9b2f5d",
)
