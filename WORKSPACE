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
    urls = ["https://github.com/google/brotli/archive/e61745a6b7add50d380cfd7d3883dd6c62fc2c71.zip"],
    sha256 = "4a79fd9fd30bae4d08dab373326cfb21ab0d6b50e0e55564043e35dde7210219",
    strip_prefix = "brotli-e61745a6b7add50d380cfd7d3883dd6c62fc2c71",
)
