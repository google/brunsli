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
    urls = ["https://github.com/google/brotli/archive/dbfebd13dcf51e1360f4f943306282b396a3e8cd.zip"],
    sha256 = "a87127241f52cf7a78f8e56f6c46ebb37e5507fc5dc94e60abf40933c8cec875",
    strip_prefix = "brotli-dbfebd13dcf51e1360f4f943306282b396a3e8cd",
)
