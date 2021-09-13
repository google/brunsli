workspace(name = "dev_brunsli")

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

local_repository(
    name = "ignore_dev_brunsli_java",
    path = "java",
)

local_repository(
    name = "brotli",
    path = "third_party/brotli",
)

local_repository(
    name = "ignore_googletest",
    path = "third_party/googletest",
)

http_archive(
    name = "com_google_googletest",
    urls = ["https://github.com/google/googletest/archive/e2239ee6043f73722e7aa812a459f54a28552929.zip"],
    sha256 = "8daa1a71395892f7c1ec5f7cb5b099a02e606be720d62f1a6a98f8f8898ec826",
    strip_prefix = "googletest-e2239ee6043f73722e7aa812a459f54a28552929",
)

new_local_repository(
    name = "highwayhash",
    path = "third_party/highwayhash",
    build_file_content = """
package(default_visibility = ["//visibility:public"])
cc_library(
    name = "highwayhash_inc",
    hdrs = glob(["highwayhash/*.h"]),
    strip_include_prefix = "",
)
    """,
)
