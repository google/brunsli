workspace(name = "dev_brunsli")

local_repository(
    name = "ignore_dev_brunsli_java",
    path = "java",
)

local_repository(
    name = "brotli",
    path = "third_party/brotli",
)

new_local_repository(
    name = "gtest",
    path = "third_party/googletest",
    build_file = "third_party/googletest/BUILD.bazel",
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
