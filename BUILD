# Description:
#   Brunsli is an compact JPEG data format.

package(
    default_visibility = ["//visibility:public"],
)

licenses(["notice"])  # MIT

STRICT_C_OPTIONS = [
    "-Wno-sign-compare",  # everywhere
]

# Real strict options:
STRICT_C_OPTIONS_ = [
    "--pedantic-errors",
    "-Wall",
    "-Wconversion",
    "-Werror",
    "-Wextra",
    "-Wlong-long",
    "-Wmissing-declarations",
    "-Wno-strict-aliasing",
    "-Wshadow",
    "-Wsign-compare",
]

filegroup(
    name = "common_headers",
    srcs = glob(["c/common/*.h"]),
)

filegroup(
    name = "common_sources",
    srcs = glob(["c/common/*.cc"]),
)

filegroup(
    name = "dec_headers",
    srcs = glob(["c/dec/*.h"]),
)

filegroup(
    name = "dec_sources",
    srcs = glob(["c/dec/*.cc"]),
)

filegroup(
    name = "enc_headers",
    srcs = glob(["c/enc/*.h"]),
)

filegroup(
    name = "enc_sources",
    srcs = glob(["c/enc/*.cc"]),
)

cc_library(
    name = "brunslicommon",
    srcs = [":common_sources"],
    hdrs = [":common_headers"],
    copts = STRICT_C_OPTIONS,
)

cc_library(
    name = "brunslidec",
    srcs = [":dec_sources"],
    hdrs = [":dec_headers"],
    copts = STRICT_C_OPTIONS,
    deps = [
        ":brunslicommon",
        "@brotli//:brotlidec",
    ],
)

cc_library(
    name = "brunslienc",
    srcs = [":enc_sources"],
    hdrs = [":enc_headers"],
    copts = STRICT_C_OPTIONS,
    deps = [
        ":brunslicommon",
        "@brotli//:brotlienc",
    ],
)

cc_binary(
    name = "cbrunsli",
    srcs = ["c/tools/cbrunsli.cc"],
    copts = STRICT_C_OPTIONS,
    deps = [
        ":brunslicommon",
        ":brunslienc",
    ],
)

cc_binary(
    name = "dbrunsli",
    srcs = ["c/tools/dbrunsli.cc"],
    copts = STRICT_C_OPTIONS,
    deps = [
        ":brunslicommon",
        ":brunslidec",
    ],
)
