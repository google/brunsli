# Description:
#   Brunsli is an compact JPEG data format.

package(
    default_visibility = ["//visibility:public"],
)

licenses(["notice"])  # MIT

exports_files(["LICENSE"])

config_setting(
    name = "darwin",
    values = {"cpu": "darwin"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "darwin_x86_64",
    values = {"cpu": "darwin_x86_64"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "windows",
    values = {"cpu": "x64_windows"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "windows_msvc",
    values = {"cpu": "x64_windows_msvc"},
    visibility = ["//visibility:public"],
)

config_setting(
    name = "windows_msys",
    values = {"cpu": "x64_windows_msys"},
    visibility = ["//visibility:public"],
)

load(":compiler_config_setting.bzl", "create_msvc_config")

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
    name = "public_headers",
    srcs = glob(["c/include/brunsli/*.h"]),
)

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
    srcs = glob(
        ["c/dec/*.cc"],
        exclude = ["c/dec/decode.cc"],
    ),
)

filegroup(
    name = "enc_headers",
    srcs = glob(["c/enc/*.h"]),
)

filegroup(
    name = "enc_sources",
    srcs = glob(
        ["c/enc/*.cc"],
        exclude = ["c/enc/encode.cc"],
    ),
)

cc_library(
    name = "brunsli_inc",
    hdrs = [":public_headers"],
    copts = STRICT_C_OPTIONS,
    strip_include_prefix = "c/include",
)

cc_library(
    name = "brunslicommon",
    srcs = [":common_sources"],
    hdrs = [":common_headers"],
    copts = STRICT_C_OPTIONS,
    deps = [":brunsli_inc"],
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

cc_library(
    name = "brunslidec_c",
    srcs = ["c/dec/decode.cc"],
    copts = STRICT_C_OPTIONS,
    deps = [
        ":brunsli_inc",
        ":brunslidec",
    ],
)

cc_library(
    name = "brunslienc_c",
    srcs = ["c/enc/encode.cc"],
    copts = STRICT_C_OPTIONS,
    deps = [
        ":brunsli_inc",
        ":brunslienc",
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

cc_test(
    name = "bit_reader_test",
    srcs = ["c/tests/bit_reader_test.cc"],
    copts = ["-Iexternal/gtest/include"],
    deps = [
        ":brunslicommon",
        ":brunslidec",
        "@gtest//:gtest_main",
    ],
)

cc_test(
    name = "huffman_tree_test",
    srcs = ["c/tests/huffman_tree_test.cc"],
    copts = ["-Iexternal/gtest/include"],
    deps = [
        ":brunslidec",
        "@gtest//:gtest_main",
    ],
)

cc_test(
    name = "lehmer_code_test",
    srcs = ["c/tests/lehmer_code_test.cc"],
    copts = ["-Iexternal/gtest/include"],
    deps = [
        ":brunslicommon",
        "@gtest//:gtest_main",
    ],
)
