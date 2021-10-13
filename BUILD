# Description:
#   Brunsli is an compact JPEG data format.

load(":compiler_config_setting.bzl", "create_msvc_config")

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

create_msvc_config()

STRICT_C_OPTIONS = select({
    ":msvc": [],
    "//conditions:default": [
        "-Wno-sign-compare",
    ],
})

# Real strict options:
# STRICT_C_OPTIONS = [
#    "--pedantic-errors",
#    "-Wall",
#    "-Wconversion",
#    "-Werror",
#    "-Wextra",
#    "-Wlong-long",
#    "-Wmissing-declarations",
#    "-Wno-strict-aliasing",
#    "-Wshadow",
#    "-Wsign-compare",
#]

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
        "@org_brotli//:brotlidec",
    ],
)

cc_library(
    name = "brunslienc",
    srcs = [":enc_sources"],
    hdrs = [":enc_headers"],
    copts = STRICT_C_OPTIONS,
    deps = [
        ":brunslicommon",
        "@org_brotli//:brotlienc",
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

config_setting(
    name = "experimental",
    define_values = {
        "brunsli_groups": "experimental",
    },
)

EXPERIMENTAL_DEPS = select({
    ":experimental": [":groups"],
    "//conditions:default": [],
})

EXPERIMENTAL_DEFINES = select({
    ":experimental": ["BRUNSLI_EXPERIMENTAL_GROUPS"],
    "//conditions:default": [],
})

EXPERIMENTAL_LINKOPTS = select({
    ":experimental": ["-pthread"],
    "//conditions:default": [],
})

cc_library(
    name = "groups",
    srcs = ["c/experimental/groups.cc"],
    hdrs = ["c/experimental/groups.h"],
    copts = STRICT_C_OPTIONS,
    deps = [
        ":brunslicommon",
        ":brunslidec",
        ":brunslienc",
    ],
)

cc_binary(
    name = "cbrunsli",
    srcs = ["c/tools/cbrunsli.cc"],
    copts = STRICT_C_OPTIONS,
    defines = EXPERIMENTAL_DEFINES,
    linkopts = EXPERIMENTAL_LINKOPTS,
    deps = [
        ":brunslicommon",
        ":brunslienc",
    ] + EXPERIMENTAL_DEPS,
)

cc_binary(
    name = "dbrunsli",
    srcs = ["c/tools/dbrunsli.cc"],
    copts = STRICT_C_OPTIONS,
    defines = EXPERIMENTAL_DEFINES,
    linkopts = EXPERIMENTAL_LINKOPTS,
    deps = [
        ":brunslicommon",
        ":brunslidec",
    ] + EXPERIMENTAL_DEPS,
)

TESTS = [
    "bit_reader",
    "build_huffman_table",
    "c_api",
    "context",
    "distributions",
    "fallback",
    "headerless",
    "huffman_tree",
    "lehmer_code",
    "quant_matrix",
    # "stream_decode", # fix brotli dependency
]

[cc_test(
    name = item + "_test",
    srcs = [
        "c/tests/" + item + "_test.cc",
        "c/tests/test_utils.cc",
        "c/tests/test_utils.h",
    ],
    copts = ["-Iexternal/gtest/include"],
    deps = [
        ":brunslicommon",
        ":brunslidec",
        ":brunslienc",
        ":brunslidec_c",
        ":brunslienc_c",
        "@com_google_googletest//:gtest_main",
    ],
) for item in TESTS]
