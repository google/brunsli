package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # MIT

filegroup(
    name = "jni_src",
    srcs = [
        "codec_jni.cc",
        "codec_jni.h",
        "codec_jni_onload.cc",
    ],
)

filegroup(
    name = "brunsli_jni",
    srcs = ["//:brunsli_jni.dll"],
)

java_library(
    name = "codec",
    srcs = glob(
        ["*.java"],
        exclude = ["*Test*.java"],
    ),
    resource_jars = ["//:license"],
)

java_library(
    name = "test_lib",
    testonly = 1,
    srcs = glob(["*Test*.java"]),
    deps = [
        ":codec",
        "@maven//:junit_junit",
    ],
)

java_test(
    name = "CodecTest",
    size = "small",
    data = [
        ":brunsli_jni",  # Bazel JNI workaround
    ],
    jvm_flags = [
        "-DBRUNSLI_JNI_LIBRARY=$(location :brunsli_jni)",
    ],
    test_class = "dev.brunsli.wrapper.CodecTest",
    runtime_deps = [":test_lib"],
)
