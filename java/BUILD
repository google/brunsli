package(
    default_visibility = ["//visibility:public"],
)

# >>> JNI headers

genrule(
    name = "copy_link_jni_header",
    srcs = ["@openjdk_jni_h//file"],
    outs = ["jni/jni.h"],
    cmd = "cp -f $< $@",
)

genrule(
    name = "copy_link_jni_md_header",
    srcs = select({
        "@dev_brunsli//:darwin": ["@openjdk_macosx_jni_md_h//file"],
        "@dev_brunsli//:darwin_x86_64": ["@openjdk_macosx_jni_md_h//file"],
        "@dev_brunsli//:windows_msys": ["@openjdk_windows_jni_md_h//file"],
        "@dev_brunsli//:windows_msvc": ["@openjdk_windows_jni_md_h//file"],
        "@dev_brunsli//:windows": ["@openjdk_windows_jni_md_h//file"],
        "//conditions:default": ["@openjdk_solaris_jni_md_h//file"],
    }),
    outs = ["jni/jni_md.h"],
    cmd = "cp -f $< $@",
)

cc_library(
    name = "jni_inc",
    hdrs = [
        ":jni/jni.h",
        ":jni/jni_md.h",
    ],
    includes = ["jni"],
)

# <<< JNI headers

genrule(
    name = "license_resource",
    srcs = ["@dev_brunsli//:LICENSE"],
    outs = ["META-INF/LICENSE"],
    cmd = "cp -f $< $@",
)

java_library(
    name = "license",
    resources = [":license_resource"],
)

########################################################
# WARNING: do not (transitively) depend on this target!
########################################################
cc_binary(
    name = "brunsli_jni.dll",
    srcs = [
        "//dev/brunsli/wrapper:jni_src",
    ],
    linkopts = [
        "-Wl,--version-script,$(location :jni.ldscript)",
    ],
    linkshared = 1,
    deps = [
        ":jni.ldscript",
        ":jni_inc",
        "@dev_brunsli//:brunsli_inc",
        "@dev_brunsli//:brunslidec",
        "@dev_brunsli//:brunslienc",
    ],
)
