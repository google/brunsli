// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

/* Macros for compiler / platform specific API declarations. */

#ifndef BRUNSLI_COMMON_PORT_H_
#define BRUNSLI_COMMON_PORT_H_

/* The following macros were borrowed from https://github.com/nemequ/hedley
 * with permission of original author - Evan Nemerson <evan@nemerson.com> */

/* >>> >>> >>> hedley macros */

#define BRUNSLI_MAKE_VERSION(major, minor, revision) \
  (((major) * 1000000) + ((minor) * 1000) + (revision))

#if defined(__GNUC__) && defined(__GNUC_PATCHLEVEL__)
#define BRUNSLI_GNUC_VERSION \
  BRUNSLI_MAKE_VERSION(__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__)
#elif defined(__GNUC__)
#define BRUNSLI_GNUC_VERSION BRUNSLI_MAKE_VERSION(__GNUC__, __GNUC_MINOR__, 0)
#endif

#if defined(BRUNSLI_GNUC_VERSION)
#define BRUNSLI_GNUC_VERSION_CHECK(major, minor, patch) \
  (BRUNSLI_GNUC_VERSION >= BRUNSLI_MAKE_VERSION(major, minor, patch))
#else
#define BRUNSLI_GNUC_VERSION_CHECK(major, minor, patch) (0)
#endif

#if defined(_MSC_FULL_VER) && (_MSC_FULL_VER >= 140000000)
#define BRUNSLI_MSVC_VERSION                               \
  BRUNSLI_MAKE_VERSION((_MSC_FULL_VER / 10000000),         \
                      (_MSC_FULL_VER % 10000000) / 100000, \
                      (_MSC_FULL_VER % 100000) / 100)
#elif defined(_MSC_FULL_VER)
#define BRUNSLI_MSVC_VERSION                             \
  BRUNSLI_MAKE_VERSION((_MSC_FULL_VER / 1000000),        \
                      (_MSC_FULL_VER % 1000000) / 10000, \
                      (_MSC_FULL_VER % 10000) / 10)
#elif defined(_MSC_VER)
#define BRUNSLI_MSVC_VERSION \
  BRUNSLI_MAKE_VERSION(_MSC_VER / 100, _MSC_VER % 100, 0)
#endif

#if !defined(_MSC_VER)
#define BRUNSLI_MSVC_VERSION_CHECK(major, minor, patch) (0)
#elif defined(_MSC_VER) && (_MSC_VER >= 1400)
#define BRUNSLI_MSVC_VERSION_CHECK(major, minor, patch) \
  (_MSC_FULL_VER >= ((major * 10000000) + (minor * 100000) + (patch)))
#elif defined(_MSC_VER) && (_MSC_VER >= 1200)
#define BRUNSLI_MSVC_VERSION_CHECK(major, minor, patch) \
  (_MSC_FULL_VER >= ((major * 1000000) + (minor * 10000) + (patch)))
#else
#define BRUNSLI_MSVC_VERSION_CHECK(major, minor, patch) \
  (_MSC_VER >= ((major * 100) + (minor)))
#endif

#if defined(__INTEL_COMPILER) && defined(__INTEL_COMPILER_UPDATE)
#define BRUNSLI_INTEL_VERSION                  \
  BRUNSLI_MAKE_VERSION(__INTEL_COMPILER / 100, \
                      __INTEL_COMPILER % 100,  \
                      __INTEL_COMPILER_UPDATE)
#elif defined(__INTEL_COMPILER)
#define BRUNSLI_INTEL_VERSION \
  BRUNSLI_MAKE_VERSION(__INTEL_COMPILER / 100, __INTEL_COMPILER % 100, 0)
#endif

#if defined(BRUNSLI_INTEL_VERSION)
#define BRUNSLI_INTEL_VERSION_CHECK(major, minor, patch) \
  (BRUNSLI_INTEL_VERSION >= BRUNSLI_MAKE_VERSION(major, minor, patch))
#else
#define BRUNSLI_INTEL_VERSION_CHECK(major, minor, patch) (0)
#endif

#if defined(__PGI) && \
    defined(__PGIC__) && defined(__PGIC_MINOR__) && defined(__PGIC_PATCHLEVEL__)
#define BRUNSLI_PGI_VERSION \
  BRUNSLI_MAKE_VERSION(__PGIC__, __PGIC_MINOR__, __PGIC_PATCHLEVEL__)
#endif

#if defined(BRUNSLI_PGI_VERSION)
#define BRUNSLI_PGI_VERSION_CHECK(major, minor, patch) \
  (BRUNSLI_PGI_VERSION >= BRUNSLI_MAKE_VERSION(major, minor, patch))
#else
#define BRUNSLI_PGI_VERSION_CHECK(major, minor, patch) (0)
#endif

#if defined(__SUNPRO_C) && (__SUNPRO_C > 0x1000)
#define BRUNSLI_SUNPRO_VERSION                                      \
  BRUNSLI_MAKE_VERSION(                                             \
    (((__SUNPRO_C >> 16) & 0xf) * 10) + ((__SUNPRO_C >> 12) & 0xf), \
    (((__SUNPRO_C >> 8) & 0xf) * 10) + ((__SUNPRO_C >> 4) & 0xf),   \
    (__SUNPRO_C & 0xf) * 10)
#elif defined(__SUNPRO_C)
#define BRUNSLI_SUNPRO_VERSION                  \
  BRUNSLI_MAKE_VERSION((__SUNPRO_C >> 8) & 0xf, \
                      (__SUNPRO_C >> 4) & 0xf, \
                      (__SUNPRO_C) & 0xf)
#elif defined(__SUNPRO_CC) && (__SUNPRO_CC > 0x1000)
#define BRUNSLI_SUNPRO_VERSION                                        \
  BRUNSLI_MAKE_VERSION(                                               \
    (((__SUNPRO_CC >> 16) & 0xf) * 10) + ((__SUNPRO_CC >> 12) & 0xf), \
    (((__SUNPRO_CC >> 8) & 0xf) * 10) + ((__SUNPRO_CC >> 4) & 0xf),   \
    (__SUNPRO_CC & 0xf) * 10)
#elif defined(__SUNPRO_CC)
#define BRUNSLI_SUNPRO_VERSION                   \
  BRUNSLI_MAKE_VERSION((__SUNPRO_CC >> 8) & 0xf, \
                      (__SUNPRO_CC >> 4) & 0xf, \
                      (__SUNPRO_CC) & 0xf)
#endif

#if defined(BRUNSLI_SUNPRO_VERSION)
#define BRUNSLI_SUNPRO_VERSION_CHECK(major, minor, patch) \
  (BRUNSLI_SUNPRO_VERSION >= BRUNSLI_MAKE_VERSION(major, minor, patch))
#else
#define BRUNSLI_SUNPRO_VERSION_CHECK(major, minor, patch) (0)
#endif

#if defined(__CC_ARM) && defined(__ARMCOMPILER_VERSION)
#define BRUNSLI_ARM_VERSION                                      \
  BRUNSLI_MAKE_VERSION((__ARMCOMPILER_VERSION / 1000000),        \
                      (__ARMCOMPILER_VERSION % 1000000) / 10000, \
                      (__ARMCOMPILER_VERSION % 10000) / 100)
#elif defined(__CC_ARM) && defined(__ARMCC_VERSION)
#define BRUNSLI_ARM_VERSION                                \
  BRUNSLI_MAKE_VERSION((__ARMCC_VERSION / 1000000),        \
                      (__ARMCC_VERSION % 1000000) / 10000, \
                      (__ARMCC_VERSION % 10000) / 100)
#endif

#if defined(BRUNSLI_ARM_VERSION)
#define BRUNSLI_ARM_VERSION_CHECK(major, minor, patch) \
  (BRUNSLI_ARM_VERSION >= BRUNSLI_MAKE_VERSION(major, minor, patch))
#else
#define BRUNSLI_ARM_VERSION_CHECK(major, minor, patch) (0)
#endif

#if defined(__ibmxl__)
#define BRUNSLI_IBM_VERSION                   \
  BRUNSLI_MAKE_VERSION(__ibmxl_version__,     \
                      __ibmxl_release__,      \
                      __ibmxl_modification__)
#elif defined(__xlC__) && defined(__xlC_ver__)
#define BRUNSLI_IBM_VERSION \
  BRUNSLI_MAKE_VERSION(__xlC__ >> 8, __xlC__ & 0xff, (__xlC_ver__ >> 8) & 0xff)
#elif defined(__xlC__)
#define BRUNSLI_IBM_VERSION \
  BRUNSLI_MAKE_VERSION(__xlC__ >> 8, __xlC__ & 0xff, 0)
#endif

#if defined(BRUNSLI_IBM_VERSION)
#define BRUNSLI_IBM_VERSION_CHECK(major, minor, patch) \
  (BRUNSLI_IBM_VERSION >= BRUNSLI_MAKE_VERSION(major, minor, patch))
#else
#define BRUNSLI_IBM_VERSION_CHECK(major, minor, patch) (0)
#endif

#if defined(__TI_COMPILER_VERSION__)
#define BRUNSLI_TI_VERSION                                        \
  BRUNSLI_MAKE_VERSION((__TI_COMPILER_VERSION__ / 1000000),       \
                      (__TI_COMPILER_VERSION__ % 1000000) / 1000, \
                      (__TI_COMPILER_VERSION__ % 1000))
#endif

#if defined(BRUNSLI_TI_VERSION)
#define BRUNSLI_TI_VERSION_CHECK(major, minor, patch) \
  (BRUNSLI_TI_VERSION >= BRUNSLI_MAKE_VERSION(major, minor, patch))
#else
#define BRUNSLI_TI_VERSION_CHECK(major, minor, patch) (0)
#endif

#if defined(__IAR_SYSTEMS_ICC__)
#if __VER__ > 1000
#define BRUNSLI_IAR_VERSION                    \
  BRUNSLI_MAKE_VERSION((__VER__ / 1000000),    \
                      (__VER__ / 1000) % 1000, \
                      (__VER__ % 1000))
#else
#define BRUNSLI_IAR_VERSION BRUNSLI_MAKE_VERSION(VER / 100, __VER__ % 100, 0)
#endif
#endif

#if defined(BRUNSLI_IAR_VERSION)
#define BRUNSLI_IAR_VERSION_CHECK(major, minor, patch) \
  (BRUNSLI_IAR_VERSION >= BRUNSLI_MAKE_VERSION(major, minor, patch))
#else
#define BRUNSLI_IAR_VERSION_CHECK(major, minor, patch) (0)
#endif

#if defined(__TINYC__)
#define BRUNSLI_TINYC_VERSION                                    \
  BRUNSLI_MAKE_VERSION(__TINYC__ / 1000, (__TINYC__ / 100) % 10, \
                       __TINYC__ % 100)
#endif

#if defined(BRUNSLI_TINYC_VERSION)
#define BRUNSLI_TINYC_VERSION_CHECK(major, minor, patch) \
  (BRUNSLI_TINYC_VERSION >= BRUNSLI_MAKE_VERSION(major, minor, patch))
#else
#define BRUNSLI_TINYC_VERSION_CHECK(major, minor, patch) (0)
#endif

#if defined(__has_attribute)
#define BRUNSLI_GNUC_HAS_ATTRIBUTE(attribute, major, minor, patch) \
  __has_attribute(attribute)
#else
#define BRUNSLI_GNUC_HAS_ATTRIBUTE(attribute, major, minor, patch) \
  BRUNSLI_GNUC_VERSION_CHECK(major, minor, patch)
#endif

#if defined(__has_builtin)
#define BRUNSLI_GNUC_HAS_BUILTIN(builtin, major, minor, patch) \
  __has_builtin(builtin)
#else
#define BRUNSLI_GNUC_HAS_BUILTIN(builtin, major, minor, patch) \
  BRUNSLI_GNUC_VERSION_CHECK(major, minor, patch)
#endif

#if defined(_WIN32) || defined(__CYGWIN__)
#define BRUNSLI_PUBLIC
#elif BRUNSLI_GNUC_VERSION_CHECK(3, 3, 0) ||                        \
    BRUNSLI_TI_VERSION_CHECK(8, 0, 0) ||                            \
    BRUNSLI_INTEL_VERSION_CHECK(16, 0, 0) ||                        \
    BRUNSLI_ARM_VERSION_CHECK(4, 1, 0) ||                           \
    BRUNSLI_IBM_VERSION_CHECK(13, 1, 0) ||                          \
    BRUNSLI_SUNPRO_VERSION_CHECK(5, 11, 0) ||                       \
    (BRUNSLI_TI_VERSION_CHECK(7, 3, 0) &&                           \
     defined(__TI_GNU_ATTRIBUTE_SUPPORT__) && defined(__TI_EABI__))
#define BRUNSLI_PUBLIC __attribute__ ((visibility ("default")))
#else
#define BRUNSLI_PUBLIC
#endif

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L) && \
    !defined(__STDC_NO_VLA__) && !defined(__cplusplus) &&         \
    !defined(__PGI) && !defined(__PGIC__) && !defined(__TINYC__)
#define BRUNSLI_ARRAY_PARAM(name) (name)
#else
#define BRUNSLI_ARRAY_PARAM(name)
#endif

/* <<< <<< <<< end of hedley macros. */

#if defined(BRUNSLI_SHARED_COMPILATION)
#if defined(_WIN32)
#if defined(BRUNSLICOMMON_SHARED_COMPILATION)
#define BRUNSLI_COMMON_API __declspec(dllexport)
#else
#define BRUNSLI_COMMON_API __declspec(dllimport)
#endif  /* BRUNSLICOMMON_SHARED_COMPILATION */
#if defined(BRUNSLIDEC_SHARED_COMPILATION)
#define BRUNSLI_DEC_API __declspec(dllexport)
#else
#define BRUNSLI_DEC_API __declspec(dllimport)
#endif  /* BRUNSLIDEC_SHARED_COMPILATION */
#if defined(BRUNSLIENC_SHARED_COMPILATION)
#define BRUNSLI_ENC_API __declspec(dllexport)
#else
#define BRUNSLI_ENC_API __declspec(dllimport)
#endif  /* BRUNSLIENC_SHARED_COMPILATION */
#else  /* _WIN32 */
#define BRUNSLI_COMMON_API BRUNSLI_PUBLIC
#define BRUNSLI_DEC_API BRUNSLI_PUBLIC
#define BRUNSLI_ENC_API BRUNSLI_PUBLIC
#endif  /* _WIN32 */
#else  /* BRUNSLI_SHARED_COMPILATION */
#define BRUNSLI_COMMON_API
#define BRUNSLI_DEC_API
#define BRUNSLI_ENC_API
#endif

#endif  /* BRUNSLI_COMMON_PORT_H_ */
