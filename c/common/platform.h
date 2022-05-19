// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

/* Macros for compiler / platform specific features and build options.

   Build options are:
    * BRUNSLI_BUILD_32_BIT disables 64-bit optimizations
    * BRUNSLI_BUILD_64_BIT forces to use 64-bit optimizations
    * BRUNSLI_BUILD_BIG_ENDIAN forces to use big-endian optimizations
    * BRUNSLI_BUILD_ENDIAN_NEUTRAL disables endian-aware optimizations
    * BRUNSLI_BUILD_LITTLE_ENDIAN forces to use little-endian optimizations
    * BRUNSLI_DEBUG enables "asserts" and extensive logging
    * BRUNSLI_DISABLE_LOG disables logging (useful for fuzzing)
*/

#ifndef BRUNSLI_COMMON_PLATFORM_H_
#define BRUNSLI_COMMON_PLATFORM_H_

#include <cstring>  /* memcpy */
#include <iomanip>
#include <ios>
#include <iostream>
#include <vector>

// Implicitly enable BRUNSLI_DEBUG when sanitizers are on.
#if !defined(BRUNSLI_DEBUG) && (BRUNSLI_SANITIZED || !defined(NDEBUG))
#define BRUNSLI_DEBUG 1
#endif

#if defined(BRUNSLI_USE_LOGGING)
#include "base/logging.h"
#else  // defined(BRUNSLI_USE_LOGGING)
#include <stdio.h>
#endif  // defined(BRUNSLI_USE_LOGGING)

#include "./port.h"
#include <brunsli/types.h>

#if defined(OS_LINUX) || defined(OS_CYGWIN)
#include <endian.h>
#elif defined(OS_FREEBSD)
#include <machine/endian.h>
#elif defined(OS_MACOSX)
#include <machine/endian.h>
/* Let's try and follow the Linux convention */
#define BRUNSLI_X_BYTE_ORDER BYTE_ORDER
#define BRUNSLI_X_LITTLE_ENDIAN LITTLE_ENDIAN
#define BRUNSLI_X_BIG_ENDIAN BIG_ENDIAN
#endif

/* The following macros were borrowed from https://github.com/nemequ/hedley
 * with permission of original author - Evan Nemerson <evan@nemerson.com> */

/* >>> >>> >>> hedley macros */

/* Define "BRUNSLI_PREDICT_TRUE" and "BRUNSLI_PREDICT_FALSE" macros for capable
   compilers.

To apply compiler hint, enclose the branching condition into macros, like this:

  if (BRUNSLI_PREDICT_TRUE(zero == 0)) {
    // main execution path
  } else {
    // compiler should place this code outside of main execution path
  }

OR:

  if (BRUNSLI_PREDICT_FALSE(something_rare_or_unexpected_happens)) {
    // compiler should place this code outside of main execution path
  }

*/
#if BRUNSLI_GNUC_HAS_BUILTIN(__builtin_expect, 3, 0, 0) || \
    BRUNSLI_INTEL_VERSION_CHECK(16, 0, 0) ||               \
    BRUNSLI_SUNPRO_VERSION_CHECK(5, 12, 0) ||              \
    BRUNSLI_ARM_VERSION_CHECK(4, 1, 0) ||                  \
    BRUNSLI_IBM_VERSION_CHECK(10, 1, 0) ||                 \
    BRUNSLI_TI_VERSION_CHECK(7, 3, 0) ||                   \
    BRUNSLI_TINYC_VERSION_CHECK(0, 9, 27)
#define BRUNSLI_PREDICT_TRUE(x) (__builtin_expect(!!(x), 1))
#define BRUNSLI_PREDICT_FALSE(x) (__builtin_expect(x, 0))
#else
#define BRUNSLI_PREDICT_FALSE(x) (x)
#define BRUNSLI_PREDICT_TRUE(x) (x)
#endif

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L) && \
    !defined(__cplusplus)
#define BRUNSLI_RESTRICT restrict
#elif BRUNSLI_GNUC_VERSION_CHECK(3, 1, 0) ||                         \
    BRUNSLI_MSVC_VERSION_CHECK(14, 0, 0) ||                          \
    BRUNSLI_INTEL_VERSION_CHECK(16, 0, 0) ||                         \
    BRUNSLI_ARM_VERSION_CHECK(4, 1, 0) ||                            \
    BRUNSLI_IBM_VERSION_CHECK(10, 1, 0) ||                           \
    BRUNSLI_PGI_VERSION_CHECK(17, 10, 0) ||                          \
    BRUNSLI_TI_VERSION_CHECK(8, 0, 0) ||                             \
    BRUNSLI_IAR_VERSION_CHECK(8, 0, 0) ||                            \
    (BRUNSLI_SUNPRO_VERSION_CHECK(5, 14, 0) && defined(__cplusplus))
#define BRUNSLI_RESTRICT __restrict
#elif BRUNSLI_SUNPRO_VERSION_CHECK(5, 3, 0) && !defined(__cplusplus)
#define BRUNSLI_RESTRICT _Restrict
#else
#define BRUNSLI_RESTRICT
#endif

#if (defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || \
    (defined(__cplusplus) && (__cplusplus >= 199711L))
#define BRUNSLI_MAYBE_INLINE inline
#elif defined(__GNUC_STDC_INLINE__) || defined(__GNUC_GNU_INLINE__) || \
    BRUNSLI_ARM_VERSION_CHECK(6, 2, 0)
#define BRUNSLI_MAYBE_INLINE __inline__
#elif BRUNSLI_MSVC_VERSION_CHECK(12, 0, 0) || \
    BRUNSLI_ARM_VERSION_CHECK(4, 1, 0) || BRUNSLI_TI_VERSION_CHECK(8, 0, 0)
#define BRUNSLI_MAYBE_INLINE __inline
#else
#define BRUNSLI_MAYBE_INLINE
#endif

#if BRUNSLI_GNUC_HAS_ATTRIBUTE(always_inline, 4, 0, 0) ||                      \
    BRUNSLI_INTEL_VERSION_CHECK(16, 0, 0) ||                                   \
    BRUNSLI_SUNPRO_VERSION_CHECK(5, 11, 0) ||                                  \
    BRUNSLI_ARM_VERSION_CHECK(4, 1, 0) ||                                      \
    BRUNSLI_IBM_VERSION_CHECK(10, 1, 0) ||                                     \
    BRUNSLI_TI_VERSION_CHECK(8, 0, 0) ||                                       \
    (BRUNSLI_TI_VERSION_CHECK(7, 3, 0) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__))
#define BRUNSLI_INLINE BRUNSLI_MAYBE_INLINE __attribute__((__always_inline__))
#elif BRUNSLI_MSVC_VERSION_CHECK(12, 0, 0)
#define BRUNSLI_INLINE BRUNSLI_MAYBE_INLINE __forceinline
#elif BRUNSLI_TI_VERSION_CHECK(7, 0, 0) && defined(__cplusplus)
#define BRUNSLI_INLINE BRUNSLI_MAYBE_INLINE _Pragma("FUNC_ALWAYS_INLINE;")
#elif BRUNSLI_IAR_VERSION_CHECK(8, 0, 0)
#define BRUNSLI_INLINE BRUNSLI_MAYBE_INLINE _Pragma("inline=forced")
#else
#define BRUNSLI_INLINE BRUNSLI_MAYBE_INLINE
#endif

#if BRUNSLI_GNUC_HAS_ATTRIBUTE(noinline, 4, 0, 0) ||                           \
    BRUNSLI_INTEL_VERSION_CHECK(16, 0, 0) ||                                   \
    BRUNSLI_SUNPRO_VERSION_CHECK(5, 11, 0) ||                                  \
    BRUNSLI_ARM_VERSION_CHECK(4, 1, 0) ||                                      \
    BRUNSLI_IBM_VERSION_CHECK(10, 1, 0) ||                                     \
    BRUNSLI_TI_VERSION_CHECK(8, 0, 0) ||                                       \
    (BRUNSLI_TI_VERSION_CHECK(7, 3, 0) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__))
#define BRUNSLI_NOINLINE __attribute__((__noinline__))
#elif BRUNSLI_MSVC_VERSION_CHECK(13, 10, 0)
#define BRUNSLI_NOINLINE __declspec(noinline)
#elif BRUNSLI_PGI_VERSION_CHECK(10, 2, 0)
#define BRUNSLI_NOINLINE _Pragma("noinline")
#elif BRUNSLI_TI_VERSION_CHECK(6, 0, 0) && defined(__cplusplus)
#define BRUNSLI_NOINLINE _Pragma("FUNC_CANNOT_INLINE;")
#elif BRUNSLI_IAR_VERSION_CHECK(8, 0, 0)
#define BRUNSLI_NOINLINE _Pragma("inline=never")
#else
#define BRUNSLI_NOINLINE
#endif

/* BRUNSLI_INTERNAL could be defined to override visibility, e.g. for tests. */
#if !defined(BRUNSLI_INTERNAL)
#if defined(_WIN32) || defined(__CYGWIN__)
#define BRUNSLI_INTERNAL
#elif BRUNSLI_GNUC_VERSION_CHECK(3, 3, 0) ||                        \
    BRUNSLI_TI_VERSION_CHECK(8, 0, 0) ||                            \
    BRUNSLI_INTEL_VERSION_CHECK(16, 0, 0) ||                        \
    BRUNSLI_ARM_VERSION_CHECK(4, 1, 0) ||                           \
    BRUNSLI_IBM_VERSION_CHECK(13, 1, 0) ||                          \
    BRUNSLI_SUNPRO_VERSION_CHECK(5, 11, 0) ||                       \
    (BRUNSLI_TI_VERSION_CHECK(7, 3, 0) &&                           \
     defined(__TI_GNU_ATTRIBUTE_SUPPORT__) && defined(__TI_EABI__))
#define BRUNSLI_INTERNAL __attribute__ ((visibility ("hidden")))
#else
#define BRUNSLI_INTERNAL
#endif
#endif

/* <<< <<< <<< end of hedley macros. */

#if BRUNSLI_GNUC_HAS_ATTRIBUTE(unused, 2, 7, 0) || \
    BRUNSLI_INTEL_VERSION_CHECK(16, 0, 0)
#define BRUNSLI_UNUSED_FUNCTION static BRUNSLI_INLINE __attribute__ ((unused))
#else
#define BRUNSLI_UNUSED_FUNCTION static BRUNSLI_INLINE
#endif

#if BRUNSLI_GNUC_HAS_ATTRIBUTE(aligned, 2, 7, 0)
#define BRUNSLI_ALIGNED(N) __attribute__((aligned(N)))
#else
#define BRUNSLI_ALIGNED(N)
#endif

#if (defined(__ARM_ARCH) && (__ARM_ARCH == 7)) || \
    (defined(M_ARM) && (M_ARM == 7))
#define BRUNSLI_TARGET_ARMV7
#endif  /* ARMv7 */

#if (defined(__ARM_ARCH) && (__ARM_ARCH == 8)) || \
    defined(__aarch64__) || defined(__ARM64_ARCH_8__)
#define BRUNSLI_TARGET_ARMV8_ANY

#if defined(__ARM_32BIT_STATE)
#define BRUNSLI_TARGET_ARMV8_32
#elif defined(__ARM_64BIT_STATE)
#define BRUNSLI_TARGET_ARMV8_64
#endif

#endif  /* ARMv8 */

#if defined(__ARM_NEON__) || defined(__ARM_NEON)
#define BRUNSLI_TARGET_NEON
#endif

#if defined(__i386) || defined(_M_IX86)
#define BRUNSLI_TARGET_X86
#endif

#if defined(__x86_64__) || defined(_M_X64)
#define BRUNSLI_TARGET_X64
#endif

#if defined(__PPC64__)
#define BRUNSLI_TARGET_POWERPC64
#endif

#if defined(__riscv) && defined(__riscv_xlen) && __riscv_xlen == 64
#define BRUNSLI_TARGET_RISCV64
#endif

#if defined(BRUNSLI_BUILD_64_BIT)
#define BRUNSLI_64_BITS 1
#elif defined(BRUNSLI_BUILD_32_BIT)
#define BRUNSLI_64_BITS 0
#elif defined(BRUNSLI_TARGET_X64) || defined(BRUNSLI_TARGET_ARMV8_64) || \
    defined(BRUNSLI_TARGET_POWERPC64) || defined(BRUNSLI_TARGET_RISCV64)
#define BRUNSLI_64_BITS 1
#else
#define BRUNSLI_64_BITS 0
#endif

#if (BRUNSLI_64_BITS)
#define brunsli_reg_t uint64_t
#else
#define brunsli_reg_t uint32_t
#endif

#if defined(BRUNSLI_BUILD_BIG_ENDIAN)
#define BRUNSLI_BIG_ENDIAN 1
#elif defined(BRUNSLI_BUILD_LITTLE_ENDIAN)
#define BRUNSLI_LITTLE_ENDIAN 1
#elif defined(BRUNSLI_BUILD_ENDIAN_NEUTRAL)
/* Just break elif chain. */
#elif defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
#define BRUNSLI_LITTLE_ENDIAN 1
#elif defined(_WIN32) || defined(BRUNSLI_TARGET_X64)
/* Win32 & x64 can currently always be assumed to be little endian */
#define BRUNSLI_LITTLE_ENDIAN 1
#elif defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
#define BRUNSLI_BIG_ENDIAN 1
#elif defined(BRUNSLI_X_BYTE_ORDER)
#if BRUNSLI_X_BYTE_ORDER == BRUNSLI_X_LITTLE_ENDIAN
#define BRUNSLI_LITTLE_ENDIAN 1
#elif BRUNSLI_X_BYTE_ORDER == BRUNSLI_X_BIG_ENDIAN
#define BRUNSLI_BIG_ENDIAN 1
#endif
#endif  /* BRUNSLI_X_BYTE_ORDER */

#if !defined(BRUNSLI_LITTLE_ENDIAN)
#define BRUNSLI_LITTLE_ENDIAN 0
#endif

#if !defined(BRUNSLI_BIG_ENDIAN)
#define BRUNSLI_BIG_ENDIAN 0
#endif

#if defined(BRUNSLI_X_BYTE_ORDER)
#undef BRUNSLI_X_BYTE_ORDER
#undef BRUNSLI_X_LITTLE_ENDIAN
#undef BRUNSLI_X_BIG_ENDIAN
#endif

/* Portable unaligned memory access: read / write values via memcpy. */
static BRUNSLI_INLINE uint16_t BrunsliUnalignedRead16(const void* p) {
  uint16_t t;
  memcpy(&t, p, sizeof t);
  return t;
}
/* Portable unaligned memory access: read / write values via memcpy. */
static BRUNSLI_INLINE void BrunsliUnalignedWrite16(void* p, uint16_t v) {
  memcpy(p, &v, sizeof v);
}
static BRUNSLI_INLINE uint32_t BrunsliUnalignedRead32(const void* p) {
  uint32_t t;
  memcpy(&t, p, sizeof t);
  return t;
}
static BRUNSLI_INLINE uint64_t BrunsliUnalignedRead64(const void* p) {
  uint64_t t;
  memcpy(&t, p, sizeof t);
  return t;
}
static BRUNSLI_INLINE void BrunsliUnalignedWrite64(void* p, uint64_t v) {
  memcpy(p, &v, sizeof v);
}

#if BRUNSLI_LITTLE_ENDIAN
/* Straight endianness. Just read / write values. */
#define BRUNSLI_UNALIGNED_LOAD16LE BrunsliUnalignedRead16
#define BRUNSLI_UNALIGNED_STORE16LE BrunsliUnalignedWrite16
#define BRUNSLI_UNALIGNED_LOAD32LE BrunsliUnalignedRead32
#define BRUNSLI_UNALIGNED_LOAD64LE BrunsliUnalignedRead64
#define BRUNSLI_UNALIGNED_STORE64LE BrunsliUnalignedWrite64
#elif BRUNSLI_BIG_ENDIAN  /* BRUNSLI_LITTLE_ENDIAN */
/* Explain compiler to byte-swap values. */
#define BRUNSLI_BSWAP16_(V) ((uint16_t)( \
  (((V) & 0xFFU) << 8) | \
  (((V) >> 8) & 0xFFU)))
static BRUNSLI_INLINE uint16_t BRUNSLI_UNALIGNED_LOAD16LE(const void* p) {
  uint16_t value = BrunsliUnalignedRead16(p);
  return BRUNSLI_BSWAP16_(value);
}
static BRUNSLI_INLINE void BRUNSLI_UNALIGNED_STORE16LE(void* p, uint16_t v) {
  uint16_t value = BRUNSLI_BSWAP16_(v);
  BrunsliUnalignedWrite16(p, value);
}
#define BRUNSLI_BSWAP32_(V) ( \
  (((V) & 0xFFU) << 24) | (((V) & 0xFF00U) << 8) | \
  (((V) >> 8) & 0xFF00U) | (((V) >> 24) & 0xFFU))
static BRUNSLI_INLINE uint32_t BRUNSLI_UNALIGNED_LOAD32LE(const void* p) {
  uint32_t value = BrunsliUnalignedRead32(p);
  return BRUNSLI_BSWAP32_(value);
}
#define BRUNSLI_BSWAP64_(V) ( \
  (((V) & 0xFFU) << 56) | (((V) & 0xFF00U) << 40) | \
  (((V) & 0xFF0000U) << 24) | (((V) & 0xFF000000U) << 8) | \
  (((V) >> 8) & 0xFF000000U) | (((V) >> 24) & 0xFF0000U) | \
  (((V) >> 40) & 0xFF00U) | (((V) >> 56) & 0xFFU))
static BRUNSLI_INLINE uint64_t BRUNSLI_UNALIGNED_LOAD64LE(const void* p) {
  uint64_t value = BrunsliUnalignedRead64(p);
  return BRUNSLI_BSWAP64_(value);
}
static BRUNSLI_INLINE void BRUNSLI_UNALIGNED_STORE64LE(void* p, uint64_t v) {
  uint64_t value = BRUNSLI_BSWAP64_(v);
  BrunsliUnalignedWrite64(p, value);
}
#else  /* BRUNSLI_LITTLE_ENDIAN */
/* Read / store values byte-wise; hopefully compiler will understand. */
static BRUNSLI_INLINE uint16_t BRUNSLI_UNALIGNED_LOAD16LE(const void* p) {
  const uint8_t* in = (const uint8_t*)p;
  return (uint16_t)(in[0] | (in[1] << 8));
}
static BRUNSLI_INLINE void BRUNSLI_UNALIGNED_STORE16LE(void* p, uint16_t v) {
  uint8_t* out = (uint8_t*)p;
  out[0] = (uint8_t)v;
  out[1] = (uint8_t)(v >> 8);
}
static BRUNSLI_INLINE uint32_t BRUNSLI_UNALIGNED_LOAD32LE(const void* p) {
  const uint8_t* in = (const uint8_t*)p;
  uint32_t value = (uint32_t)(in[0]);
  value |= (uint32_t)(in[1]) << 8;
  value |= (uint32_t)(in[2]) << 16;
  value |= (uint32_t)(in[3]) << 24;
  return value;
}
static BRUNSLI_INLINE uint64_t BRUNSLI_UNALIGNED_LOAD64LE(const void* p) {
  const uint8_t* in = (const uint8_t*)p;
  uint64_t value = (uint64_t)(in[0]);
  value |= (uint64_t)(in[1]) << 8;
  value |= (uint64_t)(in[2]) << 16;
  value |= (uint64_t)(in[3]) << 24;
  value |= (uint64_t)(in[4]) << 32;
  value |= (uint64_t)(in[5]) << 40;
  value |= (uint64_t)(in[6]) << 48;
  value |= (uint64_t)(in[7]) << 56;
  return value;
}
static BRUNSLI_INLINE void BRUNSLI_UNALIGNED_STORE64LE(void* p, uint64_t v) {
  uint8_t* out = (uint8_t*)p;
  out[0] = (uint8_t)v;
  out[1] = (uint8_t)(v >> 8);
  out[2] = (uint8_t)(v >> 16);
  out[3] = (uint8_t)(v >> 24);
  out[4] = (uint8_t)(v >> 32);
  out[5] = (uint8_t)(v >> 40);
  out[6] = (uint8_t)(v >> 48);
  out[7] = (uint8_t)(v >> 56);
}
#endif  /* BRUNSLI_LITTLE_ENDIAN */

/* BRUNSLI_IS_CONSTANT macros returns true for compile-time constants. */
#if BRUNSLI_GNUC_HAS_BUILTIN(__builtin_constant_p, 3, 0, 1) || \
    BRUNSLI_INTEL_VERSION_CHECK(16, 0, 0)
#define BRUNSLI_IS_CONSTANT(x) (!!__builtin_constant_p(x))
#else
#define BRUNSLI_IS_CONSTANT(x) (!!0)
#endif

#if defined(BRUNSLI_TARGET_ARMV7) || defined(BRUNSLI_TARGET_ARMV8_ANY)
#define BRUNSLI_HAS_UBFX (!!1)
#else
#define BRUNSLI_HAS_UBFX (!!0)
#endif

// "else" branch is never evaluated, but provides the sink.
#define BRUNSLI_VOID_LOG() if (true) {} else std::cerr

// This macro allows easy logging engine replacement.
#if defined(BRUNSLI_USE_LOGGING)
#define BRUNSLI_LOG_(LEVEL) LOG(LEVEL)
#define BRUNSLI_ENDL() ""
#else  // defined(BRUNSLI_USE_LOGGING)
#define BRUNSLI_LOG_(LEVEL) std::cerr
#define BRUNSLI_ENDL() std::endl
#endif  // defined(BRUNSLI_USE_LOGGING)

#if defined(BRUNSLI_DISABLE_LOG)
#define BRUNSLI_LOG_DEBUG() BRUNSLI_VOID_LOG()
#define BRUNSLI_LOG_INFO() BRUNSLI_VOID_LOG()
#define BRUNSLI_LOG_WARNING() BRUNSLI_VOID_LOG()
#define BRUNSLI_LOG_ERROR() BRUNSLI_VOID_LOG()
#else  // defined(BRUNSLI_DISABLE_LOG)
// TODO(eustas): get rid of base/logging.h dependency
#if defined(BRUNSLI_ENABLE_LOG)
#define BRUNSLI_LOG_DEBUG() BRUNSLI_LOG_(INFO)
#else  //  defined(BRUNSLI_ENABLE_LOG)
#define BRUNSLI_LOG_DEBUG() BRUNSLI_VOID_LOG()
#endif  //  defined(BRUNSLI_ENABLE_LOG)
#define BRUNSLI_LOG_INFO() BRUNSLI_LOG_(INFO)
#define BRUNSLI_LOG_WARNING() BRUNSLI_LOG_(WARNING)
#define BRUNSLI_LOG_ERROR() BRUNSLI_LOG_(ERROR)
#endif  // defined(BRUNSLI_DISABLE_LOG)

namespace brunsli {
void BrunsliDumpAndAbort(const char* f, int l, const char* fn);
static BRUNSLI_INLINE void Append(std::vector<uint8_t>* dst,
                                  const uint8_t* begin, const uint8_t* end) {
  dst->insert(dst->end(), begin, end);
}
static BRUNSLI_INLINE void Append(std::vector<uint8_t>* dst,
                                  const uint8_t* begin, size_t length) {
  Append(dst, begin, begin + length);
}
static BRUNSLI_INLINE void Append(std::vector<uint8_t>* dst,
                                  const std::vector<uint8_t>& src) {
  Append(dst, src.data(), src.size());
}
}  // namespace brunsli

// TODO(eustas): use "predict false" to move the code out of hot path.
#define BRUNSLI_CHECK(V) \
  if (!(V)) {                                                         \
    ::brunsli::BrunsliDumpAndAbort(__FILE__, __LINE__, __FUNCTION__); \
    /* Tell the compiler, that there is no escape route. */           \
    while (true) ;                                                    \
  }

#if defined(BRUNSLI_DEBUG)
#define BRUNSLI_DCHECK(V) BRUNSLI_CHECK(V)
#else
#define BRUNSLI_DCHECK(V)
#endif

// TODO(eustas): Pick up upgrade after https://github.com/google/brotli/pull/636
//               is landed and merged.
inline int Log2FloorNonZero(uint32_t n) {
#ifdef __GNUC__
  return 31 ^ __builtin_clz(n);
#else
  unsigned int result = 0;
  while (n >>= 1) result++;
  return result;
#endif
}

#define BRUNSLI_UNUSED(X) (void)(X)

BRUNSLI_UNUSED_FUNCTION void BrunsliSuppressUnusedFunctions(void) {
  BRUNSLI_UNUSED(
      static_cast<void (*)(std::vector<uint8_t>*, const std::vector<uint8_t>&)>(
          &brunsli::Append));
  BRUNSLI_UNUSED(&BrunsliSuppressUnusedFunctions);
  BRUNSLI_UNUSED(&BrunsliUnalignedRead16);
  BRUNSLI_UNUSED(&BrunsliUnalignedWrite16);
  BRUNSLI_UNUSED(&BrunsliUnalignedRead32);
  BRUNSLI_UNUSED(&BrunsliUnalignedRead64);
  BRUNSLI_UNUSED(&BrunsliUnalignedWrite64);
  BRUNSLI_UNUSED(&BRUNSLI_UNALIGNED_LOAD16LE);
  BRUNSLI_UNUSED(&BRUNSLI_UNALIGNED_STORE16LE);
  BRUNSLI_UNUSED(&BRUNSLI_UNALIGNED_LOAD32LE);
  BRUNSLI_UNUSED(&BRUNSLI_UNALIGNED_LOAD64LE);
  BRUNSLI_UNUSED(&BRUNSLI_UNALIGNED_STORE64LE);
}

#endif  // BRUNSLI_COMMON_PLATFORM_H_
