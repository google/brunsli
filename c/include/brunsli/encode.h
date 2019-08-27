/*
Copyright (c) Google LLC 2019

Use of this source code is governed by an MIT-style
license that can be found in the LICENSE file or at
https://opensource.org/licenses/MIT.
*/

#ifndef BRUNSLI_ENC_ENCODE_H_
#define BRUNSLI_ENC_ENCODE_H_

#include <string.h>

/* C API for brunsli encoder */

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

/*
Encodes brunsli file to JPEG. Returns 1 on success, 0 on error.
The input data must be complete (encodes in one shot).
Outputs to outfun, outfun must return amount of consumed bytes, any return value
not equal to the input size is considered an error. It will pass on the outdata
to outfun.
*/
int EncodeBrunsli(size_t insize, const unsigned char* in, void* outdata,
    size_t (*outfun)(void* outdata, const unsigned char* buf, size_t size));

#if defined(__cplusplus) || defined(c_plusplus)
}  /* extern "C" */
#endif

#endif  /* BRUNSLI_ENC_ENCODE_H_ */
