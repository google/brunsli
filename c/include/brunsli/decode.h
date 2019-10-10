/*
Copyright (c) Google LLC 2019

Use of this source code is governed by an MIT-style
license that can be found in the LICENSE file or at
https://opensource.org/licenses/MIT.
*/

#ifndef BRUNSLI_DEC_DECODE_H_
#define BRUNSLI_DEC_DECODE_H_

#include <brunsli/types.h>

/* C API for brunsli decoder */

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

typedef size_t (*DecodeBrunsliSink)(void* out_data, const uint8_t* buf,
                                    size_t size);

/*
Decodes brunsli file to JPEG. Returns 1 on success, 0 on error.
The input data must be complete (decodes in one shot).
Outputs to out_fun, out_fun must return amount of consumed bytes, any return
value not equal to the input size is considered an error. It will pass on the
out_data to out_fun.
*/
int DecodeBrunsli(size_t in_size, const uint8_t* in, void* out_data,
                  DecodeBrunsliSink out_fun);

#if defined(__cplusplus) || defined(c_plusplus)
} /* extern "C" */
#endif

#endif /* BRUNSLI_DEC_DECODE_H_ */
