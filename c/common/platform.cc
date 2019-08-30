// Copyright (c) Google LLC 2019
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "./platform.h"

#include <cstdio>
#include <cstdlib>  // for abort

namespace brunsli {

void BrunsliDumpAndAbort(const char* f, int l, const char* fn) {
  fprintf(stderr, "%s:%d (%s)\n", f, l, fn);
  fflush(stderr);
  abort();
}

}  // namespace brunsli
