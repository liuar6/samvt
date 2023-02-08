/* The MIT License (MIT)

   Copyright (c) 2023 Anrui Liu <liuar6@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   “Software”), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
 */

#ifndef FA_LIBRARY_H
#define FA_LIBRARY_H

#include <stdint.h>
#include <stdio.h>
#include "khash.h"
#include "vector.h"

typedef struct fai_record_s{
    char *chrom;
    int char_num;
    int base_num;
    int32_t chrom_len;
    uint64_t offset;
} fai_record_t;

VEC_DECLARE(fai, fai_record_t *)
KHASH_DECLARE(fai, char *, fai_record_t *)

typedef struct fai_s{
    vec_t (fai) *fv;
    khash_t (fai) *fh;
} fai_t;

typedef struct fa_s{
    FILE *fp;
    fai_t *fai;
}fa_t;

fa_t *fa_open(const char *fname, const char *fai_name);
void fa_close(fa_t *fa);
char *extract_sequence(fa_t *fa, char *chrom, int32_t start, int32_t end, char strand, char *seq, int32_t *msize);

#endif //FA_LIBRARY_H