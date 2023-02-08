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

#ifndef SAMVT_SAM_H
#define SAMVT_SAM_H

#include "htslib/bgzf.h"
#include "htslib/sam.h"

#include "htslib_bgzf.h"
#include "mt.h"

typedef struct bt_bam_s{
    char *fn;
    samFile *fp;
    bam_hdr_t *hdr;
    mt_server *s;
} bt_bam_t;

bt_bam_t *bt_bam_open(const char* fn, int n_threads);
int bt_bam_close(bt_bam_t *s);
int bt_bam_next(bt_bam_t *s, bam1_t *b);
int bt_bam_next2(bt_bam_t *s, bam1_t *b1, bam1_t *b2);
mt_server *bt_bam_mt_server(bt_bam_t *s);


#endif //SAMVT_SAM_H
