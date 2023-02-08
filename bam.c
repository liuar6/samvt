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

#include "bam.h"

bt_bam_t *bt_bam_open(const char* fn, int n_threads){
    bt_bam_t *s = malloc(sizeof(bt_bam_t));
    s->fn = strdup(fn);
    s->fp = sam_open(fn, "r");
    s->s = NULL;
    if (s->fp->is_bgzf) bgzf_mt(s->fp->fp.bgzf, n_threads, 128);
    else s->s = mt_server_init(n_threads);
    s->hdr = sam_hdr_read(s->fp);
    return s;
}

int bt_bam_close(bt_bam_t *s) {
    free(s->fn);
    bam_hdr_destroy(s->hdr);
    sam_close(s->fp);
    if (s->s) mt_server_destroy(s->s);
    free(s);
    return 0;
}

int bt_bam_next(bt_bam_t *s, bam1_t *b) {
    if (sam_read1(s->fp, s->hdr, b) >= 0) return 0;
    else return -1;
}

int bt_bam_next2(bt_bam_t *s, bam1_t *b1, bam1_t *b2) {
    /* each time, return the pair end records */
    return 0;
}

mt_server *bt_bam_mt_server(bt_bam_t *s) {
    if (s->fp->is_bgzf) return (mt_server *) s->fp->fp.bgzf->mt->pool;
    else return s->s;
}