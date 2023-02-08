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

typedef struct coverage_s coverage_t;
coverage_t *coverage_init(int32_t n_targets, char **target_name, uint32_t *target_len, uint32_t coverage_block_shift);
coverage_t *coverage_mt(coverage_t *cov);
int coverage_destroy(coverage_t * cov);
int coverage_update(coverage_t *cov, int32_t target, uint32_t start, uint32_t end);
typedef struct coverage2_s coverage2_t;
coverage2_t *coverage2_init(int32_t n_targets, char **target_name, uint32_t *target_len, uint32_t coverage_block_shift);
coverage2_t *coverage2_mt(coverage2_t *cov);
int coverage2_destroy(coverage2_t * cov);
int coverage2_update(coverage2_t *cov, int32_t target, uint32_t start, uint32_t end, char strand, uint8_t *seq, int rpos);
int output_bw(coverage_t *cov, char *fn, mt_server *s);