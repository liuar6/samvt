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

#include <stdint.h>
#include <stdio.h>

#include "khash.h"
#include "vector.h"

#include "seq.h"

typedef struct fai_record_s{
    char *chrom;
    int char_num;
    int base_num;
    int32_t chrom_len;
    uint64_t offset;
} fai_record_t;

KHASH_MAP_INIT_STR(fai, fai_record_t *)
VEC_INIT(fai, fai_record_t *)

typedef struct fai_s{
    vec_t (fai) *fv;
    khash_t (fai) *fh;
} fai_t;


typedef struct fa_s{
    FILE *fp;
    fai_t *fai;
}fa_t;

static int strsplit(char * line, char ** results, int length, char c){
    char *start=line;
    char *end=NULL;
    int i=0, j;
    while ((end=strchr(start, c))!=NULL && i<length){
        end[0]='\0';
        results[i]=start;
        start=end+1;
        i=i+1;
    }
    if (i<length && start[0]!='\0') {
        results[i]=start;
        i=i+1;
    }
    for (j = i;j < length;++j) results[j]=NULL;
    return i;
}

fai_t *fai_open(const char *fname){
    char buffer[8192];
    char *field[6];
    fai_t *fai = malloc(sizeof(*fai));
    fai_record_t *r;
    vec_t (fai) *fv = vec_init(fai);

    FILE *fp = fopen(fname, "r");
    while (fgets(buffer, 8192, fp)){
        strsplit(buffer, field, 6, '\t');
        r = malloc(sizeof(*r));
        r->chrom = strdup(field[0]);
        r->chrom_len = strtol(field[1], NULL, 10);
        r->offset = strtol(field[2], NULL, 10);
        r->base_num = strtol(field[3], NULL, 10);
        r->char_num = strtol(field[4], NULL, 10);
        vec_add(fai, fv, r);
    }
    fclose(fp);

    khash_t (fai) *fh = kh_init(fai);
    for (int i = 0, ret; i < fv->size; ++i){
        khiter_t idx = kh_put(fai, fh, fv->data[i]->chrom, &ret);
        kh_val(fh, idx) = fv->data[i];
    }

    fai->fv = fv;
    fai->fh = fh;
    return fai;
}

void fai_close(fai_t *fai){
    vec_t (fai) *fv = fai->fv;
    for (int i = 0; i < fv->size; ++i){
        free(fv->data[i]->chrom);
        free(fv->data[i]);
    }
    vec_destroy(fai, fai->fv);
    kh_destroy(fai, fai->fh);
    free(fai);
}

fa_t *fa_open(const char *fname, const char *fai_name){
    fa_t *fa = malloc(sizeof(*fa));
    fa->fp = fopen(fname, "r");
    fa->fai = fai_open(fai_name);
    return fa;
}

void fa_close(fa_t *fa){
    fclose(fa->fp);
    fai_close(fa->fai);
    free(fa);
}

char *extract_sequence(fa_t *fa, char *chrom, int32_t start, int32_t end, char strand, char *seq, int32_t *msize){
    if (*msize < end - start + 1 || seq == NULL) {
        *msize = end - start + 1;
        seq = realloc(seq, *msize);
    }

    char *buf_start;
    uint64_t offset;
    int n, base_left;
    char sep[2];

    fai_record_t *r = kh_val(fa->fai->fh, kh_get(fai, fa->fai->fh, chrom));
    offset = r->offset + (start / r->base_num) * r->char_num + start % r->base_num;
    fseek(fa->fp, offset, SEEK_SET);

    buf_start = seq;
    base_left = end - start;

    n = fread(seq, 1, r->base_num - start % r->base_num < base_left ? r->base_num - start % r->base_num : base_left, fa->fp);
    base_left -= n;
    buf_start += n;
    while (base_left){
        fread(sep, 1, r->char_num - r->base_num, fa->fp);
        n=fread(buf_start, 1, base_left < r->base_num ? base_left : r->base_num, fa->fp);
        base_left -= n;
        buf_start += n;
        if (n == 0) break;
    }
    buf_start[0] = '\0';
    if (strand == '-') {
        seq_complement(seq);
        seq_reverse(seq);
    }
    return seq;
}