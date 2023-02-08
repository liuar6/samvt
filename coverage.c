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

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <pthread.h>

#include "htslib/bgzf.h"
#include "htslib/sam.h"

#include "bigWig.h"

#include "bam.h"
#include "mt.h"
#include "mt_buffer.h"
#include "coverage.h"

#define cov_val_t double
typedef struct coverage_s{
    int32_t n_targets;
    char **target_name;
    uint32_t *target_len;
    uint32_t bin_size;
    uint32_t coverage_block_shift;
    uint32_t coverage_block_size;
    cov_val_t ***coverage_blocks;
    int is_mt;
    uint32_t coverage_mutex_shift;
    pthread_mutex_t **coverage_block_mutexes;
} coverage_t;

coverage_t *coverage_init(int32_t n_targets, char **target_name, uint32_t *target_len, uint32_t coverage_block_shift){
    coverage_t *cov = calloc(1, sizeof(coverage_t));
    cov->n_targets = n_targets;
    cov->target_name = calloc(n_targets, sizeof(char *));
    cov->target_len = calloc(n_targets, sizeof(uint32_t));
    cov->bin_size = 1;
    cov->coverage_block_shift =  coverage_block_shift;
    cov->coverage_block_size = 1u<<cov->coverage_block_shift;
    cov->coverage_blocks = calloc(n_targets, sizeof(cov_val_t *));
    for (int i = 0; i < cov->n_targets; ++i) {
        cov->target_name[i] = strdup(target_name[i]);
        cov->target_len[i] = target_len[i];
        int32_t bin_count = (cov->target_len[i]-1)/cov->bin_size+1;
        int32_t block_count =  ((bin_count-1)/cov->coverage_block_size)+1;
        cov->coverage_blocks[i] =  calloc(block_count, sizeof(cov_val_t*));
    }
    return cov;
}

coverage_t *coverage_mt(coverage_t *cov){
    cov->is_mt = 1;
    cov->coverage_mutex_shift = 1;
    cov->coverage_block_mutexes = calloc(cov->n_targets, sizeof(pthread_mutex_t *));
    for (int i = 0; i < cov->n_targets; ++i) {
        int32_t block_count =  ((cov->target_len[i]-1)/cov->coverage_block_size)+1;
        int32_t block_mutex_count = ((block_count-1u)>>cov->coverage_mutex_shift)+1;
        cov->coverage_block_mutexes[i] =  calloc(block_mutex_count, sizeof(pthread_mutex_t));
        for (int j = 0; j < block_mutex_count; ++j) pthread_mutex_init(&cov->coverage_block_mutexes[i][j], NULL);
    }
    return cov;
}

int coverage_destroy(coverage_t * cov){
    for (int i = 0; i < cov->n_targets; ++i) {
        uint32_t block_count = ((cov->target_len[i]-1)>>cov->coverage_block_shift)+1;
        for (int j = 0; j < block_count; ++j) if (cov->coverage_blocks[i][j]!=NULL) free(cov->coverage_blocks[i][j]);
        free(cov->coverage_blocks[i]);
        if (cov->is_mt) {
            uint32_t block_mutex_count = ((block_count-1)>>cov->coverage_mutex_shift)+1;
            for (int j = 0; j < block_mutex_count; ++j) pthread_mutex_destroy(&cov->coverage_block_mutexes[i][j]);
            free(cov->coverage_block_mutexes[i]);
        }
    }
    free(cov->coverage_blocks);
    if (cov->is_mt) free(cov->coverage_block_mutexes);
    for (int i=0; i < cov->n_targets; ++i) free(cov->target_name[i]);
    free(cov->target_name);
    free(cov->target_len);
    free(cov);
    return 0;
}

/* starts are 0-based and ends are 1-based */
int coverage_update(coverage_t *cov, int32_t target, uint32_t start, uint32_t end){
    /* this function simply trust its arguments without checking whether target is present and coordinate is in valid range */
    uint32_t block_index_start, block_index_end, block_index, block_start, block_end, new_start, new_end;
    cov_val_t *coverage_block;
    cov_val_t **coverage_block_target;
    pthread_mutex_t *mutex;
    block_index_start = start/cov->coverage_block_size;
    block_index_end = (end-1)/cov->coverage_block_size;
    while (block_index_start <= block_index_end){
        block_index = block_index_start;
        block_start = block_index<<cov->coverage_block_shift;
        new_start = (start > block_start) ? start - block_start : 0;
        new_end = (end - block_start > cov->coverage_block_size) ? cov->coverage_block_size : end - block_start;
        coverage_block_target = cov->coverage_blocks[target];
        if (cov->is_mt) {
            mutex=&cov->coverage_block_mutexes[target][block_index>>cov->coverage_mutex_shift];
            pthread_mutex_lock(mutex);
        }
        coverage_block = coverage_block_target[block_index];
        if (coverage_block == NULL) {
            uint32_t target_len = cov->target_len[target];
            int needed=cov->coverage_block_size > target_len - block_start ? target_len - block_start : cov->coverage_block_size ;
            coverage_block = calloc(needed, sizeof(cov_val_t));
            cov->coverage_blocks[target][block_index] = coverage_block;
        }
        while (new_start < new_end) coverage_block[new_start++]++;
        if (cov->is_mt) pthread_mutex_unlock(mutex);
        block_index_start++;
    }
    return 0;
}

typedef struct {
    double A;
    double C;
    double G;
    double T;
    double N;
} cov2_val_t;

typedef struct coverage2_s{
    int32_t n_targets;
    char **target_name;
    uint32_t *target_len;
    uint32_t coverage_block_shift;
    uint32_t coverage_block_size;
    cov2_val_t ***coverage_blocks;
    int is_mt;
    uint32_t coverage_mutex_shift;
    pthread_mutex_t **coverage_block_mutexes;
} coverage2_t;

coverage2_t *coverage2_init(int32_t n_targets, char **target_name, uint32_t *target_len, uint32_t coverage_block_shift){
    coverage2_t *cov = calloc(1, sizeof(coverage2_t));
    cov->n_targets = n_targets;
    cov->target_name = calloc(n_targets, sizeof(char *));
    cov->target_len = calloc(n_targets, sizeof(uint32_t));
    cov->coverage_block_shift =  coverage_block_shift;
    cov->coverage_block_size = 1u<<cov->coverage_block_shift;
    cov->coverage_blocks = calloc(n_targets * 2, sizeof(cov2_val_t *));
    for (int i = 0; i < cov->n_targets; ++i) {
        cov->target_name[i] = strdup(target_name[i]);
        cov->target_len[i] = target_len[i];
        int32_t bin_count = (cov->target_len[i]-1)+1;
        int32_t block_count =  ((bin_count-1)/cov->coverage_block_size)+1;
        cov->coverage_blocks[i] =  calloc(block_count, sizeof(cov_val_t*));
    }
    for (int i = cov->n_targets; i < cov->n_targets * 2; ++i){
        int32_t bin_count = (cov->target_len[i - cov->n_targets]-1)+1;
        int32_t block_count =  ((bin_count-1)/cov->coverage_block_size)+1;
        cov->coverage_blocks[i] =  calloc(block_count, sizeof(cov2_val_t*));
    }
    return cov;
}

coverage2_t *coverage2_mt(coverage2_t *cov){
    cov->is_mt = 1;
    cov->coverage_mutex_shift = 1;
    cov->coverage_block_mutexes = calloc(cov->n_targets * 2, sizeof(pthread_mutex_t *));
    for (int i = 0; i < cov->n_targets; ++i) {
        int32_t block_count =  ((cov->target_len[i]-1)/cov->coverage_block_size)+1;
        int32_t block_mutex_count = ((block_count-1u)>>cov->coverage_mutex_shift)+1;
        cov->coverage_block_mutexes[i] =  calloc(block_mutex_count, sizeof(pthread_mutex_t));
        for (int j = 0; j < block_mutex_count; ++j) pthread_mutex_init(&cov->coverage_block_mutexes[i][j], NULL);
    }
    for (int i = cov->n_targets; i < cov->n_targets * 2; ++i){
        int32_t block_count =  ((cov->target_len[i - cov->n_targets]-1)/cov->coverage_block_size)+1;
        int32_t block_mutex_count = ((block_count-1u)>>cov->coverage_mutex_shift)+1;
        cov->coverage_block_mutexes[i] =  calloc(block_mutex_count, sizeof(pthread_mutex_t));
        for (int j = 0; j < block_mutex_count; ++j) pthread_mutex_init(&cov->coverage_block_mutexes[i][j], NULL);
    }
    return cov;
}

int coverage2_destroy(coverage2_t * cov){
    for (int i = 0; i < cov->n_targets; ++i) {
        uint32_t block_count = ((cov->target_len[i]-1)>>cov->coverage_block_shift)+1;
        for (int j = 0; j < block_count; ++j) if (cov->coverage_blocks[i][j]!=NULL) free(cov->coverage_blocks[i][j]);
        free(cov->coverage_blocks[i]);
        if (cov->is_mt) {
            uint32_t block_mutex_count = ((block_count-1)>>cov->coverage_mutex_shift)+1;
            for (int j = 0; j < block_mutex_count; ++j) pthread_mutex_destroy(&cov->coverage_block_mutexes[i][j]);
            free(cov->coverage_block_mutexes[i]);
        }
    }
    for (int i = cov->n_targets; i < cov->n_targets * 2; ++i){
        uint32_t block_count = ((cov->target_len[i - cov->n_targets]-1)>>cov->coverage_block_shift)+1;
        for (int j = 0; j < block_count; ++j) if (cov->coverage_blocks[i][j]!=NULL) free(cov->coverage_blocks[i][j]);
        free(cov->coverage_blocks[i]);
        if (cov->is_mt) {
            uint32_t block_mutex_count = ((block_count-1)>>cov->coverage_mutex_shift)+1;
            for (int j = 0; j < block_mutex_count; ++j) pthread_mutex_destroy(&cov->coverage_block_mutexes[i][j]);
            free(cov->coverage_block_mutexes[i]);
        }
    }
    free(cov->coverage_blocks);
    if (cov->is_mt) free(cov->coverage_block_mutexes);
    for (int i=0; i < cov->n_targets; ++i) free(cov->target_name[i]);
    free(cov->target_name);
    free(cov->target_len);
    free(cov);
    return 0;
}

/* starts are 0-based and ends are 1-based */
static int base2index[16] = {-1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, 4};
int coverage2_update(coverage2_t *cov, int32_t target, uint32_t start, uint32_t end, char strand, uint8_t *read, int read_pos){
    /* this function simply trust its arguments without checking whether target is present and coordinate is in valid range */
    uint32_t block_index_start, block_index_end, block_index, block_start, block_end, new_start, new_end;
    cov2_val_t *coverage_block;
    cov2_val_t **coverage_block_target;
    pthread_mutex_t *mutex;
    int32_t target_index = target + (strand == '-' ? cov->n_targets: 0);
    block_index_start = start/cov->coverage_block_size;
    block_index_end = (end-1)/cov->coverage_block_size;
    while (block_index_start <= block_index_end){
        block_index = block_index_start;
        block_start = block_index<<cov->coverage_block_shift;
        new_start = (start > block_start) ? start - block_start : 0;
        new_end = (end - block_start > cov->coverage_block_size) ? cov->coverage_block_size : end - block_start;
        coverage_block_target = cov->coverage_blocks[target_index];
        if (cov->is_mt) {
            mutex=&cov->coverage_block_mutexes[target_index][block_index>>cov->coverage_mutex_shift];
            pthread_mutex_lock(mutex);
        }
        coverage_block = coverage_block_target[block_index];
        if (coverage_block == NULL) {
            uint32_t target_len = cov->target_len[target];
            int needed=cov->coverage_block_size > target_len - block_start ? target_len - block_start : cov->coverage_block_size ;
            coverage_block = calloc(needed, sizeof(cov2_val_t));
            cov->coverage_blocks[target_index][block_index] = coverage_block;
        }
        while (new_start < new_end) {
            (*(((double *)(&coverage_block[new_start++]))+base2index[bam_seqi(read, read_pos)]))++;
            read_pos++;
        }
        if (cov->is_mt) pthread_mutex_unlock(mutex);
        block_index_start++;
    }
    return 0;
}

typedef struct interval_s{
    char *target;
    uint32_t target_len;
    int capacity;
    int size;
    uint32_t *start;
    uint32_t *end;
    float *value;
} interval_t;

static void *interval_init() {
    interval_t *i = calloc(1, sizeof(interval_t));
    return i;
}

static int interval_destroy(interval_t *i) {
    free(i->start);
    free(i->end);
    free(i->value);
    free(i);
    return 0;
}

static int interval_expand(interval_t *i) {
    int new_capacity = i->size<<1;
    if (new_capacity < 1000) new_capacity = 1000;
    uint32_t *new_start = realloc(i->start, new_capacity * sizeof(uint32_t));
    if (!new_start) return -1;
    uint32_t *new_end = realloc(i->end, new_capacity * sizeof(uint32_t));
    if (!new_end) {
        free(new_start);
        return -1;
    }
    float *new_value= realloc(i->value, new_capacity * sizeof(float));
    if (!new_value){
        free(new_start);
        free(new_end);
        return -1;
    }
    i->capacity = new_capacity;
    i->start = new_start;
    i->end = new_end;
    i->value = new_value;
    return 0;
}

#define interval_add(i, start, end, value) {\
int index = i->size++; \
if (index == i->capacity)  interval_expand(i); \
i->start[index] = start; \
i->end[index] = end; \
i->value[index] = value; \
}

struct extract_interval_arg{
    cov_val_t **coverage_blocks;
    uint32_t block_index_start;
    uint32_t block_index_end;
    uint32_t bin_size;
    uint32_t block_size;
    interval_t *itv;
};
struct extract_interval_arg *extract_interval_arg_init(){
    struct extract_interval_arg *arg;
    arg = malloc(sizeof(struct extract_interval_arg));
    if (!arg) return NULL;
    arg->itv = interval_init();
    if (!arg->itv) {
        free(arg);
        return NULL;
    }
    return arg;
}
void extract_interval_arg_destroy(void *_arg){
    struct extract_interval_arg *arg = _arg;
    interval_destroy(arg->itv);
    free(arg);
}
void *extract_interval(void *_args){
    struct extract_interval_arg *args = _args;
    interval_t *itv = args->itv;
    cov_val_t ** coverage_blocks = args->coverage_blocks;
    uint32_t block_index_start = args->block_index_start;
    uint32_t block_index_end = args->block_index_end;
    uint32_t bin_size = args->bin_size;
    uint32_t block_size = args->block_size;

    uint32_t target_len = itv->target_len;
    int32_t bin_count = (target_len-1)/bin_size+1;
    int32_t block_count =  ((bin_count-1)/block_size)+1;

    uint32_t block_bin_count;
    uint32_t block_index;

    block_index = block_index_start;
    cov_val_t *coverage = coverage_blocks[block_index];
    if (block_index == block_count -1) block_bin_count = bin_count - (block_index) * block_size;
    else block_bin_count = block_size;


    int bin_index = 0;

    int start =  bin_size * (block_index * block_size);
    int end = 0;

    uint32_t range_end = block_index_end * block_size * bin_size;
    if (range_end > target_len) range_end = target_len;

    while (start < range_end){
        cov_val_t value = coverage?coverage[bin_index]:0;
        while (1){
            if (!coverage){
                if (value == 0){
                    block_index++;
                    bin_index=0;
                    if (block_index < block_index_end) {
                        coverage = coverage_blocks[block_index];
                        if (block_index == block_count -1) block_bin_count = bin_count - (block_index) * block_size;
                        else block_bin_count = block_size;
                    } else break;
                } else break;
            } else {
                if (value == coverage[bin_index]){
                    bin_index++;
                    if (bin_index == block_bin_count){
                        if (block_index == block_count -1) break;
                        block_index++;
                        bin_index = 0;
                        if (block_index < block_index_end) {
                            coverage = coverage_blocks[block_index];
                            if (block_index == block_count -1) block_bin_count = bin_count - (block_index) * block_size;
                            else block_bin_count = block_size;
                        } else break;
                    }
                } else break;
            }
        }

        end = (block_index * block_size + bin_index) * bin_size;
        interval_add(itv, start, end, value);
        start = end;
    }
    if (itv->end[itv->size-1] > itv->target_len) itv->end[itv->size-1] = itv->target_len;
    return _args;
}

struct output_bw_mt_writer_arg{
    mt_queue *q;
    mt_buffer *b;
    bigWigFile_t *fp;
};
void * output_bw_mt_writer(void *_arg){
    mt_queue *q = ((struct output_bw_mt_writer_arg *) _arg)->q;
    mt_buffer *b = ((struct output_bw_mt_writer_arg *) _arg)->b;
    bigWigFile_t *fp = ((struct output_bw_mt_writer_arg *) _arg)->fp;

    struct extract_interval_arg *arg;
    interval_t *itv;
    uint32_t itv_last_start;
    uint32_t itv_last_end;
    float itv_last_value;
    int block_count = 0;
    int no_last = 1;
    int init = 1;
    char *last_target = "\0";
    while (mt_queue_receive(q, (void *)&arg, 0) == 0){
        time_t t;
        itv = arg->itv;
        /* check if new target is meet */
        if (strcmp(last_target, itv->target)!=0){
            init = 1;
            no_last = 1;
            block_count =  (((itv->target_len-1)/arg->bin_size-1)/arg->block_size)+1;
            last_target = itv->target;
        }
        if (itv->size > 0) {
            /* handle last issues */
            if (!no_last) {
                if (itv_last_value == itv->value[0] && itv_last_end == itv->start[0])
                    itv->start[0] = itv_last_start;
                else {
                    if (init) {
                        init = 0;
                        bwAddIntervals(fp, &itv->target, &itv_last_start, &itv_last_end, &itv_last_value, 1);
                    } else bwAppendIntervals(fp, &itv_last_start, &itv_last_end, &itv_last_value, 1);
                }
            }

            if (arg->block_index_end < block_count) {
                itv_last_start = itv->start[itv->size - 1];
                itv_last_end = itv->end[itv->size - 1];
                itv_last_value = itv->value[itv->size - 1];
                no_last = 0;
                itv->size--;
            }
        }

        if (itv->size > 0) {
            if (init) {
                init = 0;
                bwAddIntervals(fp, &itv->target, &itv->start[0], &itv->end[0], &itv->value[0], 1);
                bwAppendIntervals(fp, itv->start + 1, itv->end + 1, itv->value + 1, itv->size - 1);
            } else bwAppendIntervals(fp, itv->start, itv->end, itv->value, itv->size);
        }

        itv->size = 0;
        mt_buffer_put(b, arg);
    }
    return NULL;
}

int output_bw(coverage_t *cov, char *fn, mt_server *s){
    /* some necessary preparation */
    bigWigFile_t *fp = NULL;
    if(bwInit(1u<<17u) != 0) return 1;
    fp = bwOpen(fn, NULL, "w");
    if(!fp) return 1;
    if(bwCreateHdr(fp, 10)) return 1;
    fp->cl = bwCreateChromList(cov->target_name, cov->target_len, cov->n_targets);
    if(!fp->cl) return 1;
    if(bwWriteHdr(fp)) return 1;
    if (s) bwMtInit(fp, s);

    /* malloc buffer*/
    uint32_t itv_last_start;
    uint32_t itv_last_end;
    float itv_last_value;
    mt_queue *q;
    mt_buffer *b;
    interval_t *itv;
    struct output_bw_mt_writer_arg mt_writer_arg;
    pthread_t mt_writer;
    struct extract_interval_arg *arg;
    if (!s) arg = extract_interval_arg_init();
    else {
        int n_thread = mt_server_n_thread(s);
        q = mt_queue_init(s, INT_MAX, INT_MAX, MT_QUEUE_MODE_SERIAL);
        b = mt_buffer_init();
        for (int i = 0; i < n_thread * 2; ++i) mt_buffer_put(b, extract_interval_arg_init());
        mt_writer_arg.q = q;
        mt_writer_arg.b = b;
        mt_writer_arg.fp = fp;
        pthread_create(&mt_writer, NULL, output_bw_mt_writer, &mt_writer_arg);
    }
    /* start writing */
    for (int i = 0; i < cov->n_targets; ++i){
        int init = 1;
        int no_last = 1;
        int block_count = ((cov->target_len[i]-1)/cov->coverage_block_size)+1;
        int block_index_start = 0, block_index_end = 0;
        int n_needed_block = (1u<<17u)/cov->coverage_block_size+1;
        int n_block = 0;
        while (block_index_end < block_count){
            block_index_start = block_index_end;

            n_block=0;
            while (block_index_end < block_count) {
                if (!cov->coverage_blocks[i][block_index_end++]) continue;
                if (n_block == n_needed_block) break;
                n_block++;
            }

            if (!s) {
                arg->block_index_start = block_index_start;
                arg->block_index_end = block_index_end;
                arg->bin_size = cov->bin_size;
                arg->block_size = cov->coverage_block_size;
                arg->coverage_blocks = cov->coverage_blocks[i];
                arg->itv->target = cov->target_name[i];
                arg->itv->target_len = cov->target_len[i];
                arg = extract_interval(arg);
                itv = arg->itv;

                if (itv->size > 0) {
                    /* handle last issues */
                    if (!no_last) {
                        if (itv_last_value == itv->value[0] && itv_last_end == itv->start[0])
                            itv->start[0] = itv_last_start;
                        else {
                            if (init) {
                                init = 0;
                                bwAddIntervals(fp, &itv->target, &itv_last_start, &itv_last_end, &itv_last_value, 1);
                            } else bwAppendIntervals(fp, &itv_last_start, &itv_last_end, &itv_last_value, 1);
                        }
                    }

                    if (block_index_end < block_count) {
                        itv_last_start = itv->start[itv->size - 1];
                        itv_last_end = itv->end[itv->size - 1];
                        itv_last_value = itv->value[itv->size - 1];
                        no_last = 0;
                        itv->size--;
                    }
                }

                if (itv->size > 0) {
                    if (init) {
                        init = 0;
                        bwAddIntervals(fp, &itv->target, &itv->start[0], &itv->end[0], &itv->value[0], 1);
                        bwAppendIntervals(fp, itv->start + 1, itv->end + 1, itv->value + 1, itv->size - 1);
                    } else bwAppendIntervals(fp, itv->start, itv->end, itv->value, itv->size);
                }

                itv->size = 0;
            } else {
                arg = mt_buffer_get(b);
                arg->block_index_start = block_index_start;
                arg->block_index_end = block_index_end;
                arg->bin_size = cov->bin_size;
                arg->block_size = cov->coverage_block_size;
                arg->coverage_blocks = cov->coverage_blocks[i];
                arg->itv->target = cov->target_name[i];
                arg->itv->target_len = cov->target_len[i];
                mt_queue_dispatch(q, extract_interval, arg, NULL, NULL, 0);
            }

        }
    }
    if (!s) {
        extract_interval_arg_destroy(arg);
    } else {
        mt_queue_dispatch_end(q);
        mt_queue_wait(q, MT_FINISH);
        pthread_join(mt_writer, NULL);
        mt_queue_destroy(q);
        mt_buffer_destroy(b, &extract_interval_arg_destroy);
    }
    bwClose(fp);
    bwCleanup();
    return 0;
}