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

#include "bam.h"
#include "mt.h"
#include "fa.h"

#include "common.h"
#include "coverage.h"


#define FR_UNSTRANDED 0
#define FR_FIRSTSTRAND 1
#define FR_SECONDSTRAND 2

KHASH_MAP_INIT_STR(target, int);

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

static struct {
    char *fn;
    char *fa;
    char *fai;
    char *out;
    char *bed;
    double count;
    double prop;
    int library_type;
    int n_threads;
} parameter;

static void parse_arg(int argc, char *argv[]);
static void usage(char *msg);
#define is_reverse(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define is_mate_reverse(b) (((b)->core.flag & BAM_FMREVERSE) != 0)
#define is_first(b) (((b)->core.flag & BAM_FREAD1) !=0u)
#define is_second(b) (((b)->core.flag & BAM_FREAD2) !=0u)

static char get_strand(bam1_t *b, int type) {
    int is_paired = 0;
    if ((((b)->core.flag & BAM_FPAIRED) != 0u)) is_paired = 1;
    if (type == FR_FIRSTSTRAND) {
        if (is_paired) {
            if ((is_first(b) && is_reverse(b)) || (is_second(b) && !is_reverse(b)))
                return '+';
            else if ((is_first(b) && !is_reverse(b)) || (is_second(b) && is_reverse(b)))
                return '-';
        } else return ((is_reverse(b) ? '+' : '-'));
    } else if (type == FR_SECONDSTRAND) {
        if (is_paired) {
            if ((is_first(b) && !is_reverse(b)) || (is_second(b) && is_reverse(b)))
                return '+';
            else if ((is_first(b) && is_reverse(b)) || (is_second(b) && !is_reverse(b)))
                return '-';
        } else return ((is_reverse(b) ? '-' : '+'));
    } else if (type == FR_UNSTRANDED) return '.';
    fprintf(stderr, "%d\t%d\n", type, b->core.flag);
    exit(1);
}

int extract_mutation(bam1_t *b, coverage2_t *cov){
    char strand = get_strand(b, parameter.library_type);
    int pos=b->core.pos, read_pos = 0;
    const uint32_t *cigar=bam_get_cigar(b);
    uint8_t *read_seq = bam_get_seq(b);
    int cigar_len=0;
    int cigar_type=0;
    for (int i=0; i<b->core.n_cigar; ++i){
        cigar_len=bam_cigar_oplen(cigar[i]);
        cigar_type=bam_cigar_type(bam_cigar_op(cigar[i]));
        if (cigar_type==1) read_pos+=cigar_len;
        if (cigar_type==2) pos+=cigar_len;
        if (cigar_type==3) {
            coverage2_update(cov, b->core.tid, pos, pos+cigar_len, strand == '-'?'-':'+', read_seq, read_pos);
            pos+=cigar_len;
            read_pos+=cigar_len;
        }
    }
    return 0;
}

struct call_mutation_arg{
    coverage2_t *cov;
    int index_start;
    int index_end;
    int32_t target_index;
    int32_t target;
    char strand;
    fa_t *fa;
    char *base2int;
    FILE *out;
};

void *call_mutation(void *_args){
    struct call_mutation_arg *args = _args;
    coverage2_t *cov =  args->cov;
    int index_start = args->index_start;
    int index_end = args->index_end;
    int32_t target_index = args->target_index;
    int32_t target = args->target;
    uint32_t coverage_block_size = cov->coverage_block_size;
    int32_t target_len = cov->target_len[target];
    char *seq = NULL;
    int seq_len = 0;
    for (int block_index = index_start; block_index < index_end; block_index++){
        if (cov->coverage_blocks[target_index][block_index] == NULL) continue;
        else {
            int32_t coverage_block_start = block_index * coverage_block_size;
            int32_t coverage_block_end = (block_index + 1) * coverage_block_size;
            if (args->fa) {
                seq = extract_sequence(args->fa, cov->target_name[target], coverage_block_start, coverage_block_end, '+', NULL, &seq_len);
            }
            if (coverage_block_end > target_len) coverage_block_end = target_len;
            for (int i = coverage_block_start; i < coverage_block_end; ++i) {
                double *counts = &(cov->coverage_blocks[target_index][block_index][i - coverage_block_start]);
                double count_sum = counts[0] + counts[1] + counts[2] + counts[3] + counts[4];
                double count_ref = 0;
                if (count_sum < parameter.count) continue;
                char base = '?';
                if (seq) base = seq[i - coverage_block_start];
                if (base != '?'){
                    count_ref = counts[args->base2int[base]];
                } else {
                    count_ref = counts[0];
                    for (int j = 1; j < 5; ++j) if (count_ref < counts[j]) count_ref = counts[j];
                }
                if ((1 - count_ref/count_sum) < parameter.prop) continue;
                fprintf(args->out, "%s\t%d\t%c\t%c\t%f\t%f\t%f\t%f\t%f\n", cov->target_name[target] ,i+1, args->strand, base, counts[0], counts[1], counts[2], counts[3], counts[4]);
            }
            }
        }
    if (seq) free(seq);
    return NULL;
}

struct {
    double A;
    double C;
    double G;
    double T;
    double N;
} test_mutation_null = {0, 0, 0, 0, 0};

struct test_mutation_arg{
    coverage2_t *cov;
    int32_t target;
    int32_t chromStart;
    int32_t chromEnd;
    char strand;
    fa_t *fa;
    char *base2int;
    FILE *out;
};

void *test_mutation(void *_args){
    struct test_mutation_arg *args = _args;
    coverage2_t *cov = args->cov;
    int32_t chromStart = args->chromStart;
    int32_t chromEnd = args->chromEnd;
    int32_t target = args->target;
    int32_t target_index = (args->strand == '-') ? target + cov->n_targets : target;
    int block_index_start = chromStart / cov->coverage_block_size;
    int block_index_end = (chromEnd - 1) / cov->coverage_block_size + 1;
    for (int block_index = block_index_start; block_index < block_index_end; ++block_index){
        int32_t coverage_block_start = block_index * cov->coverage_block_size;
        int32_t coverage_block_end = (block_index + 1) * cov->coverage_block_size;
        int32_t start = chromStart > coverage_block_start ? chromStart : coverage_block_start;
        int32_t end = chromEnd > coverage_block_end ? coverage_block_end : chromEnd;
        for (int i = start; i < end; ++i){
            double *counts;
            if (cov->coverage_blocks[target_index][block_index] == NULL) counts = &test_mutation_null;
            else counts = &(cov->coverage_blocks[target_index][block_index][i - coverage_block_start]);
            fprintf(args->out, "%s\t%d\t%c\t%c\t%f\t%f\t%f\t%f\t%f\n", cov->target_name[target] ,i+1, args->strand, '?', counts[0], counts[1], counts[2], counts[3], counts[4]);
        }
    }
    return NULL;
}

int samvt_mutation(int argc, char *argv[]){
    parse_arg(argc, argv);
    bt_bam_t *s = bt_bam_open(parameter.fn, parameter.n_threads);
    fa_t *fa = NULL;
    if (parameter.fa) fa = fa_open(parameter.fa, parameter.fai);
    coverage2_t *cov = coverage2_init(s->hdr->n_targets, s->hdr->target_name, s->hdr->target_len, 12);
    char base2int[256];
    for (int i = 0; i < 256; ++i) base2int[i] = 4;
    base2int['a'] = 0;
    base2int['A'] = 0;
    base2int['c'] = 1;
    base2int['C'] = 1;
    base2int['g'] = 2;
    base2int['G'] = 2;
    base2int['t'] = 3;
    base2int['T'] = 3;
    if (parameter.n_threads == 0){
        bam1_t *b1 = bam_init1();
        while (bt_bam_next(s, b1) == 0) extract_mutation(b1, cov);
        bam_destroy1(b1);
    }
    FILE *out = fopen(parameter.out, "w");
    if (!parameter.bed) {
        for (int i = 0; i < cov->n_targets * 2; ++i) {
            int target_index = i;
            int target = target_index >= cov->n_targets ? target_index - cov->n_targets : target_index;
            int block_count = ((cov->target_len[target] - 1) / cov->coverage_block_size) + 1;
            int block_index_start = 0, block_index_end = 0;
            int n_needed_block = (1u << 17u) / cov->coverage_block_size + 1;
            int n_block = 0;
            while (block_index_end < block_count) {
                block_index_start = block_index_end;
                n_block = 0;
                while (block_index_end < block_count) {
                    if (!cov->coverage_blocks[i][block_index_end++]) continue;
                    if (n_block == n_needed_block) break;
                    n_block++;
                }
                struct call_mutation_arg args;
                args.cov = cov;
                args.index_start = block_index_start;
                args.index_end = block_index_end;
                args.strand = (target_index == target) ? '+' : '-';
                args.target_index = target_index;
                args.target = target;
                args.fa = fa;
                args.out = out;
                args.base2int = base2int;
                call_mutation(&args);
            }

        }
    } else {
        kh_target_t *chrom2id = kh_init(target);
        for (int i  = 0; i <  cov->n_targets; ++i){
            int kh_ret;
            khiter_t key = kh_put(target, chrom2id, cov->target_name[i], &kh_ret);
            kh_val(chrom2id, key) = i;
        }
        char line[4096];
        char *item[7];
        FILE *bed = fopen(parameter.bed, "r");
        while (fgets(line, 4096, bed)){
            strsplit(line, item, 7, '\t');
            struct test_mutation_arg args;
            khiter_t key = kh_get(target, chrom2id, item[0]);
            if (key == kh_end(chrom2id)) continue;
            args.cov = cov;
            args.target = kh_val(chrom2id, key);
            args.chromStart = strtol(item[1], NULL, 10);
            args.chromEnd = strtol(item[2], NULL, 10);
            args.strand = item[5][0];
            args.out = out;
            test_mutation(&args);
        }
        kh_destroy(target, chrom2id);
        fclose(bed);
    }
    fclose(out);
    coverage2_destroy(cov);
    bt_bam_close(s);
    if (fa) fa_close(fa);
    return 0;
}

static void parse_arg(int argc, char *argv[]){
    char c;
    int show_help=0;

    parameter.fn = NULL;
    parameter.fa = NULL;
    parameter.fai = NULL;
    parameter.out = NULL;
    parameter.prop = 0.15;
    parameter.count = 50;
    parameter.library_type = FR_UNSTRANDED;
    parameter.n_threads = 0;


    if (argc == 1) usage("");
    const char *shortOptions = "ho:i:f:p:t:a:b:c:e:";
    const struct option longOptions[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "out" , required_argument , NULL, 'o' },
                    { "bam" , required_argument, NULL, 'i' },
                    { "bed" , required_argument, NULL, 'b' },
                    { "fa" , required_argument, NULL, 'a' },
                    { "library-type" , required_argument, NULL, 't' },
                    { "count" , required_argument, NULL, 'c' },
                    { "prop" , required_argument, NULL, 'e' },
                    { "threads" , required_argument, NULL, 'p' },
                    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
            };

    while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'h':
                show_help = 1;
                break;
            case 'o':
                parameter.out = optarg;
                break;
            case 'i':
                parameter.fn=optarg;
                break;
            case 'a':
                parameter.fa=optarg;
                break;
            case 'b':
                parameter.bed=optarg;
                break;
            case 'c':
                parameter.count=strtod(optarg, NULL);
                break;
            case 'e':
                parameter.prop=strtod(optarg, NULL);
                break;
            case 't':
                if (strcmp(optarg, "fr-firststrand") == 0) parameter.library_type = FR_FIRSTSTRAND;
                else if (strcmp(optarg, "fr-secondstrand") == 0) parameter.library_type = FR_SECONDSTRAND;
                else usage("Unknown value for -t/--library-type.");
                break;
            case 'p':
                parameter.n_threads = strtol(optarg, NULL, 10);
                break;
            default:
                usage("Unknown parameter.");
        }
    }
    if (parameter.fn == NULL) parameter.fn = "/dev/stdin";
    if (parameter.out == NULL) parameter.out = "/dev/stdout";
    if (parameter.fa) {
        parameter.fai = strcpy(malloc(strlen(parameter.fa)+5), parameter.fa);
        strcpy(parameter.fai + strlen(parameter.fa), ".fai");
    }
    if (argc != optind) usage("Unrecognized parameter");
    if (show_help)    usage("");
}

static void usage(char *msg)
{
    const char *usage_info = "Usage:  samvt mutation [options] --bam <alignment file> --bw <big wig file>\n \
[options]\n\
-i/--bam                       : bam alignment file. [required]\n\
-o/--out                       : tsv file for output. [required]\n\
-a/--fa                        : use the reference fasta file to determine variant bases\n\
-h/--help                      : show help informations.\n\
-t/--library-type              : library type, one of fr-firststrand or fr-secondstrand, default: fr-firststrand.\n\
-b/--bed                       : exclude the position not specified by bed file.\n\
-c/--count                     : exclude the position with lower coverage.\n\
-e/--prop                      : exclude the position with lower proportion of variants.\n\
-p/--threads                   : number of threads to use (not implemented). \n\n\
";
    if (msg==NULL || msg[0] == '\0') fprintf(stderr, "%s", usage_info);
    else fprintf(stderr, "%s\n\n%s", msg, usage_info);
    exit(1);
}