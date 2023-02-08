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

#include "bigWig.h"

#include "bam.h"
#include "mt.h"
#include "mt_buffer.h"

#include "coverage.h"

#define FR_UNSTRANDED 0
#define FR_FIRSTSTRAND 1
#define FR_SECONDSTRAND 2

#define STRAND_ALL 0
#define STRAND_FORWARD 1
#define STRAND_REVERSE 2

#define SELECT_ALL 0
#define SELECT_FIRST_FORWARD 1
#define SELECT_FIRST_REVERSE 2


static struct {
    char *fn;
    char *out;
    int bin_size;
    int library_type;
    int strand;
    int n_threads;

    int select;
} parameter;

typedef struct samvt_coverage_job_s{
    bam1_t **bam;
    int size;
    int capacity;
    mt_buffer *bf;
    coverage_t *cov;
} samvt_coverage_job_t;

samvt_coverage_job_t *samvt_coverage_job_init(int capacity){
    samvt_coverage_job_t *job;
    job = malloc(sizeof(*job));
    job->size = 0;
    job->capacity = capacity;
    job->bam = malloc(sizeof(*job->bam) * capacity);
    for (int i = 0; i < job->capacity; ++i) job->bam[i] = bam_init1();
    return job;
}

void samvt_coverage_job_destroy(void *_job){
    samvt_coverage_job_t *job = _job;
    for (int i = 0; i < job->capacity; ++i) bam_destroy1(job->bam[i]);
    free(job->bam);
    free(job);
}


int extract_coverage(bam1_t *b, coverage_t *cov){
    int select = parameter.select;
    if (select != SELECT_ALL){
        uint16_t flag = b->core.flag;
        if ((flag & BAM_FPAIRED)) {
            if (select == SELECT_FIRST_REVERSE){
                if (((flag & BAM_FREAD1) && !(flag & BAM_FREVERSE)) || ((flag & BAM_FREAD2) && (flag & BAM_FREVERSE))) return 0;
            } else if (select == SELECT_FIRST_FORWARD){
                if (((flag & BAM_FREAD1) && (flag & BAM_FREVERSE)) || ((flag & BAM_FREAD2) && !(flag & BAM_FREVERSE))) return 0;
            }
        } else if ((select == SELECT_FIRST_REVERSE && !(flag & BAM_FREVERSE)) || (select == SELECT_FIRST_FORWARD && (flag & BAM_FREVERSE))) return 0;
    }
    int pos=b->core.pos;
    const uint32_t *cigar=bam_get_cigar(b);
    int cigar_len=0;
    int cigar_type=0;
    for (int i=0; i<b->core.n_cigar; ++i){
        cigar_len=bam_cigar_oplen(cigar[i]);
        cigar_type=bam_cigar_type(bam_cigar_op(cigar[i]));
        if (cigar_type==2) pos+=cigar_len;
        if (cigar_type==3) {
            coverage_update(cov, b->core.tid, pos, pos+cigar_len);
            pos+=cigar_len;
        }
    }
    return 0;
}

void *extract_coverage_mt(void *arg){
    samvt_coverage_job_t* j = arg;
    for (int i = 0; i < j->size; ++i){
        extract_coverage(j->bam[i], j->cov);
    }
    mt_buffer_put(j->bf, j);
    return NULL;
}

static void parse_arg(int argc, char *argv[]);
static void usage(char *msg);

int samvt_coverage(int argc, char *argv[]){
    parse_arg(argc, argv);
    int select = SELECT_ALL;
    if ((parameter.library_type == FR_FIRSTSTRAND && parameter.strand == STRAND_FORWARD) || (parameter.library_type == FR_SECONDSTRAND && parameter.strand == STRAND_REVERSE)) select = SELECT_FIRST_REVERSE;
    if ((parameter.library_type == FR_FIRSTSTRAND && parameter.strand == STRAND_REVERSE) || (parameter.library_type == FR_SECONDSTRAND && parameter.strand == STRAND_FORWARD)) select = SELECT_FIRST_FORWARD;
    parameter.select=select;
    bt_bam_t *s = bt_bam_open(parameter.fn, parameter.n_threads);
    coverage_t *cov = coverage_init(s->hdr->n_targets, s->hdr->target_name, s->hdr->target_len, 12);
    if (parameter.n_threads == 0){
        bam1_t *b1 = bam_init1();
        while (bt_bam_next(s, b1) == 0) extract_coverage(b1, cov);
        bam_destroy1(b1);
    } else {
        mt_queue *q = mt_queue_init(bt_bam_mt_server(s), parameter.n_threads * 8, 0, MT_QUEUE_MODE_IGNORED);
        mt_buffer *bf = mt_buffer_init();
        for (int i = 0; i < parameter.n_threads * 5; ++i) mt_buffer_put(bf, samvt_coverage_job_init(10000));
        coverage_mt(cov);
        while(1){
            int ret1;
            samvt_coverage_job_t *job = mt_buffer_get(bf);
            job->size = 0;
            job->cov = cov;
            job->bf = bf;
            while(job->size < job->capacity && (ret1=bt_bam_next(s, job->bam[job->size]))==0) ++job->size;
            mt_queue_dispatch(q, extract_coverage_mt, job, NULL, NULL, 0);
            if (ret1 != 0) {
                mt_queue_dispatch_end(q);
                break;
            }
        }
        mt_queue_wait(q, MT_FINISH);
        mt_queue_destroy(q);
        mt_buffer_destroy(bf, &samvt_coverage_job_destroy);
    }
    output_bw(cov, parameter.out, parameter.n_threads?bt_bam_mt_server(s):NULL);
    bt_bam_close(s);
    coverage_destroy(cov);
    return 0;
}

static void parse_arg(int argc, char *argv[]){
    char c;
    int show_help=0;

    parameter.fn = NULL;
    parameter.out = NULL;
    parameter.bin_size = 1;
    parameter.library_type = FR_FIRSTSTRAND;
    parameter.strand = STRAND_ALL;
    parameter.n_threads = 0;


    if (argc == 1) usage("");
    const char *shortOptions = "ho:i:t:s:B:I:p:";
    const struct option longOptions[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "bw" , required_argument , NULL, 'o' },
                    { "bam" , required_argument, NULL, 'i' },
                    { "library-type" , required_argument, NULL, 't' },
                    { "strand" , required_argument, NULL, 's' },
                    { "bin-size" , required_argument, NULL, 'B' },
                    { "item-size" , required_argument, NULL, 'I' },
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
            case 't':
                if (strcmp(optarg, "fr-firststrand") == 0) parameter.library_type = FR_FIRSTSTRAND;
                else if (strcmp(optarg, "fr-secondstrand") == 0) parameter.library_type = FR_SECONDSTRAND;
                else usage("Unknown value for -t/--library-type.");
                break;
            case 's':
                if (strcmp(optarg, "forward") == 0) parameter.strand = STRAND_FORWARD;
                else if (strcmp(optarg, "reverse") == 0) parameter.strand = STRAND_REVERSE;
                else usage("Unknown value for -s/--strand.");
                break;
            case 'B':
                parameter.bin_size = strtol(optarg, NULL, 10);
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
    if (argc != optind) usage("Unrecognized parameter");
    if (show_help)    usage("");
}

static void usage(char *msg)
{
    const char *usage_info = "Usage:  samvt coverage [options] --bam <alignment file> --bw <big wig file>\n \
[options]\n\
-i/--bam                       : bam alignment file. [required]\n\
-o/--bw                        : bigwig file for output. [required]\n\
-h/--help                      : show help informations.\n\
-t/--library-type              : library type, one of fr-firststrand or fr-secondstrand, default: fr-firststrand.\n\
-s/--strand                    : strand on the genome used for coverage calculation.\n\
-B/--bin-size                  : bin size for coverage calculation (not implemented). \n\
-p/--threads                   : number of threads to use. \n\n\
";
    if (msg==NULL || msg[0] == '\0') fprintf(stderr, "%s", usage_info);
    else fprintf(stderr, "%s\n\n%s", msg, usage_info);
    exit(1);
}
