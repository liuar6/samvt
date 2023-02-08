/// @file htslib/bgzf.h
/// Low-level routines for direct BGZF operations.
/*
   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology
                 2011, 2012 Attractive Chaos <attractor@live.co.uk>
   Copyright (C) 2009, 2013, 2014, 2017, 2018-2019 Genome Research Ltd

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

/* The BGZF library was originally written by Bob Handsaker from the Broad
 * Institute. It was later improved by the SAMtools developers. */

/* This file is an extract from the original "bgzf.h" from htslib. */

typedef struct pool_alloc_t pool_alloc_t;
typedef struct bgzf_job bgzf_job;
typedef struct hts_tpool hts_tpool;
typedef struct hts_tpool_process hts_tpool_process;

typedef struct {
    hts_pos_t beg, end;
    int tid, is_mapped;  // args for hts_idx_push
    uint64_t offset, block_number;
} hts_idx_cache_entry;

typedef struct {
    int nentries, mentries; // used and allocated
    hts_idx_cache_entry *e; // hts_idx elements
} hts_idx_cache_t;

enum mtaux_cmd {
    NONE = 0,
    SEEK,
    SEEK_DONE,
    HAS_EOF,
    HAS_EOF_DONE,
    CLOSE,
};
typedef struct bgzf_mtaux_t {
    // Memory pool for bgzf_job structs, to avoid many malloc/free
    pool_alloc_t *job_pool;
    bgzf_job *curr_job;

    // Thread pool
    int n_threads;
    int own_pool;
    hts_tpool *pool;

    // Output queue holding completed bgzf_jobs
    hts_tpool_process *out_queue;

    // I/O thread.
    pthread_t io_task;
    pthread_mutex_t job_pool_m;
    int jobs_pending; // number of jobs waiting
    int flush_pending;
    void *free_block;
    int hit_eof;  // r/w entirely within main thread

    // Message passing to the reader thread; eg seek requests
    int errcode;
    uint64_t block_address;
    int eof;
    pthread_mutex_t command_m; // Set whenever fp is being updated
    pthread_cond_t command_c;
    enum mtaux_cmd command;

    // For multi-threaded on-the-fly indexing. See bgzf_idx_push below.
    pthread_mutex_t idx_m;
    hts_idx_t *hts_idx;
    uint64_t block_number, block_written;
    hts_idx_cache_t idx_cache;
} mtaux_t;
