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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define SAMVT_VERSION "1.0.0"

int samvt_coverage(int argc, char *argv[]);
int samvt_mutation(int argc, char *argv[]);

static void usage(char *msg){
    const char *usage_info="\
samvt: versatile tools for jobs associated with bam file.\n\
usage:    samvt <subcommand> [options]\n\n";
    if (msg==NULL || msg[0] == '\0') fprintf(stderr, "%s", usage_info);
    else fprintf(stderr, "%s\n\n%s", msg, usage_info);
    exit(1);
}

static void version(){
    fprintf(stderr, "samvt-%s\n\n", SAMVT_VERSION);
    exit(1);
}


int main(int argc, char *argv[]) {
    if (argc == 1) usage("Please provide subcommand.");
    if (strcmp(argv[1], "coverage") == 0) return samvt_coverage(argc - 1, argv + 1);
    else if (strcmp(argv[1], "mutation") == 0) return samvt_mutation(argc - 1, argv + 1);
    else usage("Unrecognized subcommand.");
    return (0);
}