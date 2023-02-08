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

static void seq_complement(char *seq) {
    char *end=seq+strlen(seq);
    for (char *p=seq; p!=end; ++p) {
        switch (*p) {
            case 'A':  *p = 'T';break;
            case 'U':  *p = 'A';break;
            case 'T':  *p = 'A';break;
            case 'C':  *p = 'G';break;
            case 'G':  *p = 'C';break;
            case 'a':  *p = 't';break;
            case 'u':  *p = 'a';break;
            case 't':  *p = 'a';break;
            case 'c':  *p = 'g';break;
            case 'g':  *p = 'c';break;
        }
    }
}

static void seq_upper(char *seq){
    char *end=seq+strlen(seq);
    for (char *p=seq; p!=end; ++p) {
        switch (*p) {
            case 'a':  *p = 'A';break;
            case 'u':  *p = 'U';break;
            case 't':  *p = 'T';break;
            case 'c':  *p = 'C';break;
            case 'g':  *p = 'G';break;
        }
    }
}

static void seq_reverse(char *seq)
{
    uint length=strlen(seq);
    int halfLen=length>>1u;
    char *p=seq+length;
    char c;
    while (--halfLen>=0)
    {
        c=*seq;
        *seq++=*--p;
        *p=c;
    }
}
