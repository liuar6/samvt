#include <string.h>
static char ** strsplit(char * line, char ** results, int length, char c){
    char *start=line;
    char *end=NULL;
    int i=0;
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
    for (;i<length;++i) results[i]=NULL;
    return results;
}