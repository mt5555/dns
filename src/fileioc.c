/*

Fortran callable wrappers to C binary I/O 

 */

#include <stdio.h>
#include <errno.h>
extern int errno;

void copen_(char *fname,char *mode,FILE **fid,int *err ) {
    char cname[80];     
    char cmode[2];
    cmode[0]=mode[0];
    cmode[1]=0;

    *fid=fopen(fname,&cmode);
    *err=0;
    if (*fid==NULL) *err=errno;
}

void cclose_(FILE **fid) {
    fclose(*fid);
}


/*
f  Write 8 byte numbers 
 */
void cwrite8_(FILE **fid,char *buf,int *len) {
    fwrite(buf,8,*len,*fid);
}

/*
  Write 8 byte numbers 
 */
void cread8_(FILE **fid, char *buf,int *len) {
    fread(buf,8,*len,*fid);
}
