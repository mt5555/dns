/*

Fortran callable wrappers to C binary I/O 

 */

#if (defined AIX || defined HPUX || defined OSF1)
#define FORTRAN(A) A
#else
#define FORTRAN(A) A##_
#endif




#include <stdio.h>
#include <errno.h>
extern int errno;

void FORTRAN(copen) (char fname[80],char mode[1],FILE **fid,int *err ) {
    int i;
    char cname[80];     
    char cmode[2];
    cmode[0]=mode[0];
    cmode[1]=0;

    strncpy(cname,fname,80);
    for (i=0; i<80; ++i) {
        if (cname[i]==' ') cname[i]=0;
    }
    
    *fid=fopen(cname,cmode);
    *err=0;
    if (*fid==NULL) *err=errno;
}

void FORTRAN(cclose) (FILE **fid) {
    fclose(*fid);
}


/*
  Write 8 byte numbers 
*/
void FORTRAN(cwrite8) (FILE **fid,char *buf,int *len) {
    fwrite(buf,8,*len,*fid);
}

/*
  Write 8 byte numbers 
 */
void FORTRAN(cread8) (FILE **fid, char *buf,int *len) {
    fread(buf,8,*len,*fid);
}
