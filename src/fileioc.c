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


static int byteswap_input=0;

void FORTRAN(set_byteswap_input)(int val) {
   byteswap_input=val;
}



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

void FORTRAN(cclose) (FILE **fid,int *err) {
    *err=fclose(*fid);
    
}


/*
  Write 8 byte numbers 
*/
void FORTRAN(cwrite8) (FILE **fid,char *buf,int *len) {
    if (*len>0) fwrite(buf,8,*len,*fid);
}




void byteswap8(char *buf,int len) {
    int i,j;
    char swap[8];
    for (i=0; i<len; ++i) {
        for (j=0; j<8; ++j) {
            swap[j]=buf[8*i + 7-j];
        }
        for (j=0; j<8; ++j) {
            buf[8*i + j]=swap[j];
        }
    }
}


/*
  Read 8 byte numbers 
 */
void FORTRAN(cread8) (FILE **fid, char *buf,int *len) {
    int n;
    n=fread(buf,8,*len,*fid);
    if (byteswap_input) byteswap8(buf,n);
}

/*
  Read 8 byte numbers , also return number of numbers read
 */
void FORTRAN(cread8e) (FILE **fid, char *buf,int *len,int *olen) {
    *olen=fread(buf,8,*len,*fid);
    if (byteswap_input) byteswap8(buf,*olen);
}

