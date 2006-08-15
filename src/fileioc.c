/*

Fortran callable wrappers to C binary I/O 

 */

#ifdef F_NO_UNDERSCORE
#define FORTRAN(A) A
#else
#define FORTRAN(A) A##_
#endif




#include <string.h>
#include <stdio.h>
#include <errno.h>
extern int errno;


static int byteswap_input=0;

void FORTRAN(set_byteswap_input)(int *val) {
   byteswap_input=*val;
}



void FORTRAN(copen) (char fname[80],char mode[1],FILE **fid,int *err ) {
    int i;
    char cname[280];     
    char cmode[2];
    cmode[0]=mode[0];
    cmode[1]=0;

    strncpy(cname,fname,280);
    for (i=0; i<280; ++i) {
        if (cname[i]==' ') cname[i]=0;
    }
    
    *fid=fopen(cname,cmode);
    *err=0;
    if (*fid==NULL) {
       *err=errno;
       perror("copen():");
       printf("Error in copen mode=%i fname=%s\n",cmode,fname);
    }

}

void FORTRAN(cclose) (FILE **fid,int *err) {
    *err=fclose(*fid);
    
}


/*
  Write 8 byte floats
*/
void FORTRAN(cwrite8) (FILE **fid,char *buf,int *len) {
    if (*len>0) fwrite(buf,8,*len,*fid);
}
/*
  Write 4 byte floats
*/
void FORTRAN(cwrite4) (FILE **fid,char *buf,int *len) {
    if (*len>0) fwrite(buf,4,*len,*fid);
}

/*
  Write 1 byte chars
*/
void FORTRAN(cwrite1) (FILE **fid,char *buf,int *len) {
    if (*len>0) fwrite(buf,1,*len,*fid);
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


void byteswap4(char *buf,int len) {
    int i,j;
    char swap[4];
    for (i=0; i<len; ++i) {
        for (j=0; j<4; ++j) {
            swap[j]=buf[4*i + 3-j];
        }
        for (j=0; j<4; ++j) {
            buf[4*i + j]=swap[j];
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



/*
  Read 4 byte numbers , also return number of numbers read
 */
void FORTRAN(cread4e) (FILE **fid, char *buf,int *len,int *olen) {
    *olen=fread(buf,4,*len,*fid);
    if (byteswap_input) byteswap4(buf,*olen);
}

void FORTRAN(cread1e) (FILE **fid, char *buf,int *len,int *olen) {
    *olen=fread(buf,1,*len,*fid);
}



