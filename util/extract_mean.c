#include <stdio.h>
#include <errno.h>
extern int errno;

int main(int argc, char **argv) {
  char *fname;
  FILE *fid, *fid2;
  double *buf_avg,*buf;
  double sum; 
  char subname[240];
  char c='A';
  int i,j,k,ioff,joff,koff;
  long nbig,ncube,size,size_avg;

  if (argc != 4) {
    printf("USAGE:  ./extract filename  orig_size  subcube_size \n");
    exit(1);
  }
  fname=argv[1];
  fid=fopen(fname,"r");


  /* 256^3 blocks: */
  nbig=atoi(argv[2]);
  ncube=atoi(argv[3]);
  printf("input file: %s\n ",fname);
  printf("original size: %i  subcube size: %i \n",nbig,ncube); 

  size=ncube;
  size=size*size*size;    
  buf=(double *)malloc(size*4);
  if (buf==NULL) {
    printf("Error: couldn't malloc subcube \n");
    exit(1);
  }

  size_avg=nbig/ncube;
  size_avg=size_avg*size_avg*size_avg;
  buf_avg=(double *)malloc(size_avg*4);
  if (buf_avg==NULL) {
    printf("Error: couldn't malloc average subcube \n");
    exit(1);
  }


  for (ioff=0; ioff <= nbig-ncube;  ioff+=ncube ) {
    for (joff=0; joff <= nbig-ncube;  joff+=ncube ) {
      for (koff=0; koff <= nbig-ncube;  koff+=ncube ) {
	k=0;
	for (i=0; i<ncube; ++i) {
	  for (j=0; j<ncube; ++j) {
	    long pos_big=nbig*(nbig*(i+ioff) + j +joff) + k+koff;
	    long pos_small=ncube*(ncube*i + j) + k;
	    fseek(fid,4*pos_big,SEEK_SET);
	    if (ncube!=fread(&buf[pos_small],4,ncube,fid)) {
	      printf("Error on read \n"); exit(1);
	    }
	  }
	}

        sum=0;
	for (i=0; i<size; ++i) {
          sum=sum+buf[i];  
	}
        sum=sum/size;
        { long pos_small=ncube*ioff + joff + koff/ncube;
        buf_avg[pos_small]=sum;
        }
        printf("averaged subcube: %i %i %i\n",ioff/ncube,joff/ncube,koff/ncube);
      }
    }
  }
done:
  fclose(fid);

  strcpy(subname,fname);
  i=strlen(subname);
  subname[i++]='.';
  subname[i++]='a';
  subname[i++]='v';
  subname[i++]='e';
  subname[i++]=0;
  fid2=fopen(subname,"w");
  fwrite(buf_avg,8,size_avg,fid2);
  fclose(fid2);

  return 0;
}






