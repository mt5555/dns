#include <stdio.h>
#include <errno.h>
extern int errno;

int main(int argc, char **argv) {
  char *fname;
  FILE *fid, *fid2;
  double *buf_avg,*buf;
  double mn,mx,sum; 
  char subname[240];
  char c='A';
  int i,j,k,ioff,joff,koff;
  long ncube,nbig,navg,size,size_avg;
  int ws=sizeof(buf[0]);

  if (argc != 4) {
    printf("USAGE:  ./extract filename  orig_size  averaging_size \n");
    exit(1);
  }
  fname=argv[1];
  fid=fopen(fname,"r");

  printf("using real*%i\n",ws);
 
  /* 256^3 blocks: */
  nbig=atoi(argv[2]);
  navg=atoi(argv[3]);
  printf("input file: %s\n ",fname);
  printf("original size: %i  subcube size: %i \n",nbig,navg); 

  size=navg;
  size=size*size*size;    
  buf=(double *)malloc(size*ws);
  if (buf==NULL) {
    printf("Error: couldn't malloc subcube \n");
    exit(1);
  }
 
  ncube=nbig/navg;
  size_avg=nbig/navg;
  size_avg=size_avg*size_avg*size_avg;
  buf_avg=(double *)malloc(size_avg*ws);
  if (buf_avg==NULL) {
    printf("Error: couldn't malloc average subcube \n");
    exit(1);
  }


  for (ioff=0; ioff <= nbig-navg;  ioff+=navg ) {
    for (joff=0; joff <= nbig-navg;  joff+=navg ) {
      for (koff=0; koff <= nbig-navg;  koff+=navg ) {
	k=0;
	for (i=0; i<navg; ++i) {
	  for (j=0; j<navg; ++j) {
	    long pos_big=nbig*(nbig*(i+ioff) + j +joff) + k+koff;
	    long pos_small=navg*(navg*i + j) + k;
	    fseek(fid,ws*pos_big,SEEK_SET);
	    if (navg!=fread(&buf[pos_small],ws,navg,fid)) {
	      printf("Error on read \n"); exit(1);
	    }
	  }
	}

        sum=0;
        mn=buf[0]; mx=buf[0]; 
	for (i=0; i<size; ++i) {
          sum=sum+buf[i];  
          if (buf[i]<mn) mn=buf[i];
          if (buf[i]>mx) mx=buf[i];
	}

        sum=sum/size;
        {  /* index:  ioff/navg, joff/navg, koff/navg.  size: ncube  */
            long pos_small=ncube*(ncube*ioff/navg + joff/navg) + koff/navg;
            buf_avg[pos_small]=sum;
        }
        printf("averaged subcube: %i %i %i  mn=%f  mx=%f\n",ioff/navg,joff/navg,koff/navg,mn,mx);
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
  fwrite(buf_avg,ws,size_avg,fid2);
  fclose(fid2);

  return 0;
}






