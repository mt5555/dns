#include <stdio.h>
#include <errno.h>
extern int errno;

int main(int argc, char **argv) {
  char *fname;
  FILE *fid, *fid2;
  float bufx1,bufx2,bufy1,bufy2,bufz1,bufz2,mn,mx;
  char subname[240];
  char c='A';
  int i,j,k,ioff,joff,koff;
  long nbig,ncube,size;

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


  size=nbig*nbig;
  bufx1=(float *)malloc(size*4);
  bufx2=(float *)malloc(size*4);
  bufy1=(float *)malloc(size*4);
  bufy2=(float *)malloc(size*4);
  bufz1=(float *)malloc(size*4);
  bufz2=(float *)malloc(size*4);

  if (bufz2==NULL) {
    printf("Error: couldn't malloc subcube \n");
    exit(1);
  }


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

	mn=buf[0];
        mx=buf[0];
	for (i=0; i<size; ++i) {
	  if (buf[i]<mn) mn=buf[i];
	  if (buf[i]>mx) mx=buf[i];
	}

	/* compute name of subcube: */
	strcpy(subname,fname);
	i=strlen(subname);
	subname[i++]='.';
	subname[i++]=c+(ioff/ncube);
	subname[i++]=c+(joff/ncube);
	subname[i++]=c+(koff/ncube);
        subname[i++]=0;
	printf("%s min=%f max=%f\n",subname,mn,mx);
	fid2=fopen(subname,"w");
	fwrite(buf,4,nbig*nbig,fid2);
	fclose(fid2);
	goto done;
      }
    }
  }
done:
  fclose(fid);
  return 0;
}

