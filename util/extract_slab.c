#include <stdio.h>
#include <errno.h>
extern int errno;


void write_ensheader(int partno,FILE *fid) {
  char text[80];
  memset(text,0,80);   strcpy(text,"part");
  fwrite(text,1,80,fid);
  fwrite(&partno,4,1,fid);
  memset(text,0,80);   strcpy(text,"block");
  fwrite(text,1,80,fid);
}



int main(int argc, char **argv) {
  char *fname;
  FILE *fid, *fid2;
  float *buf,*bufx1,*bufx2,*bufy1,*bufy2,*bufz1,*bufz2,mn,mx;
  char subname[240];
  char c='A';
  int i,j,k,ioff,joff,koff,partno;
  long nbig,size;

  if (argc != 4) {
    printf("USAGE:  ./extract filename  orig_size part_no_start \n");
    exit(1);
  }
  fname=argv[1];
  fid=fopen(fname,"r");


  /* 256^3 blocks: */
  nbig=atoi(argv[2]);
  partno=atoi(argv[3]);
  printf("input file: %s\n ",fname);
  printf("original size: %i  \n",nbig);
  printf("part number: %i  \n",partno);


  size=nbig*nbig;
  buf=(float *)malloc(nbig*4);
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


  for (i=0; i<nbig; ++i) {
    printf("reading i=%i\n",i);
    for (j=0; j<nbig; ++j) {
      if (nbig!=fread(buf,4,nbig,fid)) {
	printf("Error on read \n"); exit(1);
      }
      for (k=0; k<nbig; ++k) {
	// data at [i,j,k] = i*nbig*nbig + j*nbig + k = buf[k]
	if (i==0) bufx1[nbig*j+k]=buf[k];
	if (i==nbig-1) bufx2[nbig*j+k]=buf[k];
	if (j==0) bufy1[nbig*i+k]=buf[k];
	if (j==nbig-1) bufy2[nbig*i+k]=buf[k];
	if (k==0) bufz1[nbig*i+j]=buf[k];
	if (k==nbig-1) bufz2[nbig*i+j]=buf[k];
      }
    }
  }
  fclose(fid);



  printf("writing slabs...\n");
  fid=fopen("sides.ens","w");


  // write out the data
  write_ensheader(partno++,fid);
  fwrite(bufx1,4,nbig*nbig,fid);

  write_ensheader(partno++,fid);
  fwrite(bufx2,4,nbig*nbig,fid);

  write_ensheader(partno++,fid);
  fwrite(bufy1,4,nbig*nbig,fid);

  write_ensheader(partno++,fid);
  fwrite(bufy2,4,nbig*nbig,fid);

  write_ensheader(partno++,fid);
  fwrite(bufz1,4,nbig*nbig,fid);

  write_ensheader(partno++,fid);
  fwrite(bufz2,4,nbig*nbig,fid);

  fclose(fid);
  return 0;
}





