#include <stdio.h>
#include <errno.h>
extern int errno;

int main(int argc, char **argv) {
  char fname[240];
  FILE *fid[3][3], *fid2, *fid3;
  double *buf,mn,mx;
  double sum[3][3],sum2[3][3]; 
  char subname[240];
  char c='A';
  int i,j,ioff,joff,koff,len,gradi,gradj;
  long nbig,ncube,size;

  if (argc != 4) {
    printf("USAGE:  ./extract_gradu filename  orig_size  subcube_size \n");
    exit(1);
  }
  /* open files of the form:  filename.ij */
  strncpy(fname,argv[1],240);
  len=strlen(fname);
  fname[len]='.';
  fname[len+3]=0;
  for (gradi=0; gradi<3; ++gradi) {
      for (gradj=0; gradj<3; ++gradj) {
          sprintf(&fname[len+1],"%1i",1+gradi);
          sprintf(&fname[len+2],"%1i",1+gradj);
          fid[gradi][gradj]=fopen(fname,"r");
          printf("filename = %i %i %s %p\n",gradi,gradj,fname,fid[gradi][gradj]);
          if (fid[gradi][gradj]==NULL) {
              printf("error opening file...\n");
              exit(1);
          }
      }
  }
  fname[len]=0;



  /* compute name of subcube: */
  strcpy(subname,fname);
  len=strlen(subname);
  strcpy(&subname[len],".gradu");
  printf("subname = %s \n",subname);
  fid2=fopen(subname,"w");

  strcpy(&subname[len],".gradu2");
  fid3=fopen(subname,"w");
  printf("subname = %s \n",subname);
  

  /* 256^3 blocks: */
  nbig=atoi(argv[2]);
  ncube=atoi(argv[3]);
  printf("input file: %s\n ",fname);
  printf("original size: %i  subcube size: %i \n",nbig,ncube); 

  size=ncube;
  size=size*size*size;    
  buf=(double *)malloc(size*sizeof(buf[1]));
  if (buf==NULL) {
    printf("Error: couldn't malloc subcube \n");
    exit(1);
  }



  for (ioff=0; ioff <= nbig-ncube;  ioff+=ncube ) {
    for (joff=0; joff <= nbig-ncube;  joff+=ncube ) {
      for (koff=0; koff <= nbig-ncube;  koff+=ncube ) {


        printf("%i %i %i \n",ioff/ncube,joff/ncube,koff/ncube);
        for (gradi=0; gradi<3; ++gradi) {
            for (gradj=0; gradj<3; ++gradj) {

                printf("u_%1i,%1i ",gradi,gradj);
                for (i=0; i<ncube; ++i) {
                    for (j=0; j<ncube; ++j) {
                        long pos_big=nbig*(nbig*(i+ioff) + j +joff) + koff;
                        long pos_small=ncube*(ncube*i + j) ;
                        pos_big = sizeof(buf[0])*pos_big;
                        fseek(fid[gradi][gradj],pos_big,SEEK_SET);
                        if (ncube!=fread(&buf[pos_small],sizeof(buf[0]),ncube,
                                         fid[gradi][gradj])) {
                            printf("Error on read \n"); 
                            exit(1);
                        }
                    }
                }
                
                mn=buf[0];
                mx=buf[0];
                sum[gradi][gradj]=0;
                sum2[gradi][gradj]=0;
                for (i=0; i<size; ++i) {
                    sum[gradi][gradj] += buf[i];  
                    sum2[gradi][gradj] += buf[i]*buf[i];  
                    if (buf[i]<mn) mn=buf[i];
                    if (buf[i]>mx) mx=buf[i];
                }
                sum[gradi][gradj]/=size;
                sum2[gradi][gradj]/=size;
                
                printf("min=%f max=%f\n",mn,mx);


            }
        }

        fprintf(fid2,"%i %i %i \n",ioff/ncube,joff/ncube,koff/ncube);
        for (i=0; i<3; ++i) {
            fprintf(fid2,"%18.10e %18.10e %18.10e \n",
                    sum[i][0],sum[i][1],sum[i][2]);
        }


        fprintf(fid3,"%i %i %i \n",ioff/ncube,joff/ncube,koff/ncube);
        for (i=0; i<3; ++i) {
            fprintf(fid3,"%18.10e %18.10e %18.10e \n",
                    sum2[i][0],sum2[i][1],sum2[i][2]);
        }


      }
    }
  }
done:
  for (gradi=0; gradi<3; ++gradi) {
      for (gradj=0; gradj<3; ++gradj) {
          fclose(fid[gradi][gradj]);
      }
  }
  fclose(fid2);

  return 0;
}






