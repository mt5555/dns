#include <stdio.h>
#include <errno.h>
#include <math.h>
extern int errno;

int main(int argc, char **argv) {
  char *fname;
  FILE *fidu,*fidv,*fidw;
  int nx=128;
  int i,j,k;
  double pi,x,y,z,buf;

  fname="test.u";
  printf("name = %s \n",fname);
  fidu=fopen(fname,"w");
  fname="test.v";
  printf("name = %s \n",fname);
  fidv=fopen(fname,"w");
  fname="test.w";
  printf("name = %s \n",fname);
  fidw=fopen(fname,"w");

  pi=1; pi=4*atan(pi);

  // write header information: (just dummy data of all zeros)
  buf=0;
  for (i=0; i<4+3*nx; ++i) {
      fwrite(&buf,8,1,fidu);
      fwrite(&buf,8,1,fidv);
      fwrite(&buf,8,1,fidw);
  }

  // write a nx^3 brick of floating point data:
  // dummy data u(i,j,k) = i+j+k
  for (k=0; k<=nx; ++k) {
      for (j=0; j<=nx; ++j) {
          for (i=0; i<=nx; ++i) {
              x = (2*pi*i)/nx;
              y = (2*pi*j)/nx;
              z = (2*pi*k)/nx;
              buf = sin(x); fwrite(&buf,8,1,fidu);
              buf = sin(y); fwrite(&buf,8,1,fidv);
              buf = sin(z); fwrite(&buf,8,1,fidw);
          }
      }
  }
  fclose(fidu);
  fclose(fidv);
  fclose(fidw);
  return 0;
}






