/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: t_fft.c,v 1.1 2001-08-06 20:33:35 mt Exp $
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "cfft.h"
#include "rfft.h"

#define Max(A,B) ((A>B)?A:B)

double mydrand(void)
{
     double d = (double)rand();
     return(2.0*(d / (double) RAND_MAX - 0.5));
}

real t_comp(int n)
{
  complex *seq,*SEQ;
  cfft *W,*Wcopy;
  real l_inf,norm;
  int i;

  seq = (complex *)calloc((size_t)n,sizeof(complex));
  SEQ = (complex *)calloc((size_t)n,sizeof(complex));
  if((seq==NULL)||(SEQ==NULL))
    {
      fprintf(stderr,"cannot allocate arrays...exiting...\n");
      exit(-1);
    }
  for(i=0;i<n;i++)
    {
      seq[i].re = mydrand();
      seq[i].im = mydrand();
      SEQ[i].re = seq[i].re;
      SEQ[i].im = seq[i].im;
    }
	
  Wcopy=cfftInit(n);
  W=cfftCopy(Wcopy);
  if((Wcopy==NULL)||(W==NULL))
    {
      fprintf(stderr,"cannot allocate handle...exiting...\n");
      exit(-1);
    }
  cfftAnalysis(seq,W); 
  cfftSynthesis(seq,W);
  cfftDelete(W);
  cfftDelete(Wcopy);
  l_inf=0.0;
  for(i=0;i<n;i++)
    {
      norm=Max(fabs(seq[i].re/(real)n-SEQ[i].re),
		   fabs(seq[i].im/(real)n-SEQ[i].im));
      l_inf=Max(l_inf,norm);
    }
  free(seq);
  free(SEQ);
  return(l_inf);
}

real t_real(int n)
{
  real *seq,*SEQ;
  rfft *W,*Wcopy;
  real l_inf,norm;
  int i;

  seq = (real *)calloc((size_t)n,sizeof(real));
  SEQ = (real *)calloc((size_t)n,sizeof(real));
  if((seq==NULL)||(SEQ==NULL))
    {
      fprintf(stderr,"cannot allocate arrays...exiting...\n");
      return(-1);
    }
  for(i=0;i<n;i++)
    {
      seq[i] = mydrand();
      SEQ[i] = seq[i];
    }
	
  Wcopy=rfftInit(n);
  W=rfftCopy(Wcopy);
  rfftAnalysis(seq,W); 
  rfftSynthesis(seq,W);
  rfftDelete(W);
  rfftDelete(Wcopy);
  l_inf=0.0;
  for(i=0;i<n;i++)
    {
      norm=fabs(seq[i]/(real)n-SEQ[i]);
      l_inf=Max(l_inf,norm);
    }
  free(seq);
  free(SEQ);
  return(l_inf);
}

real t_cfft(int n,int k,int inverse,int odd)
{
  complex *seq;
  cfft *W,*Wcopy;
  real x,y;
  real l_inf,norm;
  int i;

  seq = (complex *)calloc((size_t)n,sizeof(complex));
  if(seq==NULL)
    {
      fprintf(stderr,"cannot allocate arrays...exiting...\n");
      exit(-1);
    }
      
  if(inverse)
    {
      for(i=0;i<n;i++)
	{
	  seq[i].re=0.0;
	  seq[i].im=0.0;
	}
      if(odd)
	{
	  seq[k].im=-0.5;
	  seq[n-k].im=0.5;
	}
      else
	{
	  seq[k].re=0.5;
	  seq[n-k].re=0.5;
	}
    }
  else
    {
      for(i=0;i<n;i++)
	{
	  x = ((real)i/(real)n)*2.0*STK_PI;
	  if(odd)
	    seq[i].re = sin((real)k*x);
	  else
	    seq[i].re = cos((real)k*x);
	  seq[i].im = 0.0;
	}
    }
  
  Wcopy=cfftInit(n);
  W=cfftCopy(Wcopy);
  if(inverse)
    cfftSynthesis(seq,W); 
  else
    cfftAnalysis(seq,W); 
  cfftDelete(Wcopy);
  cfftDelete(W);
  if(inverse)
    {
      l_inf=0.0;
      for(i=0;i<n;i++)
	{
	  x=((real)i/(real)n)*2.0*STK_PI;
	  y=odd?sin((real)k*x):cos((real)k*x);
	  norm=Max(fabs(seq[i].re-y),fabs(seq[i].im));
	  l_inf=Max(l_inf,norm);
	}	  
    }
  else
    {
      if(odd)
	{
	  seq[k].im=seq[k].im+(double)n/2.0;
	  seq[n-k].im=seq[n-k].im-(double)n/2.0;
	}
      else
	{
	  seq[k].re=seq[k].re-(double)n/2.0;
	  seq[n-k].re=seq[n-k].re-(double)n/2.0;
	}
      l_inf=0.0;
      for(i=0;i<n;i++)
	{
	  norm=Max(fabs(seq[i].re),fabs(seq[i].im))/(double)n;
	  l_inf=Max(l_inf,norm);
	}
    }
  free(seq);
  return(l_inf);
}

real t_rfft(int n,int k,int inverse,int odd)
{
  real *seq;
  rfft *W,*Wcopy;
  real x,y;
  real l_inf,norm;
  int i;

  seq = (real *)calloc((size_t)n,sizeof(real));
  if(seq==NULL)
    {
      fprintf(stderr,"cannot allocate arrays...exiting...\n");
      exit(-1);
    }
  if(inverse)
    {
      for(i=0;i<n;i++)
	seq[i]=0.0;
      if(odd)
	seq[2*k]=-0.5;
      else
	seq[2*k-1]=0.5;
    }
  else
    {
      for(i=0;i<n;i++)
	{
	  x = ((real)i/(real)n)*2.0*STK_PI;
	  seq[i] = odd?sin(((real)k)*x):cos(((real)k)*x);
	}
    }
  Wcopy=rfftInit(n);
  W=rfftCopy(Wcopy);
  if(inverse)
    rfftSynthesis(seq,W);
  else
    rfftAnalysis(seq,W); 
  rfftDelete(W);
  rfftDelete(Wcopy);
  if(inverse)
    {
      l_inf=0.0;
      for(i=0;i<n;i++)
	{
	  x = ((real)i/(real)n)*2.0*STK_PI;
	  y = odd?sin((real)k*x):cos((real)k*x);
	  norm=fabs(seq[i]-y);
	  l_inf=Max(l_inf,norm);
	}
    }
  else
    {
      if(odd)
	seq[2*k]+=(double)n/2.0;
      else
	seq[2*k-1]-=(double)n/2.0;
      l_inf=0.0;
      for(i=0;i<n;i++)
	{
	  norm=fabs(seq[i])/(double)n;
	  l_inf=Max(l_inf,norm);
	}	
    }
  free(seq);
  return(l_inf);
}

real t_rfftm(int n,int m)
{
  real *seq,*c;
  rfft *W,*W_clone;
  real scale,l_inf,norm;
  int j;

  seq = (real *)calloc((size_t)(n*m),sizeof(real));
  c=(real *)calloc(m*n,sizeof(real));
  if((seq==NULL)||(c==NULL))
    {
      fprintf(stderr,"cannot allocate sequence array...exiting...\n");
      return(-1);
    }
      
  for(j=0;j<m*n;j++)
    {
      c[j]=mydrand();
      seq[j] = c[j];
    }
  
  W=rfftInit(n); 
  W_clone=rfftCopy(W);
  if(W==NULL || W_clone==NULL)
    {
      fprintf(stderr,"cannot allocate initialization arrays...exiting...\n");
      return(-1);
    }
  rfftAnalysisM(c,W_clone,n,0,m);
  rfftSynthesisM(c,W,n,0,m);
  rfftDelete(W_clone);
  rfftDelete(W);

  l_inf=0.0;
  scale=1.0/(real)n;
  for(j=0;j<m*n;j++)
    {
      norm=fabs(scale*c[j]-seq[j]);
      l_inf=Max(l_inf,norm);
    }
  free(c);
  free(seq);
  return(l_inf);
}

real t_cfftm(int n,int m)
{
  complex *seq,*c;
  cfft *W,*W_clone;
  real scale,l_inf,norm;
  int j;

  seq = (complex *)calloc((size_t)(n*m),sizeof(complex));
  c=(complex *)calloc(m*n,sizeof(complex));
  if((seq==NULL)||(c==NULL))
    {
      fprintf(stderr,"cannot allocate sequence array...exiting...\n");
      return(-1);
    }
  
  for(j=0;j<m*n;j++)
    {
      c[j].re=mydrand();
      c[j].im=mydrand();
      seq[j].re=c[j].re;
      seq[j].im=c[j].im;
    }
  
  W=cfftInit(n); 
  W_clone=cfftCopy(W);
  if(W==NULL || W_clone==NULL)
    {
      fprintf(stderr,"cannot allocate initialization arrays...exiting...\n");
      return(-1);
    }
  cfftAnalysisM(c,W_clone,n,0,m);
  cfftSynthesisM(c,W,n,0,m);
  cfftDelete(W_clone);
  cfftDelete(W);

  l_inf=0.0;
  scale=1.0/(real)n;
  for(j=0;j<m*n;j++)
    {
      norm=Max(fabs(scale*c[j].re-seq[j].re),fabs(scale*c[j].im-seq[j].im));
      l_inf=Max(l_inf,norm);
    }
  free(c);
  free(seq);
  return(l_inf);
}


void stkTest1a(int N)
{
  real e_max;
  real Er,Ec;
  int n,k;

  /* ==== 1D test, scans through all possible wave numbers ==== */

  Er=0.0;
  for(n=3;n<=N;n++)
    {
      e_max=0.0;
      for(k=1;k<n/2;k++)
	{
	  e_max=Max(t_rfft(n,k,0,1),e_max);
	  e_max=Max(t_rfft(n,k,1,1),e_max);
	  e_max=Max(t_rfft(n,k,0,0),e_max);
	  e_max=Max(t_rfft(n,k,1,0),e_max);
	} 
      Er=Max(Er,e_max);
    }
  printf("1D Test1 n=[3:%d] \t:: Real    l_inf=%le\n",N,Er);
  Ec=0.0;
  for(n=3;n<=N;n++)
    {
      e_max=0.0;
      for(k=1;k<n/2;k++)
	{
	  e_max=Max(t_cfft(n,k,0,1),e_max);
	  e_max=Max(t_cfft(n,k,1,1),e_max);
	  e_max=Max(t_cfft(n,k,0,0),e_max);
	  e_max=Max(t_cfft(n,k,1,0),e_max);
	} 
      Ec=Max(Ec,e_max);
    }
  printf("1D Test1 n=[3:%d] \t:: Complex l_inf=%le\n",N,Ec);
}

void stkTest1b(int n1)
{
  real l_inf;
  real Er,Ec;
  int n;

  /* ==== 1D test, loops over range of nx ==== */

  Er=0.0;
  for(n=3;n<=n1;n++)
    {
      l_inf=t_real(n);
      Er=Max(Er,l_inf);
    }
  printf("1D Test2 n=[3:%d] \t:: Real    l_inf=%le\n",n1,Er);
  Ec=0.0;
  for(n=3;n<=n1;n++)
    {
      l_inf=t_comp(n);
      Ec=Max(Ec,l_inf);
    }
  printf("1D Test2 n=[3:%d] \t:: Complex l_inf=%le\n",n1,Ec);

}  

void stkTest1m(int n1,int m1)
{
  real l_inf;
  real Er,Ec;
  int n;
  int m;

  m=m1;

  /* ==== 1D test, loops over range of nx ==== */

  Er=0.0;
  for(n=3;n<=n1;n++)
    {
      l_inf=t_rfftm(n,m);
      Er=Max(Er,l_inf);
    }
  printf("1D Test2 n=[3:%d] \t:: Real    l_inf=%le\n",n1,Er);
  Ec=0.0;
  for(n=3;n<=n1;n++)
    {
      l_inf=t_cfftm(n,m);
      Ec=Max(Ec,l_inf);
    }
  printf("1D Test2 n=[3:%d] \t:: Complex l_inf=%le\n",n1,Ec);

}  

int main(int argc,char **argv)
{
  int N=Max(RADIX_8_MIN,64);
  int n=Max(RADIX_4_MIN,64);
  int n2=32;
  int n3=16;
  if(argc==2)
    N=atoi(argv[1]);
  printf("radix-8 cutoff: %d  radix-4 cutoff: %d\n",RADIX_8_MIN,RADIX_4_MIN);
  printf("performing wavenumber test up to N=%d...\n",n);
  stkTest1a(n);
  printf("performing random data test up to N=%d...\n",N);
  stkTest1b(N);
  printf("preforming multiple instance tests up to N=%d...\n",n);
  stkTest1m(n,24);
  return(0);
}

