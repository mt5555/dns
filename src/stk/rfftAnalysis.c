/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftAnalysis.c,v 1.4 2001-08-23 18:49:13 mt Exp $
 */

#include "rfft.h"

void rfftAnalysis(real *c,rfft *X)
{
  int i,n,m,l,L,N,odd,r;
  real *in,*out;
  rfftNode *Y;

  n=X->n;
  N=n;
  m=n;
  L=1;
  Y=(rfftNode *)(X->top);
  l=0;
  while(m>1)
    {
      r=Y->r;
      N=N/r;
      l++;
      odd=l%2;
      if(odd)
	{
          in=c;
	  out=X->work;
	}
      else
	{
	  in=X->work;
	  out=c;
	}
      switch(r)
	{
	case 2:
	  rfftAnalysisR2(in,out,Y->omega,L,N);
	  break;
	case 3:
	  rfftAnalysisR3(in,out,Y->omega,L,N);
	  break;
	case 4:
	  rfftAnalysisR4(in,out,Y->omega,L,N);
	  break;
	case 5:
	  rfftAnalysisR5(in,out,Y->omega,L,N);
	  break;
	case 8:
	  rfftAnalysisR8(in,out,Y->omega,L,N);
	  break; 
	default:
	  rfftAnalysisRg(in,out,Y->aux,Y->aux+r,Y->omega,L,N,r);
	  break;
	}
      m=m/r;
      L=L*r;
      Y=(rfftNode *)(Y->next);
    }
    if(odd)
      {
	for(i=0;i<n;i++)
	  in[i]=out[i];
      }
}



void rfftAnalysis2(real *c,rfft *X)
{
  int i,n,m,l,L,N,odd,r;
  real *in,*out;
  rfftNode *Y;

  n=X->n;
  N=n;
  m=n;
  L=1;
  Y=(rfftNode *)(X->top);
  l=0;
  while(m>1)
    {
      r=Y->r;
      N=N/r;
      l++;
      odd=l%2;
      if(odd)
	{
	  in = (l==1) ? c : c+1;
	  out=X->work;
	}
      else
	{
	  in=X->work;
	  out=c+1;
	}
      switch(r)
	{
	case 2:
	  rfftAnalysisR2(in,out,Y->omega,L,N);
	  break;
	case 3:
	  rfftAnalysisR3(in,out,Y->omega,L,N);
	  break;
	case 4:
	  rfftAnalysisR4(in,out,Y->omega,L,N);
	  break;
	case 5:
	  rfftAnalysisR5(in,out,Y->omega,L,N);
	  break;
	case 8:
	  rfftAnalysisR8(in,out,Y->omega,L,N);
	  break; 
	default:
	  rfftAnalysisRg(in,out,Y->aux,Y->aux+r,Y->omega,L,N,r);
	  break;
	}
      m=m/r;
      L=L*r;
      Y=(rfftNode *)(Y->next);
    }
  /*
          input:    1 2 3 4 5 6 7 8 *     (* = 1 extra word of storage required)
         output:    0 1 1 2 2 3 3 4 *    
re-order output:    0 4 1 1 2 2 3 3 *       
        
   */
    if(odd)
      { /* answer is in out=X->work. copy into C, with re-order */
        c[0]=out[0];
        c[1]=out[n-1];
	for(i=2;i<n;i++)
	  c[i]=out[i-1];
      }
    else
      { /* answer is in out=c+1.  do a little re-ordering */
          c[0]=c[1] /* = out[0] */;
          c[1]=c[n] /* = out[n-1] */;  
      }
}




void rfftAnalysisM(real *c,rfft *X,int n,int p,int m)
{
  int i;
  c=c+p*n;
  for(i=0;i<m;i++)
    {
      rfftAnalysis(c,X);
      c=c+n;
    }
}

void rfftAnalysisM2(real *c,rfft *X,int n,int p,int m)
{
  int i;
  c=c+p*n;
  for(i=0;i<m;i++)
    {
      rfftAnalysis2(c,X);
      c=c+n;
    }
}
