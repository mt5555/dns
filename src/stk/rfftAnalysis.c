/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftAnalysis.c,v 1.1 2001-08-06 20:33:33 mt Exp $
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
