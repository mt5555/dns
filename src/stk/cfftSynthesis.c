/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: cfftSynthesis.c,v 1.1 2001-08-06 20:33:31 mt Exp $
 */

#include "cfft.h"

void cfftSynthesis(complex *c,cfft *X)
{
  int m,l,L,N,odd,r;
  complex *in,*out;
  cfftNode *Y;

  N=X->n;
  m=X->n;
  L=1;
  Y=(cfftNode *)(X->top);
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
	  out=(N==1)?in:X->work;
	}
      else
	{
	  in=X->work;
	  out=c;
	}
      switch(r)
	{
	case 2:
	  cfftSynthesisR2(in,out,Y->omega,L,N);
	  break;
	case 3:
	  cfftSynthesisR3(in,out,Y->omega,L,N);
	  break;
	case 4:
	  cfftSynthesisR4(in,out,Y->omega,L,N);
	  break;
	case 5:
	  cfftSynthesisR5(in,out,Y->omega,L,N);
	  break;
	case 8:
	  cfftSynthesisR8(in,out,Y->omega,L,N);
	  break;
	default:
	  cfftSynthesisRg(in,out,Y->aux,Y->aux+r,Y->omega,L,N,r);
	  break;
	}
      m=m/r;
      L=L*r;
      Y=(cfftNode *)(Y->next);
    }
}

void cfftSynthesisM(complex *c,cfft *X,int n,int p,int m)
{
  int i;
  for(i=p;i<p+m;i++)
    {
      cfftSynthesis(c+i*n,X);
    }
}
