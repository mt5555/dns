/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftSynthesis.c,v 1.2 2001-08-20 18:33:27 mt Exp $
 */

#include "rfft.h"

void rfftSynthesis(real *c,rfft *X)
{
  int m,l,L,N,odd,r,i,n;
  real *in,*out;
  rfftNode *Y;
	
  n=X->n;
  m=X->n;

  L=n;
  N=1;
  
  Y=(rfftNode *)(X->bot);
  l=0;
  while(m>1)
    {
      r=Y->r;
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
      L=L/r;
      switch(r)
	{
	case 2:
	  rfftSynthesisR2(in,out,Y->omega,L,N);
	  break;
	case 3:
	  rfftSynthesisR3(in,out,Y->omega,L,N);
	  break;
	case 4:
	  rfftSynthesisR4(in,out,Y->omega,L,N);
	  break;
	case 5:
	  rfftSynthesisR5(in,out,Y->omega,L,N);
	  break;
	case 8:
	  rfftSynthesisR8(in,out,Y->omega,L,N);
	  break;
	default:
	  rfftSynthesisRg(in,out,Y->aux,Y->aux+r,Y->omega,L,N,r);
	  break;
	}
      m=m/r;
      N=N*r;
      Y=(rfftNode *)(Y->prev);
    }
  if(odd)
    {
      for(i=0;i<n;i++)
	in[i]=out[i];
    }
}

void rfftSynthesis2(real *c,rfft *X)
{
  int m,l,L,N,odd,r,i,n;
  real *in,*out;
  rfftNode *Y;
	
  n=X->n;
  m=X->n;

  L=n;
  N=1;
  
  Y=(rfftNode *)(X->bot);
  l=0;
  while(m>1)
    {
      r=Y->r;
      l++;
      odd=l%2;
      if(odd)
	{
          if (l==1) {
              c[n]=c[1];
              c[1]=c[0];
              in=c+1;
          }else{
            in = c;
          }
	  out=X->work;
	}
      else
	{
	  in=X->work;
	  out=c;
	}
      L=L/r;
      switch(r)
	{
	case 2:
	  rfftSynthesisR2(in,out,Y->omega,L,N);
	  break;
	case 3:
	  rfftSynthesisR3(in,out,Y->omega,L,N);
	  break;
	case 4:
	  rfftSynthesisR4(in,out,Y->omega,L,N);
	  break;
	case 5:
	  rfftSynthesisR5(in,out,Y->omega,L,N);
	  break;
	case 8:
	  rfftSynthesisR8(in,out,Y->omega,L,N);
	  break;
	default:
	  rfftSynthesisRg(in,out,Y->aux,Y->aux+r,Y->omega,L,N,r);
	  break;
	}
      m=m/r;
      N=N*r;
      Y=(rfftNode *)(Y->prev);
    }
  if(odd)
    { 
      for(i=0;i<n;i++)
	c[i]=out[i];
    }

}

void rfftSynthesisM(real *c,rfft *X,int n,int p,int m)
{
  int i;
  c=c+p*n;
  for(i=p;i<p+m;i++)
    {
      rfftSynthesis(c,X);
      c=c+n;
    }
}
