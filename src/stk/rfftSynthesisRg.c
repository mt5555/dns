/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftSynthesisRg.c,v 1.1 2001-08-06 20:33:35 mt Exp $
 */

#include "rfft.h"

void rfftSynthesisRg(real *y,real *x,complex *R,complex *T,
		complex omega,int L,int N,int r)
{
  complex w,z,a,b,W;
  real s;
  int l,k,j,q;
  
  for(j=0;j<N;j++)
    {
      T[0].re=y[0+j*r*L];
      T[0].im=0.0;
      for(l=1;l<(r+1)/2;l++)
	{
	  T[l].re=y[2*l*L-1+j*r*L];
	  T[l].im=y[2*l*L  +j*r*L];
	}
      if(r%2==0)
	{
	  T[r/2].re=y[r*L-1+j*r*L];
	  T[r/2].im=0.0;
	}
      s=T[0].re;
      for(q=1;q<(r+1)/2;q++)
	s=s+2.0*T[q].re;
      if(r%2==0)
	s=s+T[r/2].re;
      x[j*L]=s;
      for(q=1;q<r;q++)
	{
	  s=T[0].re;
	  w.re=R[q].re;
	  w.im=R[q].im;
	  for(l=1;l<(r+1)/2;l++)
	    {
	      s=s+2.0*(w.re*T[l].re+w.im*T[l].im);
	      a.re=w.re*R[q].re-w.im*R[q].im;
	      a.im=w.re*R[q].im+w.im*R[q].re;
	      w.re=a.re;
	      w.im=a.im;
	    }
	  if(r%2==0)
	    if(q%2==0)
	      s=s+T[r/2].re;
	    else
	      s=s-T[r/2].re;
	  x[q*N*L+j*L]=s;
	}
      w.re=omega.re+1.0;                   
      w.im=omega.im;
      for(k=1;k<(L+1)/2;k++)
	{
	  T[0].re=y[2*k-1+j*r*L];
	  T[0].im=y[2*k  +j*r*L];
	  for(l=1;l<(r+1)/2;l++)
	    {
	      T[l].re=y[2*(k+l*L)-1+j*r*L];
	      T[l].im=y[2*(k+l*L)  +j*r*L];
	      T[l+(r+1)/2-1].re= y[2*(L-k+(l-1)*L)-1+j*r*L];
	      T[l+(r+1)/2-1].im=-y[2*(L-k+(l-1)*L)  +j*r*L];
	    }
	  if(r%2==0)
	    {
	      T[r-1].re= y[2*(L-k+(r/2-1)*L)-1+j*r*L];
	      T[r-1].im=-y[2*(L-k+(r/2-1)*L)  +j*r*L];
	    }
	  W.re=1.0;
	  W.im=0.0;
	  for(q=0;q<r;q++)
	    {
	      z.re=T[0].re;
	      z.im=T[0].im;
	      a.re=R[q].re;
	      a.im=R[q].im;
	      for(l=1;l<(r+1)/2;l++)
		{
		  z.re=z.re+T[l].re*a.re+T[l].im*a.im;
		  z.im=z.im+T[l].im*a.re-T[l].re*a.im;
		  b.re=a.re*R[q].re-a.im*R[q].im;
		  b.im=a.im*R[q].re+a.re*R[q].im;
		  a.re=b.re;
		  a.im=b.im;
		}
	      a.re=R[q].re;
	      a.im=R[q].im;
	      for(l=(r+1)/2;l<r;l++)
		{
		  z.re=z.re+T[l].re*a.re-T[l].im*a.im;
		  z.im=z.im+T[l].im*a.re+T[l].re*a.im;
		  b.re=a.re*R[q].re-a.im*R[q].im;
		  b.im=a.im*R[q].re+a.re*R[q].im;
		  a.re=b.re;
		  a.im=b.im;
		}
	      x[2*k-1+q*N*L+j*L]=z.re*W.re+z.im*W.im;
	      x[2*k  +q*N*L+j*L]=z.im*W.re-z.re*W.im;
	      z.re=W.re*w.re-W.im*w.im;
	      z.im=W.im*w.re+W.re*w.im;
	      W.re=z.re;
	      W.im=z.im;
	    }
	  z.re=omega.re*w.re-omega.im*w.im+w.re;
	  z.im=omega.re*w.im+omega.im*w.re+w.im; 
	  w.re=z.re;
	  w.im=z.im;
	}
      if(L%2==0)
	{
	  for(l=0;l<r/2;l++)
	    {
	      T[l].re=y[L+2*l*L-1+j*r*L];
	      T[l].im=y[L+2*l*L  +j*r*L];
	    }
	  if(r%2==1)
	    {
	      T[(r-1)/2].re=y[r*L-1+j*r*L];
	      T[(r-1)/2].im=0.0;
	    }
	  for(l=1;l<=r/2;l++)
	    {
	      T[r-l].re= T[l-1].re;
	      T[r-l].im=-T[l-1].im;
	    }
	  W.re=1.0;
	  W.im=0.0;
	  for(q=0;q<r;q++)
	    {
	      a.re=R[q].re;
	      a.im=R[q].im;
	      z.re=T[0].re;
	      z.im=T[0].im;
	      for(l=1;l<r;l++)
		{
		  z.re=z.re+T[l].re*a.re+T[l].im*a.im;
		  z.im=z.im+T[l].im*a.re-T[l].re*a.im;
		  b.re=a.re*R[q].re-a.im*R[q].im;
		  b.im=a.im*R[q].re+a.re*R[q].im;
		  a.re=b.re;
		  a.im=b.im;
		}
	      x[L-1+q*N*L+j*L]=z.re*W.re+z.im*W.im;
	      z.re=W.re*w.re-W.im*w.im;
	      z.im=W.im*w.re+W.re*w.im;
	      W.re=z.re;
	      W.im=z.im;
	    }
	}
    }
}
