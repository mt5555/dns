/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: cfftAnalysisRg.c,v 1.1 2001-08-06 20:33:31 mt Exp $
 */

#include "cfft.h"

void cfftAnalysisRg(complex *x,complex *y,complex *R,complex *T,
	   complex omega,int L,int N,int r)
{
  complex w,z,a,b;
  int l,k,j,i;

  for(j=0;j<N;j++)
    {
      w.re=1.0;
      w.im=0.0;
      for(i=0;i<L;i++)
	{
	  z.re=1.0;
	  z.im=0.0;
	  for(k=0;k<r;k++)
	    {
	      T[k].re=x[i+j*L+k*N*L].re*z.re-x[i+j*L+k*N*L].im*z.im;
	      T[k].im=x[i+j*L+k*N*L].im*z.re+x[i+j*L+k*N*L].re*z.im;
	      a.re=z.re*w.re-z.im*w.im;
	      a.im=z.re*w.im+z.im*w.re;
	      z.re=a.re;
	      z.im=a.im;
	    } 
	  for(k=0;k<r;k++)
	    {
	      z.re=z.im=0.0;
	      a.re=1.0;
	      a.im=0.0;
	      for(l=0;l<r;l++)
		{
		  z.re=z.re+T[l].re*a.re-T[l].im*a.im;
		  z.im=z.im+T[l].re*a.im+T[l].im*a.re;
		  b.re=a.re*R[k].re-a.im*R[k].im;
		  b.im=a.re*R[k].im+a.im*R[k].re;
		  a.re=b.re;
		  a.im=b.im;
		}
	      y[i+j*r*L+k*L].re=z.re;
	      y[i+j*r*L+k*L].im=z.im;
	    }
	  z.re=omega.re*w.re-omega.im*w.im+w.re;
	  z.im=omega.re*w.im+omega.im*w.re+w.im; 
	  w.re=z.re;
	  w.im=z.im;
	}
    }
}

