/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftAnalysisRg.c,v 1.1 2001-08-06 20:33:34 mt Exp $
 */

#include "rfft.h"

void rfftAnalysisRg(real *x,real *y,complex *R,complex *T,
		complex omega,int L,int N,int r)
{
  complex w,z,a,b;
  complex c,d,Z;
  real s,t;
  int l,k,j,q,K;

  for(j=0;j<N;j++)
    {
      for(q=0;q<r;q++)
	{
	  T[q].re=x[q*N*L+j*L];
	}
      s=T[0].re;           
      for(q=1;q<r;q++)
	{
	  s=s+T[q].re;
	}
      y[0+j*r*L]=s;
      for(l=1;l<(r+1)/2;l++)
	{
	  z.re=T[0].re;
	  z.im=0.0;
	  w.re=R[l].re;
	  w.im=R[l].im;
	  for(q=1;q<r;q++)
	    {
	      z.re=z.re+T[q].re*w.re;
	      z.im=z.im+T[q].re*w.im;
	      a.re=w.re*R[l].re-w.im*R[l].im;
	      a.im=w.re*R[l].im+w.im*R[l].re;
	      w.re=a.re;
	      w.im=a.im;
	    }
	  K=l*L;
	  y[2*K-1+j*r*L]=z.re;
	  y[2*K  +j*r*L]=z.im;
	}
      if(r%2==0)
	{
	  t=T[0].re-T[1].re;
	  for(q=2;q<r-1;q+=2)
	    {
	      t=t+T[q].re-T[q+1].re;
	    }
	  K=(r/2)*L;
	  y[2*K-1+j*r*L]=t;
	}
      w.re=omega.re+1.0;                   
      w.im=omega.im;
      for(k=1;k<(L+1)/2;k++)
	{
	  z.re=1.0;
	  z.im=0.0;
	  for(q=0;q<r;q++)
	    {
	      b.re=x[2*k-1+q*N*L+j*L];
	      b.im=x[2*k  +q*N*L+j*L];
	      T[q].re=b.re*z.re-b.im*z.im;
	      T[q].im=b.im*z.re+b.re*z.im;
	      a.re=z.re*w.re-z.im*w.im;
	      a.im=z.re*w.im+z.im*w.re;
	      z.re=a.re;
	      z.im=a.im;
	    } 
	  z.re=T[0].re;              
	  z.im=T[0].im;
	  for(q=1;q<r;q++)
	    {
	      z.re=z.re+T[q].re;
	      z.im=z.im+T[q].im;
	    }
	  K=k;
	  y[2*K-1+j*r*L]=z.re;
	  y[2*K  +j*r*L]=z.im;
	  for(l=1;l<(r+1)/2;l++)
	    {
	      z.re=z.im=0.0;
	      Z.re=Z.im=0.0;
	      a.re=1.0;
	      a.im=0.0;
	      for(q=0;q<r;q++)
		{
		  c.re=T[q].re*a.re;
		  c.im=T[q].re*a.im;
		  d.re=T[q].im*a.re;
		  d.im=T[q].im*a.im;
		  z.re=z.re+c.re-d.im;
		  z.im=z.im+c.im+d.re;
		  Z.re=Z.re+c.re+d.im;
		  Z.im=Z.im-c.im+d.re;
		  b.re=a.re*R[l].re-a.im*R[l].im;
		  b.im=a.re*R[l].im+a.im*R[l].re;
		  a.re=b.re;
		  a.im=b.im;
		}
	      K=k+l*L;
	      y[2*K-1+j*r*L] = z.re;
	      y[2*K  +j*r*L] = z.im;
	      K=L-k+(l-1)*L;
	      y[2*K-1+j*r*L] = Z.re;
	      y[2*K  +j*r*L] =-Z.im;
	    }
	  if(r%2==0)
	    {
	      Z.re=0.0;
	      Z.im=0.0;
	      for(q=0;q<r-1;q+=2)
		{
		  Z.re=Z.re+T[q].re-T[q+1].re;
		  Z.im=Z.im+T[q].im-T[q+1].im;
		}
	      K=L-k+(r/2-1)*L;
	      y[2*K-1+j*r*L] = Z.re;
	      y[2*K  +j*r*L] =-Z.im;
	    }

	  z.re=omega.re*w.re-omega.im*w.im+w.re;
	  z.im=omega.re*w.im+omega.im*w.re+w.im; 
	  w.re=z.re;
	  w.im=z.im;
	}
      if(L%2==0)
	{
	  z.re=1.0;
	  z.im=0.0;
	  for(q=0;q<r;q++)
	    {
	      T[q].re=x[L-1+q*N*L+j*L]*z.re;
	      T[q].im=x[L-1+q*N*L+j*L]*z.im;
	      a.re=z.re*w.re-z.im*w.im;
	      a.im=z.re*w.im+z.im*w.re;
	      z.re=a.re;
	      z.im=a.im;
	    } 
	  z.re=T[0].re;
	  z.im=T[0].im;
	  for(q=1;q<r;q++)
	    {
	      z.re=z.re+T[q].re;
	      z.im=z.im+T[q].im;
	    }
	  K=L/2;
	  y[2*K-1+j*r*L]=z.re;
	  y[2*K  +j*r*L]=z.im;
	  for(l=1;l<r/2;l++)
	    {
	      z.re=z.im=0.0;
	      a.re=1.0;
	      a.im=0.0;
	      for(q=0;q<r;q++)
		{
		  z.re=z.re+T[q].re*a.re-T[q].im*a.im;
		  z.im=z.im+T[q].re*a.im+T[q].im*a.re;
		  b.re=a.re*R[l].re-a.im*R[l].im;
		  b.im=a.re*R[l].im+a.im*R[l].re;
		  a.re=b.re;
		  a.im=b.im;
		}
	      K=L/2+l*L;
	      y[2*K-1+j*r*L] = z.re;
	      y[2*K  +j*r*L] = z.im;
	    }
	  if(r%2==1)
	    {
	      s=0.0;
	      a.re=1.0;
	      a.im=0.0;
	      l=(r-1)/2;
	      for(q=0;q<r;q++)
		{
		  s=s+T[q].re*a.re-T[q].im*a.im;
		  b.re=a.re*R[l].re-a.im*R[l].im;
		  b.im=a.re*R[l].im+a.im*R[l].re;
		  a.re=b.re;
		  a.im=b.im;
		}
	      K=L/2+((r-1)/2)*L;
	      y[2*K-1+j*r*L] = s;
	    }
	}
    }
}

