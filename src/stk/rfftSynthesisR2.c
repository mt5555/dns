/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftSynthesisR2.c,v 1.1 2001-08-06 20:33:34 mt Exp $
 */

#include "rfft.h"

void rfftSynthesisR2(real *y,real *x,complex omega,int L,int N)
{
  complex w,a,b;
  complex T0,T1;
  complex *x0,*x1,*y0,*y1;
  real *X0,*X1,*Y0,*Y1;
  int k,j;

  X0=x;
  X1=x+N*L;
  Y0=y;
  Y1=y+2*L-1;
  for(j=0;j<N;j++)
    {
      T0.re=*Y0;
      T1.re=*Y1;
      *X0=T0.re+T1.re;
      *X1=T0.re-T1.re;
      w.re=omega.re+1.0;                   
      w.im=omega.im;
      X0++;
      X1++;
      Y0++;
      Y1--;
      Y1--;
      x0=(complex *)X0;
      x1=(complex *)X1;
      y0=(complex *)Y0;
      y1=(complex *)Y1;
      for(k=1;k<(L+1)/2;k++)
	{
	  T0.re= y0->re;
	  T0.im= y0->im;
	  T1.re= y1->re;
	  T1.im=-y1->im;
	  a.re=T0.re+T1.re;
	  a.im=T0.im+T1.im;
	  b.re=T0.re-T1.re;
	  b.im=T0.im-T1.im;
	  x0->re=a.re;
	  x0->im=a.im;
	  x1->re=b.re*w.re+b.im*w.im;
	  x1->im=b.im*w.re-b.re*w.im;
	  a.re=omega.re*w.re-omega.im*w.im+w.re;
	  a.im=omega.re*w.im+omega.im*w.re+w.im; 
	  w.re=a.re;
	  w.im=a.im;
	  x0++;
	  x1++;
	  y0++;
	  y1--;
	}
      X0=(real *)x0;
      X1=(real *)x1;
      if(L%2==0)
	{
	  T0.re=y0->re; 
	  T0.im=y0->im; 
	  T1.re= T0.re;
	  T1.im=-T0.im;
	  a.re=T0.re+T1.re;
	  a.im=T0.im+T1.im;
	  b.re=T0.re-T1.re;
	  b.im=T0.im-T1.im;
	  *X0=a.re;
	  *X1=b.re*w.re+b.im*w.im;
	  X0++;
	  X1++;
	  Y0=((real *)y0)+1+L;
	  Y1=((real *)y1)+3*L;
	}
      else
	{
	  Y0=((real *)y0)+L;
	  Y1=((real *)y1)+1+3*L;
	}
    }
}

