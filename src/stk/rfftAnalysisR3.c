/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftAnalysisR3.c,v 1.1 2001-08-06 20:33:33 mt Exp $
 */

#include "rfft.h"

void rfftAnalysisR3(real *x,real *y,complex omega,int L,int N)
{
  register complex t0,t1,t2;
  register complex w1,w2;
  register complex z1,z2;
  complex *x0,*x1,*x2,*y0,*y1,*y2;
  real *X0,*X1,*X2,*Y0,*Y1;
  int k,j;
  int NL,L3,L2;
  complex a={-0.5L,0.5L*STK_R3};
  NL=N*L;
  L3=L*3;
  L2=L*2;
  X0=x;
  X1=X0+NL;
  X2=X1+NL;
  for(j=0;j<N;j++)
    {
      Y0=y+j*L3;
      Y1=Y0+L2-1;
      y1=(complex *)Y1;
      y2=y1-1;
      t0.re=*X0; 
      t1.re=*X1;
      t2.re=*X2;
      z1.re=t1.re+t2.re;
      z1.im=t2.re-t1.re;
      *Y0=t0.re+z1.re;
      y1->re=t0.re+a.re*z1.re;
      y1->im=      a.im*z1.im;
      w1.re=omega.re+1.0;
      w1.im=omega.im;
      X0++;
      X1++;
      X2++;
      Y0++;
      y1++;
      x0=(complex *)X0;
      x1=(complex *)X1;
      x2=(complex *)X2;
      y0=(complex *)Y0;
      for(k=1;k<(L+1)/2;k++)
	{
	  w2.re=w1.re*w1.re-w1.im*w1.im;
	  w2.im=2.0*w1.re*w1.im;
	  t0.re=x0->re;
	  t0.im=x0->im;
	  z1.re=x1->re;
	  z1.im=x1->im;
	  z2.re=x2->re;
	  z2.im=x2->im;
	  t1.re=z1.re*w1.re-z1.im*w1.im;
	  t1.im=z1.im*w1.re+z1.re*w1.im;
	  t2.re=z2.re*w2.re-z2.im*w2.im;
	  t2.im=z2.im*w2.re+z2.re*w2.im;
	  z1.re=a.re*(t1.re+t2.re);
	  z1.im=a.re*(t1.im+t2.im);
	  z2.re=a.im*(t1.re-t2.re);
	  z2.im=a.im*(t1.im-t2.im);
	  y0->re= (t0.re+t1.re+t2.re);
	  y0->im= (t0.im+t1.im+t2.im);
	  y1->re= (t0.re+z1.re+z2.im);
	  y1->im= (t0.im+z1.im-z2.re);
	  y2->re= (t0.re+z1.re-z2.im);
	  y2->im=-(t0.im+z1.im+z2.re);
	  w2.re=omega.re*w1.re-omega.im*w1.im+w1.re;
	  w2.im=omega.re*w1.im+omega.im*w1.re+w1.im; 
	  w1.re=w2.re;
	  w1.im=w2.im;
	  x0++;
	  x1++;
	  x2++;
	  y0++;
	  y1++;
	  y2--;
	}
      X0=(real *)x0;
      X1=(real *)x1;
      X2=(real *)x2;
      Y1=(real *)y1;
      if(L%2==0)
	{
	  w2.re=w1.re*w1.re-w1.im*w1.im;
	  w2.im=2.0*w1.re*w1.im;
	  t0.re=*X0;
	  t0.im=0.0;
	  t1.re=*X1*w1.re; 
	  t1.im=*X1*w1.im;
	  t2.re=*X2*w2.re; 
	  t2.im=*X2*w2.im;
	  y0->re=t0.re+t1.re+t2.re; 
	  y0->im=t0.im+t1.im+t2.im; 
	  *Y1=t0.re+a.re*(t1.re+t2.re)+a.im*(t1.im-t2.im);
	  X0++;
	  X1++;
	  X2++;
	}
    }
}

