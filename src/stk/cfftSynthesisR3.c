/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: cfftSynthesisR3.c,v 1.1 2001-08-06 20:33:32 mt Exp $
 */

#include "cfft.h"

void cfftSynthesisR3(complex *x,complex *y,complex omega,int L,int N)
{
  register complex w,w2;
  register complex t0,t1,t2;
  register complex z0,z1;
  complex a={-0.5L,0.5L*STK_R3};
  complex *x0,*x1,*x2;
  complex *y0,*y1,*y2;
  int j,i;
  x0=x;
  x1=x0+N*L;
  x2=x1+N*L;
  for(j=0;j<N;j++)
    {
      y0=y+3*L*j;
      y1=y0+L;
      y2=y1+L;
      w.re=1.0;
      w.im=0.0;
      for(i=0;i<L;i++)
	{
	  w2.re=w.re*w.re-w.im*w.im;
	  w2.im=2.0*w.re*w.im;
	  t0.re=x0->re;
	  t0.im=x0->im;
	  t1.re=w.re*x1->re-w.im*x1->im;
	  t1.im=w.re*x1->im+w.im*x1->re;
	  t2.re=w2.re*x2->re-w2.im*x2->im;
	  t2.im=w2.re*x2->im+w2.im*x2->re;
	  z0.re=t0.re+a.re*(t1.re+t2.re);
	  z0.im=a.im*(t2.im-t1.im);
	  z1.re=t0.im+a.re*(t1.im+t2.im);
	  z1.im=a.im*(t1.re-t2.re);
	  y0->re=t0.re+t1.re+t2.re;
	  y0->im=t0.im+t1.im+t2.im;
	  y1->re=z0.re+z0.im;
	  y1->im=z1.re+z1.im;
	  y2->re=z0.re-z0.im;
	  y2->im=z1.re-z1.im;
	  t1.re=omega.re*w.re+omega.im*w.im+w.re;
	  t1.im=omega.re*w.im-omega.im*w.re+w.im; 
	  w.re=t1.re;
	  w.im=t1.im;
	  x0++;
	  x1++;
	  x2++;
	  y0++;
	  y1++;
	  y2++;
	}
    }
}	  
