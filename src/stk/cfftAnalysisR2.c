/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: cfftAnalysisR2.c,v 1.1 2001-08-06 20:33:30 mt Exp $
 */

#include "cfft.h"

void cfftAnalysisR2(complex *x,complex *y,complex omega,int L,int N)
{
  register complex w;
  register complex t0,t1;
  complex *x0,*x1;
  complex *y0,*y1;
  int j,i;
  x0=x;
  x1=x+N*L;
  y0=y;
  y1=y+L;
  for(j=0;j<N;j++)
    {
      w.re=1.0;
      w.im=0.0;
      for(i=0;i<L;i++)
	{
	  t0.re=x0->re;
	  t0.im=x0->im;
	  t1.re=w.re*x1->re-w.im*x1->im;
	  t1.im=w.re*x1->im+w.im*x1->re;
	  y0->re=t0.re+t1.re;
	  y0->im=t0.im+t1.im;	
	  y1->re=t0.re-t1.re;
	  y1->im=t0.im-t1.im;	
	  t1.re=omega.re*w.re-omega.im*w.im+w.re;
	  t1.im=omega.re*w.im+omega.im*w.re+w.im; 
	  w.re=t1.re;
	  w.im=t1.im;
	  x0++;
	  x1++;
	  y0++;
	  y1++;
	}
      y0=y0+L;
      y1=y1+L;
    }
}	  

