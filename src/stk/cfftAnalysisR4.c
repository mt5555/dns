/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: cfftAnalysisR4.c,v 1.1 2001-08-06 20:33:31 mt Exp $
 */

#include "cfft.h"

void cfftAnalysisR4(complex *x,complex *y,complex omega,int L,int N)
{
  register complex w1,w2,w3;
  register complex t0,t1,t2,t3;
  complex *x0,*x1,*x2,*x3;
  complex *y0,*y1,*y2,*y3;
  int j,i;
  x0=x;
  x1=x0+L*N;
  x2=x1+L*N;
  x3=x2+L*N;
  for(j=0;j<N;j++)
    { 
      y0=y+4*L*j;
      y1=y0+L;
      y2=y1+L;
      y3=y2+L;
      t0.re=x0->re;
      t0.im=x0->im;
      t1.re=x1->re;
      t1.im=x1->im;
      t2.re=x2->re;
      t2.im=x2->im;
      t3.re=x3->re;
      t3.im=x3->im;
      y0->re = t0.re+t1.re+t2.re+t3.re;
      y0->im = t0.im+t1.im+t2.im+t3.im;
      y1->re = t0.re+t1.im-t2.re-t3.im;
      y1->im = t0.im-t1.re-t2.im+t3.re;
      y2->re = t0.re-t1.re+t2.re-t3.re;
      y2->im = t0.im-t1.im+t2.im-t3.im;
      y3->re = t0.re-t1.im-t2.re+t3.im;
      y3->im = t0.im+t1.re-t2.im-t3.re;
      w1.re=omega.re+1.0;
      w1.im=omega.im;
      x0++;
      x1++;
      x2++;
      x3++;
      y0++;
      y1++;
      y2++;
      y3++;
      for(i=1;i<L;i++)
	{
	  w2.re = w1.re*w1.re-w1.im*w1.im;
	  w2.im = 2.0*w1.re*w1.im;
	  w3.re = w2.re*w1.re-w2.im*w1.im;
	  w3.im = w2.re*w1.im+w2.im*w1.re; 
	  t0.re=x0->re;
	  t0.im=x0->im;
	  t1.re = w1.re*x1->re-w1.im*x1->im;
	  t1.im = w1.re*x1->im+w1.im*x1->re;
	  t2.re = w2.re*x2->re-w2.im*x2->im;
	  t2.im = w2.re*x2->im+w2.im*x2->re;
	  t3.re = w3.re*x3->re-w3.im*x3->im;
	  t3.im = w3.re*x3->im+w3.im*x3->re;
	  y0->re = t0.re+t1.re+t2.re+t3.re;
	  y0->im = t0.im+t1.im+t2.im+t3.im;
	  y1->re = t0.re+t1.im-t2.re-t3.im;
	  y1->im = t0.im-t1.re-t2.im+t3.re;
	  y2->re = t0.re-t1.re+t2.re-t3.re;
	  y2->im = t0.im-t1.im+t2.im-t3.im;
	  y3->re = t0.re-t1.im-t2.re+t3.im;
	  y3->im = t0.im+t1.re-t2.im-t3.re;
	  t1.re=omega.re*w1.re-w1.im*omega.im+w1.re;
	  t1.im=omega.re*w1.im+w1.re*omega.im+w1.im;
	  w1.re=t1.re;
	  w1.im=t1.im;
	  x0++;
	  x1++;
	  x2++;
	  x3++;
	  y0++;
	  y1++;
	  y2++;
	  y3++;
	}
    }
}	  









