/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: cfftAnalysisR5.c,v 1.1 2001-08-06 20:33:31 mt Exp $
 */

#include "cfft.h"

void cfftAnalysisR5(complex *x,complex *y,complex omega,int L,int N)
{
  register complex w1,w2,w3,w4;
  register complex t0,t1,t2,t3,t4;
  register complex z0,z1;
  complex *x0,*x1,*x2,*x3,*x4;
  complex *y0,*y1,*y2,*y3,*y4;
  complex a1={STK_R5C2,STK_R5S2};
  complex a2={STK_R5C4,STK_R5S4};
  int j,i;
  x0=x;
  x1=x0+L*N;
  x2=x1+L*N;
  x3=x2+L*N;
  x4=x3+L*N;
  for(j=0;j<N;j++)
    { 
      y0=y+5*L*j;
      y1=y0+L;
      y2=y1+L;
      y3=y2+L;
      y4=y3+L;
      t0.re=x0->re;
      t0.im=x0->im;
      t1.re=x1->re;
      t1.im=x1->im;
      t2.re=x2->re;
      t2.im=x2->im;
      t3.re=x3->re;
      t3.im=x3->im;
      t4.re=x4->re;
      t4.im=x4->im;
      y0->re=t0.re+t1.re+t2.re+t3.re+t4.re;
      y0->im=t0.im+t1.im+t2.im+t3.im+t4.im;
      z0.re=t0.re+a1.re*(t1.re+t4.re)+a2.re*(t2.re+t3.re);
      z0.im=      a1.im*(t4.im-t1.im)+a2.im*(t3.im-t2.im);
      z1.re=t0.im+a1.re*(t1.im+t4.im)+a2.re*(t2.im+t3.im);
      z1.im=      a1.im*(t1.re-t4.re)+a2.im*(t2.re-t3.re);
      y1->re=z0.re+z0.im;
      y4->re=z0.re-z0.im;
      y1->im=z1.re+z1.im;
      y4->im=z1.re-z1.im;
      z0.re=t0.re+a2.re*(t1.re+t4.re)+a1.re*(t2.re+t3.re);
      z0.im=      a2.im*(t4.im-t1.im)+a1.im*(t2.im-t3.im);
      z1.re=t0.im+a2.re*(t1.im+t4.im)+a1.re*(t2.im+t3.im);
      z1.im=     +a2.im*(t1.re-t4.re)+a1.im*(t3.re-t2.re);
      y2->re=z0.re+z0.im;
      y3->re=z0.re-z0.im;
      y2->im=z1.re+z1.im;
      y3->im=z1.re-z1.im;
      w1.re=omega.re+1.0;
      w1.im=omega.im;
      x0++;
      x1++;
      x2++;
      x3++;
      x4++;
      y0++;
      y1++;
      y2++;
      y3++;
      y4++;
      for(i=1;i<L;i++)
	{
	  w2.re = w1.re*w1.re-w1.im*w1.im;
	  w2.im = 2.0*w1.re*w1.im;
	  w3.re = w2.re*w1.re-w2.im*w1.im;
	  w3.im = w2.re*w1.im+w2.im*w1.re; 
	  w4.re = w3.re*w1.re-w3.im*w1.im;
	  w4.im = w3.re*w1.im+w3.im*w1.re; 
	  t0.re=x0->re;
	  t0.im=x0->im;
	  t1.re = w1.re*x1->re-w1.im*x1->im;
	  t1.im = w1.re*x1->im+w1.im*x1->re;
	  t2.re = w2.re*x2->re-w2.im*x2->im;
	  t2.im = w2.re*x2->im+w2.im*x2->re;
	  t3.re = w3.re*x3->re-w3.im*x3->im;
	  t3.im = w3.re*x3->im+w3.im*x3->re;
	  t4.re = w4.re*x4->re-w4.im*x4->im;
	  t4.im = w4.re*x4->im+w4.im*x4->re;
	  y0->re=t0.re+t1.re+t2.re+t3.re+t4.re;
	  y0->im=t0.im+t1.im+t2.im+t3.im+t4.im;
	  z0.re=t0.re+a1.re*(t1.re+t4.re)+a2.re*(t2.re+t3.re);
	  z0.im=      a1.im*(t4.im-t1.im)+a2.im*(t3.im-t2.im);
	  z1.re=t0.im+a1.re*(t1.im+t4.im)+a2.re*(t2.im+t3.im);
	  z1.im=      a1.im*(t1.re-t4.re)+a2.im*(t2.re-t3.re);
	  y1->re=z0.re+z0.im;
	  y4->re=z0.re-z0.im;
	  y1->im=z1.re+z1.im;
	  y4->im=z1.re-z1.im;
	  z0.re=t0.re+a2.re*(t1.re+t4.re)+a1.re*(t2.re+t3.re);
	  z0.im=      a2.im*(t4.im-t1.im)+a1.im*(t2.im-t3.im);
	  z1.re=t0.im+a2.re*(t1.im+t4.im)+a1.re*(t2.im+t3.im);
	  z1.im=     +a2.im*(t1.re-t4.re)+a1.im*(t3.re-t2.re);
	  y2->re=z0.re+z0.im;
	  y3->re=z0.re-z0.im;
	  y2->im=z1.re+z1.im;
	  y3->im=z1.re-z1.im;
	  t1.re=omega.re*w1.re-w1.im*omega.im+w1.re;
	  t1.im=omega.re*w1.im+w1.re*omega.im+w1.im;
	  w1.re=t1.re;
	  w1.im=t1.im;
	  x0++;
	  x1++;
	  x2++;
	  x3++;
	  x4++;
	  y0++;
	  y1++;
	  y2++;
	  y3++;
	  y4++;
	}
    }
}	  
