/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftSynthesisR4.c,v 1.1 2001-08-06 20:33:34 mt Exp $
 */

#include "rfft.h"

void rfftSynthesisR4(real *y,real *x,complex omega,int L,int N)
{
  register complex t0,t1,t2,t3;
  register complex w,w2,w3;
  register complex z0,z1,z2,z3;
  complex *x0,*x1,*x2,*x3,*y0,*y1,*y2,*y3;
  real *X0,*X1,*X2,*X3,*Y0,*Y2;
  int NL,L4;
  int k,j;
  NL=N*L;
  L4=L*4;
  X0=x;
  X1=X0+NL;
  X2=X1+NL;
  X3=X2+NL;
  for(j=0;j<N;j++)
    {
      Y0=y+j*L4;
      Y2=Y0+L4-1;
      y1=(complex *)(Y0-1)+L;
      y2=y1-1;
      y3=y2+L;
      t0.re=*Y0;
      t1.re=y1->re;
      t1.im=y1->im;
      t2.re=*Y2;
      *X0=t0.re+2.0*t1.re+t2.re;
      *X1=t0.re-2.0*t1.im-t2.re;
      *X2=t0.re-2.0*t1.re+t2.re;
      *X3=t0.re+2.0*t1.im-t2.re;
      w.re=omega.re+1.0;                   
      w.im=omega.im;
      X0++;
      X1++;
      X2++;
      X3++;
      Y0++;
      y1++;
      x0=(complex *)X0;
      x1=(complex *)X1;
      x2=(complex *)X2;
      x3=(complex *)X3;
      y0=(complex *)Y0;
      for(k=1;k<(L+1)/2;k++)
	{	      
	  w2.re=w.re*w.re-w.im*w.im;
	  w2.im=2.0*w.re*w.im;
	  w3.re=w.re*w2.re-w.im*w2.im;
	  w3.im=w.re*w2.im+w.im*w2.re;
	  t0.re= y0->re; 
	  t0.im= y0->im;
	  t1.re= y1->re;
	  t1.im= y1->im; 
	  t2.re= y2->re; 
	  t2.im=-y2->im; 
	  t3.re= y3->re; 
	  t3.im=-y3->im; 
	  z0.re=t0.re+t1.re+t2.re+t3.re;
	  z0.im=t0.im+t1.im+t2.im+t3.im;
	  z1.re=t0.re-t1.im+t2.im-t3.re;
	  z1.im=t0.im+t1.re-t2.re-t3.im;
	  z2.re=t0.re-t1.re-t2.re+t3.re;
	  z2.im=t0.im-t1.im-t2.im+t3.im;
	  z3.re=t0.re+t1.im-t2.im-t3.re;
	  z3.im=t0.im-t1.re+t2.re-t3.im;
	  x0->re=z0.re;
	  x0->im=z0.im;
	  x1->re=z1.re*w.re+z1.im*w.im;
	  x1->im=z1.im*w.re-z1.re*w.im;
	  x2->re=z2.re*w2.re+z2.im*w2.im;
	  x2->im=z2.im*w2.re-z2.re*w2.im;
	  x3->re=z3.re*w3.re+z3.im*w3.im;
	  x3->im=z3.im*w3.re-z3.re*w3.im;
	  z0.re=omega.re*w.re-omega.im*w.im+w.re;
	  z0.im=omega.re*w.im+omega.im*w.re+w.im; 
	  w.re=z0.re;
	  w.im=z0.im;
	  x0++;
	  x1++;
	  x2++;
	  x3++;
	  y0++;
	  y1++;
	  y2--;
	  y3--;  
	}
	X0=(real *)x0;
      X1=(real *)x1;
      X2=(real *)x2;
      X3=(real *)x3;
      if(L%2==0)
	{
	  w2.re=w.re*w.re-w.im*w.im;
	  w2.im=2.0*w.re*w.im;
	  w3.re=w.re*w2.re-w.im*w2.im;
	  w3.im=w.re*w2.im+w.im*w2.re;
	  t0.re=y0->re;
	  t0.im=y0->im;
	  t1.re=y1->re;
	  t1.im=y1->im;
	  t2.re= t1.re;
	  t2.im=-t1.im;
	  t3.re= t0.re;
	  t3.im=-t0.im;         
	  z0.re=t0.re+t1.re+t2.re+t3.re;
	  z0.im=t0.im+t1.im+t2.im+t3.im;
	  z1.re=t0.re-t1.im-t2.re+t3.im;
	  z1.im=t0.im+t1.re-t2.im-t3.re;
	  z2.re=t0.re-t1.re+t2.re-t3.re;
	  z2.im=t0.im-t1.im+t2.im-t3.im;
	  z3.re=t0.re+t1.im-t2.re-t3.im;
	  z3.im=t0.im-t1.re-t2.im+t3.re;
	  *X0=z0.re;
	  *X1=z1.re*w.re+z1.im*w.im;
	  *X2=z2.re*w2.re+z2.im*w2.im;
	  *X3=z3.re*w3.re+z3.im*w3.im;
	  X0++;
	  X1++;
	  X2++;
	  X3++;
	}
    }
}


