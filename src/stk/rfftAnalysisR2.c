/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftAnalysisR2.c,v 1.1 2001-08-06 20:33:33 mt Exp $
 */

#include "rfft.h"

void rfftAnalysisR2(real *x,real *y,complex omega,int L,int N)
{
  register complex w,T0,T1;
  complex *x0,*x1,*y0,*y1;
  real *X0,*X1,*Y0,*Y1;
  int j,k;

  X0=x;
  X1=x+N*L;
  Y0=y;
  Y1=y+2*L-1;
  for(j=0;j<N;j++)
    {
      T0.re=*X0;
      T1.re=*X1;
      *Y0=T0.re+T1.re;
      *Y1=T0.re-T1.re;
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
	  T0.re=x0->re;
	  T0.im=x0->im;
	  T1.re=x1->re*w.re-x1->im*w.im;
	  T1.im=x1->im*w.re+x1->re*w.im;
	  y0->re=T0.re+T1.re;              
	  y0->im=T0.im+T1.im;
	  y1->re=T0.re-T1.re;
	  y1->im=T1.im-T0.im;
	  T0.re=omega.re*w.re-omega.im*w.im+w.re;
	  T0.im=omega.re*w.im+omega.im*w.re+w.im; 
	  w.re=T0.re;
	  w.im=T0.im;
	  x0++;
	  x1++;
	  y0++;
	  y1--;
	}
      X0=(real *)x0;
      X1=(real *)x1;
      if(L%2==0)
	{
	  T0.re=*X0;
	  T1.re=*X1*w.re;
	  T1.im=*X1*w.im;
	  y0->re=T0.re+T1.re;
	  y0->im=T1.im;
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



