/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftSynthesisR8.c,v 1.1 2001-08-06 20:33:34 mt Exp $
 */

#include "rfft.h"

void rfftSynthesisR8(real *y,real *x,complex omega,int L,int N)
{
  complex w1,w2,w3,w4,w5,w6,w7;
  complex s0,s1,s2,s3,s4,s5,s6,s7;
  complex t0,t1,t2,t3,t4,t5,t6,t7;
  complex *x0,*x1,*x2,*x3,*x4,*x5,*x6,*x7;
  complex *y0,*y1,*y2,*y3,*y4,*y5,*y6,*y7;
  real *X0,*X1,*X2,*X3,*X4,*X5,*X6,*X7;
  real *Y0,*Y4;
  int NL,L8;
  real c=1.0/STK_R2;
  int j,k;
  NL=N*L;
  L8=L*8;

  X0=x;
  X1=X0+NL;
  X2=X1+NL;
  X3=X2+NL;
  X4=X3+NL;  
  X5=X4+NL;
  X6=X5+NL;
  X7=X6+NL;  
  for(j=0;j<N;j++)
    {
    
      Y0=y+j*L8;
      Y4=Y0+L8-1;
      y1=(complex *)(Y0-1)+L;
      y2=y1+L;
      y3=y2+L;
      y4=y1-1;
      y5=y4+L;
      y6=y5+L;
      y7=y6+L;
      
      t0.re=*Y0;
      t1.re=y1->re;
      t1.im=y1->im;
      t2.re=y2->re;
      t2.im=y2->im;
      t3.re=y3->re;
      t3.im=y3->im;
      t4.re=*Y4;
      
      s0.re=t0.re+2.0*(t1.re+t2.re+t3.re)+t4.re;
      s1.re=t0.re+2.0*(c*(t1.re-t1.im)-t2.im-c*(t3.re+t3.im))-t4.re;
      s2.re=t0.re+2.0*(-t1.im-t2.re+t3.im)+t4.re;
      s3.re=t0.re+2.0*(c*(-t1.re-t1.im)+t2.im+c*(t3.re-t3.im))-t4.re;
      s4.re=t0.re+2.0*(-t1.re+t2.re-t3.re)+t4.re;
      s5.re=t0.re+2.0*(c*(-t1.re+t1.im)-t2.im+c*(t3.re+t3.im))-t4.re;
      s6.re=t0.re+2.0*(t1.im-t2.re-t3.im)+t4.re;
      s7.re=t0.re+2.0*(c*(t1.re+t1.im)+t2.im+c*(-t3.re+t3.im))-t4.re;

      *X0=s0.re;
      *X1=s1.re;
      *X2=s2.re;
      *X3=s3.re;
      *X4=s4.re;
      *X5=s5.re;
      *X6=s6.re;
      *X7=s7.re;

      w1.re=omega.re+1.0;                   
      w1.im=omega.im;
      
      X0++;
      X1++;
      X2++;
      X3++;
      X4++;
      X5++;
      X6++;
      X7++;
      Y0++;
      y1++;
      y2++;
      y3++;
      x0=(complex *)X0;
      x1=(complex *)X1;
      x2=(complex *)X2;
      x3=(complex *)X3;
      x4=(complex *)X4;
      x5=(complex *)X5;
      x6=(complex *)X6;
      x7=(complex *)X7;
      y0=(complex *)Y0;

      for(k=1;k<(L+1)/2;k++)
	{
	  w2.re=w1.re*w1.re-w1.im*w1.im;
	  w2.im=2.0*w1.re*w1.im;
	  w3.re=w2.re*w1.re-w2.im*w1.im;
	  w3.im=w2.re*w1.im+w2.im*w1.re; 
	  w4.re=w3.re*w1.re-w3.im*w1.im;
	  w4.im=w3.re*w1.im+w3.im*w1.re; 
	  w5.re=w4.re*w1.re-w4.im*w1.im;
	  w5.im=w4.re*w1.im+w4.im*w1.re; 
	  w6.re=w5.re*w1.re-w5.im*w1.im;
	  w6.im=w5.re*w1.im+w5.im*w1.re; 
	  w7.re=w6.re*w1.re-w6.im*w1.im;
	  w7.im=w6.re*w1.im+w6.im*w1.re; 

	  t0.re= y0->re;
	  t0.im= y0->im;
	  t1.re= y1->re;
	  t1.im= y1->im;
	  t2.re= y2->re;
	  t2.im= y2->im;
	  t3.re= y3->re;
	  t3.im= y3->im;
	  t4.re= y4->re;
	  t4.im=-y4->im;
	  t5.re= y5->re;
	  t5.im=-y5->im;
	  t6.re= y6->re;
	  t6.im=-y6->im;
	  t7.re= y7->re;
	  t7.im=-y7->im;

	  s0.re=  ( t0.re+t2.re+t5.re+t7.re);
	  s0.im=  ( t0.im+t2.im+t5.im+t7.im);
	  s1.re=  ( t1.re+t3.re+t4.re+t6.re);
	  s1.im=  ( t1.im+t3.im+t4.im+t6.im);

	  s2.re=  ( t0.re-t2.im+t5.im-t7.re);
	  s2.im=  ( t0.im+t2.re-t5.re-t7.im);
	  s3.re=c*( t1.re-t1.im-t3.re-t3.im+t4.re+t4.im-t6.re+t6.im);
	  s3.im=c*( t1.im+t1.re-t3.im+t3.re+t4.im-t4.re-t6.im-t6.re);

	  s4.re=  ( t0.re-t2.re-t5.re+t7.re);
	  s4.im=  ( t0.im-t2.im-t5.im+t7.im);
	  s5.re=  (-t1.im+t3.im+t4.im-t6.im);
	  s5.im=  ( t1.re-t3.re-t4.re+t6.re);

	  s6.re=  ( t0.re+t2.im-t5.im-t7.re);
	  s6.im=  ( t0.im-t2.re+t5.re-t7.im);
	  s7.re=c*(-t1.re-t1.im+t3.re-t3.im-t4.re+t4.im+t6.re+t6.im);
	  s7.im=c*(-t1.im+t1.re+t3.im+t3.re-t4.im-t4.re+t6.im-t6.re);

	  t0.re=s0.re+s1.re;
	  t0.im=s0.im+s1.im;
	  t1.re=s2.re+s3.re;
	  t1.im=s2.im+s3.im;
	  t2.re=s4.re+s5.re;
	  t2.im=s4.im+s5.im;
	  t3.re=s6.re+s7.re;
	  t3.im=s6.im+s7.im;
	  t4.re=s0.re-s1.re;
	  t4.im=s0.im-s1.im;
	  t5.re=s2.re-s3.re;
	  t5.im=s2.im-s3.im;
	  t6.re=s4.re-s5.re;
	  t6.im=s4.im-s5.im;
	  t7.re=s6.re-s7.re;
	  t7.im=s6.im-s7.im;

	  x0->re=t0.re;
	  x0->im=t0.im;
	  x1->re=t1.re*w1.re+t1.im*w1.im;
	  x1->im=t1.im*w1.re-t1.re*w1.im;
	  x2->re=t2.re*w2.re+t2.im*w2.im;
	  x2->im=t2.im*w2.re-t2.re*w2.im;
	  x3->re=t3.re*w3.re+t3.im*w3.im;
	  x3->im=t3.im*w3.re-t3.re*w3.im;
	  x4->re=t4.re*w4.re+t4.im*w4.im;
	  x4->im=t4.im*w4.re-t4.re*w4.im;
	  x5->re=t5.re*w5.re+t5.im*w5.im;
	  x5->im=t5.im*w5.re-t5.re*w5.im;
	  x6->re=t6.re*w6.re+t6.im*w6.im;
	  x6->im=t6.im*w6.re-t6.re*w6.im;
	  x7->re=t7.re*w7.re+t7.im*w7.im;
	  x7->im=t7.im*w7.re-t7.re*w7.im;

	  t0.re=omega.re*w1.re-omega.im*w1.im+w1.re;
	  t0.im=omega.re*w1.im+omega.im*w1.re+w1.im; 
	  w1.re=t0.re;
	  w1.im=t0.im;
	  x0++;
	  x1++;
	  x2++;
	  x3++;
	  x4++;
	  x5++;
	  x6++;
	  x7++;
	  y0++;
	  y1++;
	  y2++;
	  y3++;
	  y4--;
	  y5--;
	  y6--;
	  y7--;
	}
      X0=(real *)x0;
	X1=(real *)x1;
	X2=(real *)x2;
	X3=(real *)x3;
	X4=(real *)x4;
	X5=(real *)x5;
	X6=(real *)x6;
	X7=(real *)x7;
      if(L%2==0)
	{
	  w2.re=w1.re*w1.re-w1.im*w1.im;
	  w2.im=2.0*w1.re*w1.im;
	  w3.re=w2.re*w1.re-w2.im*w1.im;
	  w3.im=w2.re*w1.im+w2.im*w1.re; 
	  w4.re=w3.re*w1.re-w3.im*w1.im;
	  w4.im=w3.re*w1.im+w3.im*w1.re; 
	  w5.re=w4.re*w1.re-w4.im*w1.im;
	  w5.im=w4.re*w1.im+w4.im*w1.re; 
	  w6.re=w5.re*w1.re-w5.im*w1.im;
	  w6.im=w5.re*w1.im+w5.im*w1.re; 
	  w7.re=w6.re*w1.re-w6.im*w1.im;
	  w7.im=w6.re*w1.im+w6.im*w1.re; 
	  	
	  t0.re=y0->re;
	  t0.im=y0->im;
	  t1.re=y1->re;
	  t1.im=y1->im;
	  t2.re=y2->re;
	  t2.im=y2->im;
	  t3.re=y3->re;
	  t3.im=y3->im;
	  t4.re= t3.re;
	  t4.im=-t3.im;
	  t5.re= t2.re;
	  t5.im=-t2.im;
	  t6.re= t1.re;
	  t6.im=-t1.im;	      	      
	  t7.re= t0.re;
	  t7.im=-t0.im;	  

	  s0.re= t0.re+t2.re+t4.re+t6.re;
	  s0.im= t0.im+t2.im+t4.im+t6.im;
	  s1.re= t1.re+t3.re+t5.re+t7.re;
	  s1.im= t1.im+t3.im+t5.im+t7.im;
	  
	  s2.re= t0.re-t2.im-t4.re+t6.im;
	  s2.im= t0.im+t2.re-t4.im-t6.re;
	  s3.re=c*(t1.re-t1.im-t3.re-t3.im-t5.re+t5.im+t7.re+t7.im);
	  s3.im=c*(t1.im+t1.re-t3.im+t3.re-t5.im-t5.re+t7.im-t7.re);

	  s4.re= t0.re-t2.re+t4.re-t6.re;
	  s4.im= t0.im-t2.im+t4.im-t6.im;
	  s5.re=-t1.im+t3.im-t5.im+t7.im;
	  s5.im= t1.re-t3.re+t5.re-t7.re;

	  s6.re= t0.re+t2.im-t4.re-t6.im;
	  s6.im= t0.im-t2.re-t4.im+t6.re;
	  s7.re=c*(-t1.re-t1.im+t3.re-t3.im+t5.re+t5.im-t7.re+t7.im);
	  s7.im=c*(-t1.im+t1.re+t3.im+t3.re+t5.im-t5.re-t7.im-t7.re);

	  t0.re=s0.re+s1.re;
	  t0.im=s0.im+s1.im;
	  t1.re=s2.re+s3.re;
	  t1.im=s2.im+s3.im;
	  t2.re=s4.re+s5.re;
	  t2.im=s4.im+s5.im;
	  t3.re=s6.re+s7.re;
	  t3.im=s6.im+s7.im;
	  t4.re=s0.re-s1.re;
	  t4.im=s0.im-s1.im;
	  t5.re=s2.re-s3.re;
	  t5.im=s2.im-s3.im;
	  t6.re=s4.re-s5.re;
	  t6.im=s4.im-s5.im;
	  t7.re=s6.re-s7.re;
	  t7.im=s6.im-s7.im;

	  *X0=t0.re;
	  *X1=t1.re*w1.re+t1.im*w1.im;
	  *X2=t2.re*w2.re+t2.im*w2.im;
	  *X3=t3.re*w3.re+t3.im*w3.im;	  
	  *X4=t4.re*w4.re+t4.im*w4.im;	  
	  *X5=t5.re*w5.re+t5.im*w5.im;	  
	  *X6=t6.re*w6.re+t6.im*w6.im;	  	   
	  *X7=t7.re*w7.re+t7.im*w7.im;	  
	  X0++;
	  X1++;
	  X2++;
	  X3++;
	  X4++;
	  X5++;
	  X6++;
	  X7++;
		  
	}
    }
}
