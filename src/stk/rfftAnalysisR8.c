/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftAnalysisR8.c,v 1.1 2001-08-06 20:33:34 mt Exp $
 */

#include "rfft.h"

void rfftAnalysisR8(real *x,real *y,complex omega,int L,int N)
{
  register complex w1,w2,w3,w4,w5,w6,w7;
  register complex s0,s1,s2,s3,s4,s5,s6,s7;
  register complex t0,t1,t2,t3,t4,t5,t6,t7;
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
      
      t0.re=*X0;
      t1.re=*X1; 
      t2.re=*X2; 
      t3.re=*X3; 
      t4.re=*X4; 
      t5.re=*X5; 
      t6.re=*X6; 
      t7.re=*X7; 
           
      s0.re=t0.re+t4.re;
      s1.re=t2.re+t6.re;
      s2.re=t1.re+t5.re;
      s3.re=t3.re+t7.re;
      s4.re=t0.re-t4.re;
      s5.re=t2.re-t6.re;
      s6.re=t1.re-t5.re;
      s7.re=t3.re-t7.re;    
      t0.re=  (s0.re+s1.re);
      t1.re=  (s2.re+s3.re);
      t2.re=  s4.re;
      t2.im= -s5.re;
      t3.re=c*(s6.re-s7.re);
      t3.im=c*(-s6.re-s7.re);
      t4.re=  (s0.re-s1.re); 
      t5.im= -(s2.re-s3.re);
      t6.re=  s4.re;
      t6.im=  s5.re;
      t7.re=-t3.re;
      t7.im=t3.im;
      *Y0   =t0.re+t1.re;
      y1->re=t2.re+t3.re;
      y1->im=t2.im+t3.im;
      y2->re=t4.re;
      y2->im=t5.im;
      y3->re=t6.re+t7.re;
      y3->im=t6.im+t7.im;
      *Y4   =t0.re-t1.re;
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
	  t0.re=x0->re; 
	  t0.im=x0->im; 
	  s1.re=x1->re; 
	  s1.im=x1->im; 
	  s2.re=x2->re; 
	  s2.im=x2->im; 	      
	  s3.re=x3->re;
	  s3.im=x3->im; 	      
	  s4.re=x4->re; 
	  s4.im=x4->im;      
	  s5.re=x5->re; 
	  s5.im=x5->im; 
	  s6.re=x6->re; 
	  s6.im=x6->im; 	      
	  s7.re=x7->re; 
	  s7.im=x7->im; 

	  t1.re=s1.re*w1.re-s1.im*w1.im;
	  t1.im=s1.im*w1.re+s1.re*w1.im;
	  t2.re=s2.re*w2.re-s2.im*w2.im;
	  t2.im=s2.im*w2.re+s2.re*w2.im;
	  t3.re=s3.re*w3.re-s3.im*w3.im;
	  t3.im=s3.im*w3.re+s3.re*w3.im;
	  t4.re=s4.re*w4.re-s4.im*w4.im;
	  t4.im=s4.im*w4.re+s4.re*w4.im;
	  t5.re=s5.re*w5.re-s5.im*w5.im;
	  t5.im=s5.im*w5.re+s5.re*w5.im;	      
	  t6.re=s6.re*w6.re-s6.im*w6.im;
	  t6.im=s6.im*w6.re+s6.re*w6.im;
	  t7.re=s7.re*w7.re-s7.im*w7.im;
	  t7.im=s7.im*w7.re+s7.re*w7.im;

        s0.re=t0.re+t4.re;
        s0.im=t0.im+t4.im;
        s1.re=t2.re+t6.re;
        s1.im=t2.im+t6.im;
        s2.re=t1.re+t5.re;
        s2.im=t1.im+t5.im;
        s3.re=t3.re+t7.re;
        s3.im=t3.im+t7.im;
        s4.re=t0.re-t4.re;
        s4.im=t0.im-t4.im;
        s5.re=t2.re-t6.re;
        s5.im=t2.im-t6.im;
        s6.re=t1.re-t5.re;
        s6.im=t1.im-t5.im;
        s7.re=t3.re-t7.re;
        s7.im=t3.im-t7.im;
	t0.re=  (s0.re+s1.re);
	t0.im=  (s0.im+s1.im);
	t1.re=  (s2.re+s3.re);
	t1.im=  (s2.im+s3.im);
	t2.re=  (s4.re+s5.im);
	t2.im=  (s4.im-s5.re);
	t3.re=c*(s6.re+s6.im-s7.re+s7.im);
	t3.im=c*(s6.im-s6.re-s7.im-s7.re);
	t4.re=  (s0.re-s1.re); 
	t4.im=  (s0.im-s1.im); 
	t5.re=  (s2.im-s3.im); 
	t5.im= -(s2.re-s3.re); 
	t6.re=  (s4.re-s5.im); 
	t6.im=  (s4.im+s5.re); 
	t7.re=c*(-s6.re+s6.im+s7.re+s7.im);
	t7.im=c*(-s6.im-s6.re+s7.im-s7.re);
	y0->re= (t0.re+t1.re);
	y0->im= (t0.im+t1.im);
	y1->re= (t2.re+t3.re);
	y1->im= (t2.im+t3.im);
	y2->re= (t4.re+t5.re);
	y2->im= (t4.im+t5.im);
	y3->re= (t6.re+t7.re);
	y3->im= (t6.im+t7.im);
	
	y4->re= (t6.re-t7.re);
	y4->im=-(t6.im-t7.im);
	y5->re= (t4.re-t5.re);
	y5->im=-(t4.im-t5.im);
	y6->re= (t2.re-t3.re);
	y6->im=-(t2.im-t3.im);
	y7->re= (t0.re-t1.re);
	y7->im=-(t0.im-t1.im);
	
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
	  t0.re=*X0;
	  t0.im=0.0;
	  t1.re=*X1*w1.re;
	  t1.im=*X1*w1.im;
	  t2.re=*X2*w2.re;
	  t2.im=*X2*w2.im;
	  t3.re=*X3*w3.re;
	  t3.im=*X3*w3.im;	
	  t4.re=*X4*w4.re;
	  t4.im=*X4*w4.im;
	  t5.re=*X5*w5.re;
	  t5.im=*X5*w5.im;	
	  t6.re=*X6*w6.re;
	  t6.im=*X6*w6.im;
	  t7.re=*X7*w7.re;
	  t7.im=*X7*w7.im;
	  
	  s0.re=t0.re+t4.re;
	  s0.im=t0.im+t4.im;
	  s1.re=t2.re+t6.re;
	  s1.im=t2.im+t6.im;
	  s2.re=t1.re+t5.re;
	  s2.im=t1.im+t5.im;
	  s3.re=t3.re+t7.re;
	  s3.im=t3.im+t7.im;
	  s4.re=t0.re-t4.re;
	  s4.im=t0.im-t4.im;
	  s5.re=t2.re-t6.re;
	  s5.im=t2.im-t6.im;
	  s6.re=t1.re-t5.re;
	  s6.im=t1.im-t5.im;
	  s7.re=t3.re-t7.re;
	  s7.im=t3.im-t7.im;
	  
	  t0.re=  ( s0.re+s1.re);
	  t0.im=  ( s0.im+s1.im);
	  t1.re=  ( s2.re+s3.re);
	  t1.im=  ( s2.im+s3.im);
	  t2.re=  ( s4.re+s5.im);
	  t2.im=  ( s4.im-s5.re);
	  t3.re=c*( s6.re+s6.im-s7.re+s7.im);
	  t3.im=c*( s6.im-s6.re-s7.im-s7.re);
	  t4.re=  ( s0.re-s1.re);
	  t4.im=  ( s0.im-s1.im);
	  t5.re=  ( s2.im-s3.im);
	  t5.im=  (-s2.re+s3.re);
	  t6.re=  ( s4.re-s5.im);
	  t6.im=  ( s4.im+s5.re);
	  t7.re=c*(-s6.re+s6.im+s7.re+s7.im);
	  t7.im=c*(-s6.im-s6.re+s7.im-s7.re);
	  y0->re = (t0.re+t1.re);
	  y0->im = (t0.im+t1.im);
	  y1->re = (t2.re+t3.re);
	  y1->im = (t2.im+t3.im);
	  y2->re = (t4.re+t5.re);
	  y2->im = (t4.im+t5.im);
	  y3->re = (t6.re+t7.re);
	  y3->im = (t6.im+t7.im);
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

