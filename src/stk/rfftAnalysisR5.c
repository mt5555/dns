/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfftAnalysisR5.c,v 1.1 2001-08-06 20:33:33 mt Exp $
 */

#include "rfft.h"

void rfftAnalysisR5(real *x,real *y,complex omega,int L,int N)
{
  register complex t0,t1,t2,t3,t4;
  register complex w1,w2,w3,w4;
  register complex z1,z2,z3,z4;
  complex *x0,*x1,*x2,*x3,*x4,*y0,*y1,*y2,*y3,*y4;
  real *X0,*X1,*X2,*X3,*X4;
  real *Y0;
  complex a1={STK_R5C2,STK_R5S2};
  complex a2={STK_R5C4,STK_R5S4};
  int j,k;
  int NL,L5;
  NL=N*L;
  L5=L*5;
  X0=x;
  X1=X0+NL;
  X2=X1+NL;
  X3=X2+NL;
  X4=X3+NL;
  for(j=0;j<N;j++)
    {
      Y0=y+j*L5;
      y1=(complex *)(Y0-1)+L;
      y2=y1+L;
      y3=y1-1;
      y4=y2-1;
      t0.re=*X0;
      t1.re=*X1; 
      t2.re=*X2; 
      t3.re=*X3; 
      t4.re=*X4; 
      z1.re=t1.re+t4.re;
      z2.re=t1.re-t4.re;
      z3.re=t2.re+t3.re;
      z4.re=t2.re-t3.re;
      *Y0=t0.re+z1.re+z3.re;      
      y1->re=t0.re+a1.re*z1.re+a2.re*z3.re;
      y1->im=      a1.im*z2.re+a2.im*z4.re;	  
      y2->re=t0.re+a2.re*z1.re+a1.re*z3.re;
      y2->im=      a2.im*z2.re-a1.im*z4.re;
      w1.re=omega.re+1.0;                   
      w1.im=omega.im;
      X0++;
      X1++;
      X2++;
      X3++;
      X4++;
      
      Y0++;
      y1++;
      y2++;
      
      y0=(complex *)Y0;
      x0=(complex *)X0; 
      x1=(complex *)X1; 
      x2=(complex *)X2; 
      x3=(complex *)X3; 
      x4=(complex *)X4; 

      for(k=1;k<(L+1)/2;k++)
	{
	  w2.re = w1.re*w1.re-w1.im*w1.im;
	  w2.im = 2.0*w1.re*w1.im;
	  w3.re = w2.re*w1.re-w2.im*w1.im;
	  w3.im = w2.re*w1.im+w2.im*w1.re; 
	  w4.re = w3.re*w1.re-w3.im*w1.im;
	  w4.im = w3.re*w1.im+w3.im*w1.re;
	  t0.re=x0->re;
	  t0.im=x0->im;
	  z1.re=x1->re;
	  z1.im=x1->im;
	  z2.re=x2->re;
	  z2.im=x2->im;
	  z3.re=x3->re;
	  z3.im=x3->im;
	  z4.re=x4->re;
	  z4.im=x4->im;
	  t1.re=z1.re*w1.re-z1.im*w1.im;
	  t1.im=z1.im*w1.re+z1.re*w1.im;
	  t2.re=z2.re*w2.re-z2.im*w2.im;
	  t2.im=z2.im*w2.re+z2.re*w2.im;
	  t3.re=z3.re*w3.re-z3.im*w3.im;
	  t3.im=z3.im*w3.re+z3.re*w3.im;
	  t4.re=z4.re*w4.re-z4.im*w4.im;
	  t4.im=z4.im*w4.re+z4.re*w4.im;
	  
	  z1.re=t1.re+t4.re;
	  z1.im=t1.im+t4.im;
        z2.re=t1.re-t4.re;
        z2.im=t1.im-t4.im;
        z3.re=t2.re+t3.re;
        z3.im=t2.im+t3.im;
        z4.re=t2.re-t3.re;
        z4.im=t2.im-t3.im;
        
	  y0->re=t0.re+z1.re+z3.re;
	  y0->im=t0.im+z1.im+z3.im;
	  
	  y1->re= (t0.re+a1.re*z1.re-a1.im*z2.im+a2.re*z3.re-a2.im*z4.im);
	    
	  y1->im= (t0.im+a1.im*z2.re+a1.re*z1.im+a2.im*z4.re+a2.re*z3.im);
	
	  y2->re= (t0.re+a2.re*z1.re-a2.im*z2.im+a1.re*z3.re+a1.im*z4.im);
		 
	  y2->im= (t0.im+a2.im*z2.re+a2.re*z1.im-a1.im*z4.re+a1.re*z3.im);
		  
	  y3->re= (t0.re+a1.re*z1.re+a1.im*z2.im+a2.re*z3.re+a2.im*z4.im);
		
	  y3->im=-(t0.im-a1.im*z2.re+a1.re*z1.im-a2.im*z4.re+a2.re*z3.im);
		         
	  y4->re= (t0.re+a2.re*z1.re+a2.im*z2.im+a1.re*z3.re-a1.im*z4.im);
		         
	  y4->im=-(t0.im-a2.im*z2.re+a2.re*z1.im+a1.im*z4.re+a1.re*z3.im);
		        
		          
	  z1.re=omega.re*w1.re-omega.im*w1.im+w1.re;
	  z1.im=omega.re*w1.im+omega.im*w1.re+w1.im; 
	  w1.re=z1.re;
	  w1.im=z1.im;
	  x0++;
	  x1++;
	  x2++;
	  x3++;
	  x4++;
	  y0++;
	  y1++;
	  y2++;
	  y3--;
	  y4--;
	}
	X0=(real *)x0;
	X1=(real *)x1;
	X2=(real *)x2;
	X3=(real *)x3;
	X4=(real *)x4;
      if(L%2==0)
	{
	  w2.re = w1.re*w1.re-w1.im*w1.im;
	  w2.im = 2.0*w1.re*w1.im;
	  w3.re = w2.re*w1.re-w2.im*w1.im;
	  w3.im = w2.re*w1.im+w2.im*w1.re; 
	  w4.re = w3.re*w1.re-w3.im*w1.im;
	  w4.im = w3.re*w1.im+w3.im*w1.re;
	  t0.re=*X0;
	  t1.re=*X1*w1.re;
	  t1.im=*X1*w1.im;
	  t2.re=*X2*w2.re;
	  t2.im=*X2*w2.im;
	  t3.re=*X3*w3.re;
	  t3.im=*X3*w3.im;
	  t4.re=*X4*w4.re;
	  t4.im=*X4*w4.im;
	  z1.re=t1.re+t4.re;
	  z1.im=t1.im+t4.im;
        z2.re=t1.re-t4.re;
        z2.im=t1.im-t4.im;
        z3.re=t2.re+t3.re;
        z3.im=t2.im+t3.im;
        z4.re=t2.re-t3.re;
        z4.im=t2.im-t3.im;
        
	  y0->re=t0.re+z1.re+z3.re;
	  
	  y0->im=z1.im+z3.im;
	  
	  y1->re=(t0.re+a1.re*z1.re-a1.im*z2.im+a2.re*z3.re-a2.im*z4.im);
	  
	  y1->im=(      a1.im*z2.re+a1.re*z1.im+a2.im*z4.re+a2.re*z3.im);
	  
	  y2->re=(t0.re+a2.re*z1.re-a2.im*z2.im+a1.re*z3.re+a1.im*z4.im);
	  
	  X0++;
	  X1++;
	  X2++;
	  X3++;
	  X4++;
	    
	}
    }
}

