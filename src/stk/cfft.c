/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: cfft.c,v 1.1 2001-08-06 20:33:30 mt Exp $
 */

#include <stdlib.h>
#include <math.h>
#include "cfft.h"

int cfftFactor(int n,int m)
{
  int f;
  if((n%8==0)&&(m>=RADIX_8_MIN))
    return(8); 
  else if((n%4==0)&&(m>=RADIX_4_MIN))
    return(4);
  else if(n%2==0)
    return(2);
  else if(n%3==0)
    return(3);
  else if(n%5==0)
    return(5);
  else
    for(f=7;f*f<=n;f+=2)
      if(n%f==0)
	return(f);
  return(n);
}

cfft *cfftInit(int n)
{
  cfft *X;
  cfftNode *Y;
  real angle;
  real pi=STK_PI;
  complex *R;
  int k,r,m,L=1;
  m=n;
  X=(cfft *)malloc(sizeof(cfft));
  if(X==NULL)
    return(NULL);
  X->n=n;
  X->work=(complex *)calloc(n,sizeof(complex));
  if(X->work==NULL)
    return(NULL);
  Y=(cfftNode *)malloc(sizeof(cfftNode));
  if(Y==NULL)
    return(NULL);
  X->top=Y;
  while(m>1)
    {
      r=cfftFactor(m,n);
      m=m/r;
      angle=-2.0*pi/((real)(L*r));
      Y->omega.re=-2.0*sin(0.5*angle)*sin(0.5*angle);
      Y->omega.im=sin(angle);
      Y->r=r;
      switch(r)
	{
	case 2:
	  break;
	case 3:
	  break;
	case 4:
	  break;
	case 5:
	  break;
	case 8:
	  break;
	default:
	  Y->aux=(complex *)calloc(2*r,sizeof(complex));
	  if(Y->aux==NULL)
	    return(NULL);
	  R=Y->aux;
	  for(k=0;k<r;k++)
	    {
	      angle=-((real)k)*(2.0*pi/(real)r);
	      R[k].re=cos(angle);
	      R[k].im=sin(angle);
	    }
	  break;
	}
      L=L*r;
      if(m>1)
	{
	  Y->next=malloc(sizeof(cfftNode));
	  if(Y->next==NULL)
	    return(NULL);
	  Y=(cfftNode *)(Y->next);
	}
      else
	{
	  Y->next=NULL;
	}
    }
  return(X);
}

void *cfftCopyNode(cfftNode *in,cfftNode *out)
{
  int r,k;
  cfftNode *M,*N;
  void *p;

  r=in->r;
  out->r=r;
  out->omega=in->omega;
  switch(r)
    {
    case 2:
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      break;
    case 8:
      break;
    default:
      out->aux=(complex *)calloc(2*r,sizeof(complex));
      if(out->aux==NULL)
	return(NULL);
      for(k=0;k<r;k++)
	{
	  out->aux[k].re=in->aux[k].re;
	  out->aux[k].im=in->aux[k].im;
	}
      break;
    }
  p=out;
  if(in->next!=NULL)
    {
      N=(cfftNode *)malloc(sizeof(cfftNode));
      if(N==NULL)
	return(NULL);
      out->next=N;
      M=(cfftNode *)(in->next);
      p=cfftCopyNode(M,N);
      if(p==NULL)
	return(NULL);
    }
  else
    out->next=NULL;
  return(p);
}

cfft *cfftCopy(cfft *X)
{
  cfft *Y;
  cfftNode *S,*T;
  void *P;
  Y=(cfft *)malloc(sizeof(cfft));
  if(Y==NULL)
    return(NULL);
  Y->n=X->n;
  Y->work=(complex *)calloc((Y->n),sizeof(complex));
  if(Y->work==NULL)
    return(NULL);
  T=(cfftNode *)malloc(sizeof(cfftNode));  
  if(T==NULL)
    return(NULL);
  Y->top=T;
  S=(cfftNode *)(X->top);
  P=cfftCopyNode(S,T);
  if(P==NULL)
    return(NULL);
  return(Y);
}

void cfftDeleteNode(cfftNode *X,int n)
{
  int r;
  r=X->r;
  if(n/r>1)
    cfftDeleteNode((cfftNode *)(X->next),n/(X->r));
  if(!((r==2)||(r==3)||(r==4)||(r==5)||(r==8)))
    free(X->aux);
  free(X);
}

void cfftDelete(cfft *X)
{
  cfftDeleteNode((cfftNode *)(X->top),X->n);
  free(X->work);
  free(X);
}


