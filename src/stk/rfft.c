/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfft.c,v 1.1 2001-08-06 20:33:33 mt Exp $
 */

#include <stdlib.h>
#include <math.h>
#include "rfft.h"

int rfftFactor(int n,int m)
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

rfft *rfftInit(int n)
{
  rfft *X;
  rfftNode *Y,*Z;
  real angle;
  real pi=STK_PI;
  complex *R;
  int k,r,m,L=1;
  m=n;
  X=(rfft *)malloc(sizeof(rfft));
  X->n=n;
  X->work=(real *)calloc(n,sizeof(real));
  Y=(rfftNode *)malloc(sizeof(rfftNode));
  X->top=Y;
  Y->prev=NULL;
  while(m>1)
    {
      r=rfftFactor(m,n);
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
	  R=Y->aux;
	  for(k=0;k<r;k++)
	    {
	      angle=-((real)(k))*(2.0*pi/(real)r);
	      R[k].re=cos(angle);
	      R[k].im=sin(angle);
	    }
	  break;
	}
      L=L*r;
      if(m>1)
	{
	  Z=(rfftNode *)malloc(sizeof(rfftNode));
	  Y->next=(void *)Z;
	  Z->prev=(void *)Y;
	  Y=(rfftNode *)Y->next;
	}
      else
	{
	  Y->next=NULL;
	}
    }
  X->bot=(void *)Y;
  return(X);
}

void *rfftCopyNode(rfftNode *in,rfftNode *out)
{
  int r,k;
  rfftNode *M,*N;
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
      N=(rfftNode *)malloc(sizeof(rfftNode));
      if(N==NULL)
	return(NULL);
      out->next=N;
      N->prev=out;
      M=(rfftNode *)(in->next);
      p=rfftCopyNode(M,N);
      if(p==NULL)
	return(NULL);
    }
  else
    out->next=NULL;
  return(p);
}


rfft *rfftCopy(rfft *X)
{
  rfft *Y;
  rfftNode *S,*T;
  void *P;
  Y=(rfft *)malloc(sizeof(rfft));
  if(Y==NULL)
    return(NULL);
  Y->n=X->n;
  Y->work=(real *)calloc((Y->n),sizeof(real));
  if(Y->work==NULL)
    return(NULL);
  T=(rfftNode *)malloc(sizeof(rfftNode));  
  if(T==NULL)
    return(NULL);
  Y->top=T;
  T->prev=NULL;
  S=(rfftNode *)(X->top);
  P=rfftCopyNode(S,T);
  if(P==NULL)
    return(NULL);
  Y->bot=P;
  return(Y);
}

void rfftDeleteNode(rfftNode *X,int n)
{
  int r;
  r=X->r;
  if(n/r>1)
    rfftDeleteNode((rfftNode *)(X->next),n/(X->r));
  if(!((r==2)||(r==3)||(r==4)||(r==5)||(r==8)))
    free(X->aux);
  free(X);
}

void rfftDelete(rfft *X)
{
  rfftDeleteNode((rfftNode *)(X->top),X->n);
  free(X->work);
  free(X);
}


