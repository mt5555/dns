/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: fortran.c,v 1.1 2001-08-06 20:33:32 mt Exp $
 */


#include "cfft.h"
#include "rfft.h"
#include <stdlib.h>
#include <stdio.h>

/* fortran wrapper macro, appends correct suffix (e.q. "_") */

#if defined (__hpux)
#define _Fortran(A) A
#elif defined (__sun)
#define _Fortran(A) A##_
#elif defined (linux)
#define _Fortran(A) A##_
#elif defined (__alpha)
#define _Fortran(A) A
#elif defined (__sgi)
#define _Fortran(A) A##_
#elif defined (_AIX)
#define _Fortran(A) A
#else
#define _Fortran(A) A
#endif

/* fortran api error handler, modify as needed */

void fortranError(char *description)
{
  fprintf(stderr,"%s\n",description);
  exit(-1);
}

void _Fortran(rfft_init) (int *n,rfft **X)
{
  *X=rfftInit(*n);
  if(*X==NULL)
    fortranError("rfft_init:: cannot allocate object");
}

void _Fortran(cfft_init) (int *n,cfft **X)
{
  *X=cfftInit(*n);
  if(*X==NULL)
    fortranError("cfft_init:: cannot allocate object");
}

void _Fortran(rfft_delete) (rfft **X)
{
  rfftDelete(*X);
}

void _Fortran(cfft_delete) (cfft **X)
{
  cfftDelete(*X);
}

void _Fortran(rfft_copy) (rfft **Y,rfft **X)
{
  *X=rfftCopy(*Y);
  if(*X==NULL)
    fortranError("rfft_copy:: cannot allocate object");
}

void _Fortran(cfft_copy) (cfft **Y,cfft **X)
{
  *X=cfftCopy(*Y);
  if(*X==NULL)
    fortranError("cfft_copy:: cannot allocate object");
}

void _Fortran(cfft_analysis) (complex *c,cfft **X) 
{ 
  cfftAnalysis(c,*X); 
}

void _Fortran(cfft_synthesis) (complex *c,cfft **X) 
{ 
  cfftSynthesis(c,*X); 
}


void _Fortran(rfft_analysis) (real *c,rfft **X) 
{ 
  rfftAnalysis(c,*X); 
}

void _Fortran(rfft_synthesis) (real *c,rfft **X) 
{ 
  rfftSynthesis(c,*X); 
}

void _Fortran(cfft_analysis_m) (complex *c,cfft **X,int *m,int *n) 
{ 
  int i;
  if((*m<=0)||(*n<(*X)->n))
    fortranError("cfft_analysis_m:: argument out of range");
  cfftAnalysisM(c,*X,*n,0,*m);
}

void _Fortran(cfft_synthesis_m) (complex *c,cfft **X,int *m,int *n) 
{ 
  int i;
  if((*m<=0)||(*n<(*X)->n))
    fortranError("cfft_synthesis_m:: argument out of range");
  cfftSynthesisM(c,*X,*n,0,*m);
}


void _Fortran(rfft_analysis_m) (real *c,rfft **X,int *m,int *n) 
{ 
  int i;
  if((*m<=0)||(*n<(*X)->n))
    fortranError("rfft_analysis_m:: argument out of range");
  rfftAnalysisM(c,*X,*n,0,*m);
}

void _Fortran(rfft_synthesis_m) (real *c,rfft **X,int *m,int *n) 
{ 
  int i;
  if((*m<=0)||(*n<(*X)->n))
    fortranError("rfft_synthesis_m:: argument out of range");
  rfftSynthesisM(c,*X,*n,0,*m);
}

