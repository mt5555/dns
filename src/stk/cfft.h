/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: cfft.h,v 1.1 2001-08-06 20:33:30 mt Exp $
 */

#ifndef _CFFT_H
#define _CFFT_H

#include "real.h"
#include "constants.h"

#define RADIX_8_MIN 512 /* min n to use radix-8 */
#define RADIX_4_MIN 0   /* min n to use radix-4 */

typedef struct
{
  complex *work;
  int n;
  void *top;
} cfft;

typedef struct
{
  complex omega;
  complex *aux;
  int r;
  void *next;
} cfftNode;

#ifdef __cplusplus
#define _C_EXTERN "C"
#else
#define _C_EXTERN
#endif

/*-------------------------- public prototypes ---------------------------- */

extern _C_EXTERN cfft *cfftInit(int);
extern _C_EXTERN void cfftDelete(cfft *);
extern _C_EXTERN cfft *cfftCopy(cfft *);
extern _C_EXTERN void cfftAnalysis(complex *,cfft *);
extern _C_EXTERN void cfftSynthesis(complex *,cfft *);
extern _C_EXTERN void cfftAnalysisM(complex *,cfft *,int,int,int);
extern _C_EXTERN void cfftSynthesisM(complex *,cfft *,int,int,int);

/*------------------------- private prototypes ----------------------------- */

extern _C_EXTERN int cfftFactor(int,int);
extern _C_EXTERN void *cfftCopyNode(cfftNode *,cfftNode *);
extern _C_EXTERN void cfftDeleteNode(cfftNode *,int);
extern _C_EXTERN void cfftAnalysisR2(complex *,complex *,complex,int,int);
extern _C_EXTERN void cfftAnalysisR3(complex *,complex *,complex,int,int);
extern _C_EXTERN void cfftAnalysisR4(complex *,complex *,complex,int,int);
extern _C_EXTERN void cfftAnalysisR5(complex *,complex *,complex,int,int);
extern _C_EXTERN void cfftAnalysisR8(complex *,complex *,complex,int,int);
extern _C_EXTERN void cfftAnalysisRg(complex *,complex *,complex *,
			   complex *,complex,int,int,int);

extern _C_EXTERN void cfftSynthesisR2(complex *,complex *,complex,int,int);
extern _C_EXTERN void cfftSynthesisR3(complex *,complex *,complex,int,int);
extern _C_EXTERN void cfftSynthesisR4(complex *,complex *,complex,int,int);
extern _C_EXTERN void cfftSynthesisR5(complex *,complex *,complex,int,int);
extern _C_EXTERN void cfftSynthesisR8(complex *,complex *,complex,int,int);
extern _C_EXTERN void cfftSynthesisRg(complex *,complex *,complex *,
			    complex *,complex,int,int,int);

#endif

