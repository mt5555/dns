/*
 * Copyright (C) 2001 University Corporation for Atmospheric Research 
 *       
 * $Id: rfft.h,v 1.1 2001-08-06 20:33:33 mt Exp $
 */

#ifndef _RFFT_H
#define _RFFT_H

#include "real.h"
#include "constants.h"

#define RADIX_8_MIN 512 /* min n to use radix-8 */
#define RADIX_4_MIN 0   /* min n to use radix-4 */

typedef struct
{
  real *work;
  int n;
  void *top;
  void *bot;
} rfft;

typedef struct
{
  complex omega;
  complex *aux;
  int r;
  void *next;
  void *prev;
} rfftNode;

#ifdef __cplusplus
#define _C_EXTERN "C"
#else
#define _C_EXTERN
#endif

/*-------------------------- public prototypes ---------------------------- */

extern _C_EXTERN rfft *rfftInit(int);
extern _C_EXTERN void rfftDelete(rfft *);
extern _C_EXTERN rfft *rfftCopy(rfft *);
extern _C_EXTERN void rfftAnalysis(real *,rfft *);
extern _C_EXTERN void rfftSynthesis(real *,rfft *);
extern _C_EXTERN void rfftAnalysisM(real *,rfft *,int,int,int);
extern _C_EXTERN void rfftSynthesisM(real *,rfft *,int,int,int);

/*------------------------- private prototypes ----------------------------- */

extern _C_EXTERN int rfftFactor(int,int);
extern _C_EXTERN void *rfftCopyNode(rfftNode *,rfftNode *);
extern _C_EXTERN void rfftDeleteNode(rfftNode *,int);
extern _C_EXTERN void rfftAnalysisR2(real *,real *,complex,int,int);
extern _C_EXTERN void rfftAnalysisR3(real *,real *,complex,int,int);
extern _C_EXTERN void rfftAnalysisR4(real *,real *,complex,int,int);
extern _C_EXTERN void rfftAnalysisR5(real *,real *,complex,int,int);
extern _C_EXTERN void rfftAnalysisR8(real *,real *,complex,int,int);
extern _C_EXTERN void rfftAnalysisRg(real *,real *,complex *,complex *,
			   complex,int,int,int);

extern _C_EXTERN void rfftSynthesisR2(real *,real *,complex,int,int);
extern _C_EXTERN void rfftSynthesisR3(real *,real *,complex,int,int);
extern _C_EXTERN void rfftSynthesisR4(real *,real *,complex,int,int);
extern _C_EXTERN void rfftSynthesisR5(real *,real *,complex,int,int);
extern _C_EXTERN void rfftSynthesisR8(real *,real *,complex,int,int);
extern _C_EXTERN void rfftSynthesisRg(real *,real *,complex *,complex *,
			    complex,int,int,int);

#endif

