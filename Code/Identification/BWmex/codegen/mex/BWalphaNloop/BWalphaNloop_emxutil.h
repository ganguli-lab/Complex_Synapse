/*
 * BWalphaNloop_emxutil.h
 *
 * Code generation for function 'BWalphaNloop_emxutil'
 *
 * C source code generated on: Fri Jun 27 13:42:24 2014
 *
 */

#ifndef __BWALPHANLOOP_EMXUTIL_H__
#define __BWALPHANLOOP_EMXUTIL_H__
/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "BWalphaNloop_types.h"

/* Function Declarations */
extern void b_emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray, int32_T numDimensions, const emlrtRTEInfo *srcLocation, boolean_T doPush);
extern void c_emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray, int32_T numDimensions, const emlrtRTEInfo *srcLocation, boolean_T doPush);
extern void emxEnsureCapacity(const emlrtStack *sp, emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize, const emlrtRTEInfo *srcLocation);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_int32_T(const emlrtStack *sp, emxArray_int32_T **pEmxArray, int32_T numDimensions, const emlrtRTEInfo *srcLocation, boolean_T doPush);
extern void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray, int32_T numDimensions, const emlrtRTEInfo *srcLocation, boolean_T doPush);
#endif
/* End of code generation (BWalphaNloop_emxutil.h) */
