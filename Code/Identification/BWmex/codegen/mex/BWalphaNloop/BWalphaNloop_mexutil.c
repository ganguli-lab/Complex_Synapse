/*
 * BWalphaNloop_mexutil.c
 *
 * Code generation for function 'BWalphaNloop_mexutil'
 *
 * C source code generated on: Thu Jun 26 11:57:26 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWalphaNloop.h"
#include "BWalphaNloop_mexutil.h"

/* Function Definitions */
void error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(sp, 0, NULL, 1, &pArray, "error", TRUE, location);
}

/* End of code generation (BWalphaNloop_mexutil.c) */
