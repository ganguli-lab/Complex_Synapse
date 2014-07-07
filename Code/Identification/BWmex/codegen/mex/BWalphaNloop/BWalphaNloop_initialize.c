/*
 * BWalphaNloop_initialize.c
 *
 * Code generation for function 'BWalphaNloop_initialize'
 *
 * C source code generated on: Mon Jul 07 15:22:55 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWalphaNloop.h"
#include "BWalphaNloop_initialize.h"

/* Variable Definitions */
static const volatile char_T *emlrtBreakCheckR2012bFlagVar;

/* Function Definitions */
void BWalphaNloop_initialize(emlrtContext *aContext)
{
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, FALSE, 0U, 0);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (BWalphaNloop_initialize.c) */
