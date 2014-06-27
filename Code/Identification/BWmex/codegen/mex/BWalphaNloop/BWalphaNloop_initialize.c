/*
 * BWalphaNloop_initialize.c
 *
 * Code generation for function 'BWalphaNloop_initialize'
 *
 * C source code generated on: Thu Jun 26 11:57:26 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWalphaNloop.h"
#include "BWalphaNloop_initialize.h"
#include "BWalphaNloop_data.h"

/* Function Definitions */
void BWalphaNloop_initialize(emlrtStack *sp, emlrtContext *aContext)
{
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  sp->tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(sp, FALSE, 0U, 0);
  emlrtEnterRtStackR2012b(sp);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (BWalphaNloop_initialize.c) */
