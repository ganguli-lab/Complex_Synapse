/*
 * BWbetaNloop_initialize.c
 *
 * Code generation for function 'BWbetaNloop_initialize'
 *
 * C source code generated on: Thu Jun 26 11:57:47 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWbetaNloop.h"
#include "BWbetaNloop_initialize.h"
#include "BWbetaNloop_data.h"

/* Function Definitions */
void BWbetaNloop_initialize(emlrtStack *sp, emlrtContext *aContext)
{
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  sp->tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(sp, FALSE, 0U, 0);
  emlrtEnterRtStackR2012b(sp);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (BWbetaNloop_initialize.c) */
