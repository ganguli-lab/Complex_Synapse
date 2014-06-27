/*
 * BWalphaNloop_terminate.c
 *
 * Code generation for function 'BWalphaNloop_terminate'
 *
 * C source code generated on: Fri Jun 27 13:42:24 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWalphaNloop.h"
#include "BWalphaNloop_terminate.h"

/* Function Definitions */
void BWalphaNloop_atexit(emlrtStack *sp)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  sp->tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(sp);
  emlrtLeaveRtStackR2012b(sp);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void BWalphaNloop_terminate(emlrtStack *sp)
{
  emlrtLeaveRtStackR2012b(sp);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (BWalphaNloop_terminate.c) */
