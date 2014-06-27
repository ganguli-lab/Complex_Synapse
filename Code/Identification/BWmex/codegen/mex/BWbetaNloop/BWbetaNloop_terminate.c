/*
 * BWbetaNloop_terminate.c
 *
 * Code generation for function 'BWbetaNloop_terminate'
 *
 * C source code generated on: Thu Jun 26 11:57:47 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWbetaNloop.h"
#include "BWbetaNloop_terminate.h"

/* Function Definitions */
void BWbetaNloop_atexit(emlrtStack *sp)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  sp->tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(sp);
  emlrtLeaveRtStackR2012b(sp);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void BWbetaNloop_terminate(emlrtStack *sp)
{
  emlrtLeaveRtStackR2012b(sp);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (BWbetaNloop_terminate.c) */
