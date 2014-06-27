/*
 * eml_error.c
 *
 * Code generation for function 'eml_error'
 *
 * C source code generated on: Thu Jun 26 11:57:47 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWbetaNloop.h"
#include "eml_error.h"

/* Variable Definitions */
static emlrtRTEInfo d_emlrtRTEI = { 20, 5, "eml_error",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_error.m" };

/* Function Definitions */
void eml_error(const emlrtStack *sp)
{
  emlrtErrorWithMessageIdR2012b(sp, &d_emlrtRTEI,
    "Coder:toolbox:reshape_emptyReshapeLimit", 0);
}

/* End of code generation (eml_error.c) */
