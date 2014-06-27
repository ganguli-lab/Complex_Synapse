/*
 * BWbetaNloop_mex.c
 *
 * Code generation for function 'BWbetaNloop'
 *
 * C source code generated on: Thu Jun 26 11:57:47 2014
 *
 */

/* Include files */
#include "mex.h"
#include "BWbetaNloop_api.h"
#include "BWbetaNloop_initialize.h"
#include "BWbetaNloop_terminate.h"

/* Function Declarations */
static void BWbetaNloop_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "BWbetaNloop", NULL, false, {2045744189U,2170104910U,2743257031U,4284093946U}, NULL };
void *emlrtRootTLSGlobal = NULL;

/* Function Definitions */
static void BWbetaNloop_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *outputs[1];
  mxArray *inputs[2];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  int nInputs = nrhs;
  emlrtStack stack={0,0,0}; /* Root of the run-time stack. */
  /* Module initialization. */
  BWbetaNloop_initialize(&stack, &emlrtContextGlobal);
  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:WrongNumberOfInputs", 5, mxINT32_CLASS, 2, mxCHAR_CLASS, 11, "BWbetaNloop");
  } else if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:TooManyOutputArguments", 3, mxCHAR_CLASS, 11, "BWbetaNloop");
  }
  /* Temporary copy for mex inputs. */
  for (n = 0; n < nInputs; ++n) {
    inputs[n] = (mxArray *)prhs[n];
  }
  /* Call the function. */
  BWbetaNloop_api(&stack, (const mxArray**)inputs, (const mxArray**)outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  BWbetaNloop_terminate(&stack);
}

void BWbetaNloop_atexit_wrapper(void)
{
  emlrtStack stack={0,0,0}; /* Root of the run-time stack. */
   BWbetaNloop_atexit(&stack);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(BWbetaNloop_atexit_wrapper);
  /* Dispatch the entry-point. */
  BWbetaNloop_mexFunction(nlhs, plhs, nrhs, prhs);
}
/* End of code generation (BWbetaNloop_mex.c) */
