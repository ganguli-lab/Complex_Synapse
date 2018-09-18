/*
 * BWalphaNloop_mex.c
 *
 * Code generation for function 'BWalphaNloop'
 *
 * C source code generated on: Mon Jul 07 15:22:55 2014
 *
 */

/* Include files */
#include "mex.h"
#include "BWalphaNloop_api.h"
#include "BWalphaNloop_initialize.h"
#include "BWalphaNloop_terminate.h"

/* Function Declarations */
static void BWalphaNloop_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "BWalphaNloop", NULL, false, {2045744189U,2170104910U,2743257031U,4284093946U}, NULL };
void *emlrtRootTLSGlobal = NULL;

/* Function Definitions */
static void BWalphaNloop_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *outputs[3];
  mxArray *inputs[3];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  int nInputs = nrhs;
  /* Module initialization. */
  BWalphaNloop_initialize(&emlrtContextGlobal);
  /* Check for proper number of arguments. */
  if (nrhs != 3) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:WrongNumberOfInputs", 5, mxINT32_CLASS, 3, mxCHAR_CLASS, 12, "BWalphaNloop");
  } else if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:TooManyOutputArguments", 3, mxCHAR_CLASS, 12, "BWalphaNloop");
  }
  /* Temporary copy for mex inputs. */
  for (n = 0; n < nInputs; ++n) {
    inputs[n] = (mxArray *)prhs[n];
  }
  /* Call the function. */
  BWalphaNloop_api((const mxArray**)inputs, (const mxArray**)outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  BWalphaNloop_terminate();
}

void BWalphaNloop_atexit_wrapper(void)
{
   BWalphaNloop_atexit();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(BWalphaNloop_atexit_wrapper);
  /* Dispatch the entry-point. */
  BWalphaNloop_mexFunction(nlhs, plhs, nrhs, prhs);
}
/* End of code generation (BWalphaNloop_mex.c) */
