/*
 * BWalphaNloop_api.c
 *
 * Code generation for function 'BWalphaNloop_api'
 *
 * C source code generated on: Wed Jun 25 13:48:35 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWalphaNloop.h"
#include "BWalphaNloop_api.h"
#include "BWalphaNloop_emxutil.h"

/* Variable Definitions */
static emlrtRTEInfo b_emlrtRTEI = { 1, 1, "BWalphaNloop_api", "" };

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void b_emlrt_marshallOut(emxArray_real_T *u, const mxArray *y);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *eta, const
  char_T *identifier, emxArray_real_T *y);
static void c_emlrt_marshallOut(emxArray_real_T *u, const mxArray *y);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *updater,
  const char_T *identifier, emxArray_real_T *y);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *alpha, const
  char_T *identifier, emxArray_real_T *y);
static void emlrt_marshallOut(emxArray_real_T *u, const mxArray *y);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  g_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void b_emlrt_marshallOut(emxArray_real_T *u, const mxArray *y)
{
  mxSetData((mxArray *)y, (void *)u->data);
  mxSetDimensions((mxArray *)y, u->size, 1);
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *eta, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  d_emlrt_marshallIn(sp, emlrtAlias(eta), &thisId, y);
  emlrtDestroyArray(&eta);
}

static void c_emlrt_marshallOut(emxArray_real_T *u, const mxArray *y)
{
  mxSetData((mxArray *)y, (void *)u->data);
  mxSetDimensions((mxArray *)y, u->size, 3);
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  h_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *updater,
  const char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  f_emlrt_marshallIn(sp, emlrtAlias(updater), &thisId, y);
  emlrtDestroyArray(&updater);
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *alpha, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  b_emlrt_marshallIn(sp, emlrtAlias(alpha), &thisId, y);
  emlrtDestroyArray(&alpha);
}

static void emlrt_marshallOut(emxArray_real_T *u, const mxArray *y)
{
  mxSetData((mxArray *)y, (void *)u->data);
  mxSetDimensions((mxArray *)y, u->size, 2);
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  i_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv2[2];
  boolean_T bv0[2];
  int32_T i;
  int32_T iv3[2];
  for (i = 0; i < 2; i++) {
    iv2[i] = -1;
    bv0[i] = TRUE;
  }

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", FALSE, 2U, iv2, bv0, iv3);
  ret->size[0] = iv3[0];
  ret->size[1] = iv3[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv4[1];
  boolean_T bv1[1];
  int32_T iv5[1];
  iv4[0] = -1;
  bv1[0] = TRUE;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", FALSE, 1U, iv4, bv1, iv5);
  ret->size[0] = iv5[0];
  ret->allocatedSize = ret->size[0];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv6[3];
  boolean_T bv2[3];
  int32_T i;
  int32_T iv7[3];
  for (i = 0; i < 3; i++) {
    iv6[i] = -1;
    bv2[i] = TRUE;
  }

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", FALSE, 3U, iv6, bv2, iv7);
  ret->size[0] = iv7[0];
  ret->size[1] = iv7[1];
  ret->size[2] = iv7[2];
  ret->allocatedSize = ret->size[0] * ret->size[1] * ret->size[2];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

void BWalphaNloop_api(emlrtStack *sp, const mxArray *prhs[3], const mxArray
                      *plhs[3])
{
  emxArray_real_T *alpha;
  emxArray_real_T *eta;
  emxArray_real_T *updater;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &alpha, 2, &b_emlrtRTEI, TRUE);
  b_emxInit_real_T(sp, &eta, 1, &b_emlrtRTEI, TRUE);
  c_emxInit_real_T(sp, &updater, 3, &b_emlrtRTEI, TRUE);
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, TRUE, -1);
  prhs[1] = emlrtProtectR2012b(prhs[1], 1, TRUE, -1);
  prhs[2] = emlrtProtectR2012b(prhs[2], 2, TRUE, -1);

  /* Marshall function inputs */
  emlrt_marshallIn(sp, emlrtAlias(prhs[0]), "alpha", alpha);
  c_emlrt_marshallIn(sp, emlrtAlias(prhs[1]), "eta", eta);
  e_emlrt_marshallIn(sp, emlrtAlias(prhs[2]), "updater", updater);

  /* Invoke the target function */
  BWalphaNloop(sp, alpha, eta, updater);

  /* Marshall function outputs */
  emlrt_marshallOut(alpha, prhs[0]);
  plhs[0] = prhs[0];
  b_emlrt_marshallOut(eta, prhs[1]);
  plhs[1] = prhs[1];
  c_emlrt_marshallOut(updater, prhs[2]);
  plhs[2] = prhs[2];
  updater->canFreeData = FALSE;
  emxFree_real_T(&updater);
  eta->canFreeData = FALSE;
  emxFree_real_T(&eta);
  alpha->canFreeData = FALSE;
  emxFree_real_T(&alpha);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (BWalphaNloop_api.c) */
