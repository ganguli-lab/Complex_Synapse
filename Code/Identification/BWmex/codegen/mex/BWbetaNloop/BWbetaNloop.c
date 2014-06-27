/*
 * BWbetaNloop.c
 *
 * Code generation for function 'BWbetaNloop'
 *
 * C source code generated on: Fri Jun 27 13:42:43 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWbetaNloop.h"
#include "BWbetaNloop_emxutil.h"
#include "BWbetaNloop_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 7, "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m"
};

static emlrtRSInfo b_emlrtRSI = { 64, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo c_emlrtRSI = { 54, "eml_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m" };

static emlrtRTEInfo emlrtRTEI = { 1, 21, "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m"
};

static emlrtRTEInfo c_emlrtRTEI = { 6, 1, "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m"
};

static emlrtBCInfo emlrtBCI = { -1, -1, 7, 31, "updater", "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m",
  0 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 7, 45, "beta", "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m",
  0 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 7, 12, "beta", "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m",
  0 };

static emlrtECInfo emlrtECI = { -1, 7, 5, "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m"
};

/* Function Definitions */
void BWbetaNloop(const emlrtStack *sp, emxArray_real_T *beta, const
                 emxArray_real_T *updater)
{
  int32_T i0;
  int32_T i;
  emxArray_int32_T *r0;
  emxArray_real_T *C;
  emxArray_real_T *a;
  emxArray_real_T *b;
  emxArray_int32_T *r1;
  real_T b_i;
  int32_T loop_ub;
  int32_T i1;
  int32_T i2;
  int32_T b_loop_ub;
  int32_T c_i;
  boolean_T guard1 = FALSE;
  real_T alpha1;
  real_T beta1;
  char_T TRANSB;
  char_T TRANSA;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  double * alpha1_t;
  double * Aia0_t;
  double * Bib0_t;
  double * beta1_t;
  double * Cic0_t;
  int32_T iv4[1];
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /* beta=BWBETANLOOP(eta,updater) normalised backward variables for Baum-Welch algorithm */
  /*    beta     = normalised backward variables */
  /*    updater  = updater matrices (M*outProj*eta) */
  i0 = beta->size[1];
  emlrtForLoopVectorCheckR2012b(beta->size[1], -1.0, 2.0, mxDOUBLE_CLASS,
    (int32_T)-(2.0 + (-1.0 - (real_T)beta->size[1])), &c_emlrtRTEI, sp);
  i = 0;
  emxInit_int32_T(sp, &r0, 1, &emlrtRTEI, TRUE);
  c_emxInit_real_T(sp, &C, 1, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &a, 2, &emlrtRTEI, TRUE);
  c_emxInit_real_T(sp, &b, 1, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r1, 1, &emlrtRTEI, TRUE);
  while (i <= (int32_T)-(2.0 + (-1.0 - (real_T)i0)) - 1) {
    b_i = (real_T)i0 + -(real_T)i;
    loop_ub = beta->size[0];
    i1 = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r0, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      r0->data[i1] = i1;
    }

    i1 = beta->size[1];
    i2 = (int32_T)(b_i - 1.0);
    emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &c_emlrtBCI, sp);
    st.site = &emlrtRSI;
    loop_ub = updater->size[0];
    b_loop_ub = updater->size[1];
    i1 = updater->size[2];
    i2 = (int32_T)(b_i - 1.0);
    c_i = emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &emlrtBCI, &st);
    i1 = a->size[0] * a->size[1];
    a->size[0] = loop_ub;
    a->size[1] = b_loop_ub;
    emxEnsureCapacity(&st, (emxArray__common *)a, i1, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        a->data[i2 + a->size[0] * i1] = updater->data[(i2 + updater->size[0] *
          i1) + updater->size[0] * updater->size[1] * (c_i - 1)];
      }
    }

    loop_ub = beta->size[0];
    i1 = beta->size[1];
    i2 = (int32_T)b_i;
    c_i = emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &b_emlrtBCI, &st);
    i1 = b->size[0];
    b->size[0] = loop_ub;
    emxEnsureCapacity(&st, (emxArray__common *)b, i1, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      b->data[i1] = beta->data[i1 + beta->size[0] * (c_i - 1)];
    }

    i1 = updater->size[1];
    guard1 = FALSE;
    if (i1 == 1) {
      guard1 = TRUE;
    } else {
      i1 = beta->size[0];
      if (i1 == 1) {
        guard1 = TRUE;
      } else {
        i1 = updater->size[0];
        b_st.site = &b_emlrtRSI;
        c_st.site = &c_emlrtRSI;
        i2 = C->size[0];
        C->size[0] = i1;
        emxEnsureCapacity(&c_st, (emxArray__common *)C, i2, (int32_T)sizeof
                          (real_T), &emlrtRTEI);
        loop_ub = i1;
        for (i1 = 0; i1 < loop_ub; i1++) {
          C->data[i1] = 0.0;
        }

        i1 = updater->size[0];
        if (i1 < 1) {
        } else {
          i1 = updater->size[1];
          if (i1 < 1) {
          } else {
            alpha1 = 1.0;
            beta1 = 0.0;
            TRANSB = 'N';
            TRANSA = 'N';
            i1 = updater->size[0];
            m_t = (ptrdiff_t)(i1);
            n_t = (ptrdiff_t)(1);
            i1 = updater->size[1];
            k_t = (ptrdiff_t)(i1);
            i1 = updater->size[0];
            lda_t = (ptrdiff_t)(i1);
            i1 = updater->size[1];
            ldb_t = (ptrdiff_t)(i1);
            i1 = updater->size[0];
            ldc_t = (ptrdiff_t)(i1);
            alpha1_t = (double *)(&alpha1);
            Aia0_t = (double *)(&a->data[0]);
            Bib0_t = (double *)(&b->data[0]);
            beta1_t = (double *)(&beta1);
            Cic0_t = (double *)(&C->data[0]);
            dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, alpha1_t, Aia0_t, &lda_t,
                  Bib0_t, &ldb_t, beta1_t, Cic0_t, &ldc_t);
          }
        }
      }
    }

    if (guard1 == TRUE) {
      i1 = C->size[0];
      C->size[0] = a->size[0];
      emxEnsureCapacity(&st, (emxArray__common *)C, i1, (int32_T)sizeof(real_T),
                        &emlrtRTEI);
      loop_ub = a->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        C->data[i1] = 0.0;
        b_loop_ub = a->size[1];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          C->data[i1] += a->data[i1 + a->size[0] * i2] * b->data[i2];
        }
      }
    }

    iv4[0] = r0->size[0];
    emlrtSubAssignSizeCheckR2012b(iv4, 1, *(int32_T (*)[1])C->size, 1, &emlrtECI,
      sp);
    i1 = r1->size[0];
    r1->size[0] = r0->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)r1, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = r0->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r1->data[i1] = r0->data[i1];
    }

    loop_ub = C->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      beta->data[r1->data[i1] + beta->size[0] * ((int32_T)(b_i - 1.0) - 1)] =
        C->data[i1];
    }

    i++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, sp);
  }

  emxFree_int32_T(&r1);
  emxFree_real_T(&b);
  emxFree_real_T(&a);
  emxFree_real_T(&C);
  emxFree_int32_T(&r0);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (BWbetaNloop.c) */
