/*
 * BWalphaNloop.c
 *
 * Code generation for function 'BWalphaNloop'
 *
 * C source code generated on: Fri Jun 27 13:42:24 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWalphaNloop.h"
#include "BWalphaNloop_emxutil.h"
#include "BWalphaNloop_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 7, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRSInfo b_emlrtRSI = { 8, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRSInfo c_emlrtRSI = { 9, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRSInfo d_emlrtRSI = { 10, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRSInfo e_emlrtRSI = { 11, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRSInfo f_emlrtRSI = { 64, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo g_emlrtRSI = { 54, "eml_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m" };

static emlrtRSInfo cb_emlrtRSI = { 1, "mrdivide",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mrdivide.p" };

static emlrtRSInfo db_emlrtRSI = { 15, "rdivide",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/rdivide.m" };

static emlrtRTEInfo emlrtRTEI = { 1, 34, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtBCInfo emlrtBCI = { -1, -1, 8, 24, "alpha", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 8, 45, "updater", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 8, 11, "alpha", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtECInfo emlrtECI = { -1, 8, 5, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtBCInfo d_emlrtBCI = { -1, -1, 9, 28, "alpha", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 10, 24, "alpha", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo f_emlrtBCI = { -1, -1, 10, 31, "eta", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 10, 11, "alpha", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtECInfo b_emlrtECI = { -1, 10, 5, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtBCInfo h_emlrtBCI = { -1, -1, 11, 34, "updater", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo i_emlrtBCI = { -1, -1, 11, 39, "eta", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo j_emlrtBCI = { -1, -1, 11, 17, "updater", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtECInfo c_emlrtECI = { -1, 11, 5, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtBCInfo k_emlrtBCI = { -1, -1, 9, 5, "eta", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

/* Function Definitions */
void BWalphaNloop(const emlrtStack *sp, emxArray_real_T *alpha, emxArray_real_T *
                  eta, emxArray_real_T *updater)
{
  int32_T i0;
  int32_T i;
  emxArray_int32_T *r0;
  emxArray_real_T *C;
  emxArray_int32_T *r1;
  emxArray_real_T *a;
  emxArray_real_T *b;
  emxArray_real_T *r2;
  emxArray_int32_T *r3;
  emxArray_int32_T *r4;
  emxArray_int32_T *r5;
  emxArray_int32_T *r6;
  int32_T i1;
  int32_T k;
  int32_T loop_ub;
  int32_T b_loop_ub;
  int32_T i2;
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
  int32_T iv6[2];
  int32_T exitg1;
  int32_T iv7[2];
  int32_T iv8[2];
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

  /* [alpha,eta]=BWALPHANLOOP(alpha,eta,updater) normalised forward variables for Baum-Welch algorithm */
  /*    alpha    = normalised forward variables */
  /*    eta      = normalisation factor */
  /*    updater  = updater matrices (M*outProj*eta) */
  st.site = &emlrtRSI;
  i0 = eta->size[0];
  i = 0;
  emxInit_int32_T(sp, &r0, 1, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &C, 2, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r1, 1, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &a, 2, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &b, 2, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &r2, 2, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r3, 1, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r4, 1, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r5, 1, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r6, 1, &emlrtRTEI, TRUE);
  while (i <= (int32_T)((real_T)i0 + -1.0) - 1) {
    i1 = alpha->size[0];
    k = (int32_T)(2.0 + (real_T)i);
    emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &c_emlrtBCI, sp);
    loop_ub = alpha->size[1];
    i1 = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r0, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      r0->data[i1] = i1;
    }

    st.site = &b_emlrtRSI;
    loop_ub = alpha->size[1];
    i1 = alpha->size[0];
    k = (int32_T)((2.0 + (real_T)i) - 1.0);
    i1 = emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &emlrtBCI, &st);
    k = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = loop_ub;
    emxEnsureCapacity(&st, (emxArray__common *)a, k, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (k = 0; k < loop_ub; k++) {
      a->data[a->size[0] * k] = alpha->data[(i1 + alpha->size[0] * k) - 1];
    }

    loop_ub = updater->size[0];
    b_loop_ub = updater->size[1];
    i1 = updater->size[2];
    k = (int32_T)((2.0 + (real_T)i) - 1.0);
    i1 = emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &b_emlrtBCI, &st);
    k = b->size[0] * b->size[1];
    b->size[0] = loop_ub;
    b->size[1] = b_loop_ub;
    emxEnsureCapacity(&st, (emxArray__common *)b, k, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (k = 0; k < b_loop_ub; k++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        b->data[i2 + b->size[0] * k] = updater->data[(i2 + updater->size[0] * k)
          + updater->size[0] * updater->size[1] * (i1 - 1)];
      }
    }

    i1 = alpha->size[1];
    guard1 = FALSE;
    if (i1 == 1) {
      guard1 = TRUE;
    } else {
      i1 = updater->size[0];
      if (i1 == 1) {
        guard1 = TRUE;
      } else {
        i1 = updater->size[1];
        k = r2->size[0] * r2->size[1];
        r2->size[0] = 1;
        emxEnsureCapacity(&st, (emxArray__common *)r2, k, (int32_T)sizeof(real_T),
                          &emlrtRTEI);
        k = r2->size[0] * r2->size[1];
        r2->size[1] = i1;
        emxEnsureCapacity(&st, (emxArray__common *)r2, k, (int32_T)sizeof(real_T),
                          &emlrtRTEI);
        loop_ub = i1;
        for (k = 0; k < loop_ub; k++) {
          r2->data[k] = 0.0;
        }

        b_st.site = &f_emlrtRSI;
        c_st.site = &g_emlrtRSI;
        k = C->size[0] * C->size[1];
        C->size[0] = 1;
        emxEnsureCapacity(&c_st, (emxArray__common *)C, k, (int32_T)sizeof
                          (real_T), &emlrtRTEI);
        k = C->size[0] * C->size[1];
        C->size[1] = i1;
        emxEnsureCapacity(&c_st, (emxArray__common *)C, k, (int32_T)sizeof
                          (real_T), &emlrtRTEI);
        loop_ub = i1;
        for (i1 = 0; i1 < loop_ub; i1++) {
          C->data[i1] = 0.0;
        }

        i1 = updater->size[1];
        if (i1 < 1) {
        } else {
          i1 = alpha->size[1];
          if (i1 < 1) {
          } else {
            alpha1 = 1.0;
            beta1 = 0.0;
            TRANSB = 'N';
            TRANSA = 'N';
            m_t = (ptrdiff_t)(1);
            i1 = updater->size[1];
            n_t = (ptrdiff_t)(i1);
            i1 = alpha->size[1];
            k_t = (ptrdiff_t)(i1);
            lda_t = (ptrdiff_t)(1);
            i1 = alpha->size[1];
            ldb_t = (ptrdiff_t)(i1);
            ldc_t = (ptrdiff_t)(1);
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
      i1 = C->size[0] * C->size[1];
      C->size[0] = 1;
      C->size[1] = b->size[1];
      emxEnsureCapacity(&st, (emxArray__common *)C, i1, (int32_T)sizeof(real_T),
                        &emlrtRTEI);
      loop_ub = b->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        C->data[C->size[0] * i1] = 0.0;
        b_loop_ub = a->size[1];
        for (k = 0; k < b_loop_ub; k++) {
          C->data[C->size[0] * i1] += a->data[a->size[0] * k] * b->data[k +
            b->size[0] * i1];
        }
      }
    }

    iv6[0] = 1;
    iv6[1] = r0->size[0];
    emlrtSubAssignSizeCheckR2012b(iv6, 2, *(int32_T (*)[2])C->size, 2, &emlrtECI,
      sp);
    i1 = r3->size[0];
    r3->size[0] = r0->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)r3, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = r0->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r3->data[i1] = r0->data[i1];
    }

    loop_ub = C->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      alpha->data[(i + alpha->size[0] * r3->data[i1]) + 1] = C->data[C->size[0] *
        i1];
    }

    st.site = &c_emlrtRSI;
    i1 = alpha->size[0];
    k = (int32_T)(2.0 + (real_T)i);
    emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &d_emlrtBCI, &st);
    i1 = alpha->size[1];
    if (i1 == 0) {
      alpha1 = 0.0;
    } else {
      alpha1 = alpha->data[i + 1];
      k = 2;
      do {
        exitg1 = 0;
        i1 = alpha->size[1];
        if (k <= i1) {
          alpha1 += alpha->data[(i + alpha->size[0] * (k - 1)) + 1];
          k++;
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
    }

    st.site = &c_emlrtRSI;
    b_st.site = &cb_emlrtRSI;
    c_st.site = &db_emlrtRSI;
    i1 = eta->size[0];
    k = (int32_T)(2.0 + (real_T)i);
    eta->data[emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &k_emlrtBCI, sp) - 1] =
      1.0 / alpha1;
    i1 = alpha->size[0];
    k = (int32_T)(2.0 + (real_T)i);
    emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &g_emlrtBCI, sp);
    loop_ub = alpha->size[1];
    i1 = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r0, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      r0->data[i1] = i1;
    }

    st.site = &d_emlrtRSI;
    i1 = alpha->size[0];
    k = (int32_T)(2.0 + (real_T)i);
    emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &e_emlrtBCI, &st);
    i1 = eta->size[0];
    k = (int32_T)(2.0 + (real_T)i);
    emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &f_emlrtBCI, &st);
    loop_ub = alpha->size[1];
    alpha1 = eta->data[i + 1];
    i1 = r2->size[0] * r2->size[1];
    r2->size[0] = 1;
    r2->size[1] = loop_ub;
    emxEnsureCapacity(&st, (emxArray__common *)r2, i1, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      r2->data[r2->size[0] * i1] = alpha->data[(i + alpha->size[0] * i1) + 1] *
        alpha1;
    }

    iv7[0] = 1;
    iv7[1] = r0->size[0];
    emlrtSubAssignSizeCheckR2012b(iv7, 2, *(int32_T (*)[2])r2->size, 2,
      &b_emlrtECI, sp);
    i1 = r4->size[0];
    r4->size[0] = r0->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)r4, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = r0->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r4->data[i1] = r0->data[i1];
    }

    loop_ub = r2->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      alpha->data[(i + alpha->size[0] * r4->data[i1]) + 1] = r2->data[r2->size[0]
        * i1];
    }

    loop_ub = updater->size[0];
    i1 = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r0, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      r0->data[i1] = i1;
    }

    loop_ub = updater->size[1];
    i1 = r1->size[0];
    r1->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r1, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      r1->data[i1] = i1;
    }

    i1 = updater->size[2];
    k = (int32_T)((2.0 + (real_T)i) - 1.0);
    emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &j_emlrtBCI, sp);
    st.site = &e_emlrtRSI;
    i1 = updater->size[2];
    k = (int32_T)((2.0 + (real_T)i) - 1.0);
    emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &h_emlrtBCI, &st);
    i1 = eta->size[0];
    k = (int32_T)(2.0 + (real_T)i);
    emlrtDynamicBoundsCheckFastR2012b(k, 1, i1, &i_emlrtBCI, &st);
    loop_ub = updater->size[0];
    b_loop_ub = updater->size[1];
    alpha1 = eta->data[i + 1];
    i1 = b->size[0] * b->size[1];
    b->size[0] = loop_ub;
    b->size[1] = b_loop_ub;
    emxEnsureCapacity(&st, (emxArray__common *)b, i1, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (k = 0; k < loop_ub; k++) {
        b->data[k + b->size[0] * i1] = updater->data[(k + updater->size[0] * i1)
          + updater->size[0] * updater->size[1] * ((int32_T)((2.0 + (real_T)i) -
          1.0) - 1)] * alpha1;
      }
    }

    iv8[0] = r0->size[0];
    iv8[1] = r1->size[0];
    emlrtSubAssignSizeCheckR2012b(iv8, 2, *(int32_T (*)[2])b->size, 2,
      &c_emlrtECI, sp);
    i1 = r5->size[0];
    r5->size[0] = r1->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)r5, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = r1->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r5->data[i1] = r1->data[i1];
    }

    i1 = r6->size[0];
    r6->size[0] = r0->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)r6, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = r0->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r6->data[i1] = r0->data[i1];
    }

    loop_ub = b->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_loop_ub = b->size[0];
      for (k = 0; k < b_loop_ub; k++) {
        updater->data[(r6->data[k] + updater->size[0] * r5->data[i1]) +
          updater->size[0] * updater->size[1] * ((int32_T)((2.0 + (real_T)i) -
          1.0) - 1)] = b->data[k + b->size[0] * i1];
      }
    }

    i++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, sp);
  }

  emxFree_int32_T(&r6);
  emxFree_int32_T(&r5);
  emxFree_int32_T(&r4);
  emxFree_int32_T(&r3);
  emxFree_real_T(&r2);
  emxFree_real_T(&b);
  emxFree_real_T(&a);
  emxFree_int32_T(&r1);
  emxFree_real_T(&C);
  emxFree_int32_T(&r0);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (BWalphaNloop.c) */
