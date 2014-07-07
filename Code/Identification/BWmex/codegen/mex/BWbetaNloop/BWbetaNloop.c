/*
 * BWbetaNloop.c
 *
 * Code generation for function 'BWbetaNloop'
 *
 * C source code generated on: Mon Jul 07 15:23:20 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWbetaNloop.h"
#include "BWbetaNloop_emxutil.h"

/* Function Definitions */
void BWbetaNloop(emxArray_real_T *beta, const emxArray_real_T *updater)
{
  int32_T i0;
  int32_T i;
  emxArray_real_T *C;
  emxArray_real_T *a;
  emxArray_real_T *b;
  real_T b_i;
  int32_T loop_ub;
  int32_T b_loop_ub;
  int32_T i1;
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
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);

  /* beta=BWBETANLOOP(eta,updater) normalised backward variables for Baum-Welch algorithm */
  /*    beta     = normalised backward variables */
  /*    updater  = updater matrices (M*outProj*eta) */
  i0 = beta->size[1];
  i = 0;
  c_emxInit_real_T(&C, 1, TRUE);
  emxInit_real_T(&a, 2, TRUE);
  c_emxInit_real_T(&b, 1, TRUE);
  while (i <= (int32_T)-(2.0 + (-1.0 - (real_T)i0)) - 1) {
    b_i = (real_T)i0 + -(real_T)i;
    loop_ub = updater->size[0];
    b_loop_ub = updater->size[1];
    i1 = a->size[0] * a->size[1];
    a->size[0] = loop_ub;
    a->size[1] = b_loop_ub;
    emxEnsureCapacity((emxArray__common *)a, i1, (int32_T)sizeof(real_T));
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        a->data[i2 + a->size[0] * i1] = updater->data[(i2 + updater->size[0] *
          i1) + updater->size[0] * updater->size[1] * ((int32_T)(b_i - 1.0) - 1)];
      }
    }

    loop_ub = beta->size[0];
    i1 = b->size[0];
    b->size[0] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b, i1, (int32_T)sizeof(real_T));
    for (i1 = 0; i1 < loop_ub; i1++) {
      b->data[i1] = beta->data[i1 + beta->size[0] * ((int32_T)b_i - 1)];
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
        i2 = C->size[0];
        C->size[0] = i1;
        emxEnsureCapacity((emxArray__common *)C, i2, (int32_T)sizeof(real_T));
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
      emxEnsureCapacity((emxArray__common *)C, i1, (int32_T)sizeof(real_T));
      loop_ub = a->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        C->data[i1] = 0.0;
        b_loop_ub = a->size[1];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          C->data[i1] += a->data[i1 + a->size[0] * i2] * b->data[i2];
        }
      }
    }

    loop_ub = C->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      beta->data[i1 + beta->size[0] * ((int32_T)(b_i - 1.0) - 1)] = C->data[i1];
    }

    i++;
  }

  emxFree_real_T(&b);
  emxFree_real_T(&a);
  emxFree_real_T(&C);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (BWbetaNloop.c) */
