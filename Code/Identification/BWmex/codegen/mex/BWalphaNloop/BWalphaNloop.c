/*
 * BWalphaNloop.c
 *
 * Code generation for function 'BWalphaNloop'
 *
 * C source code generated on: Mon Jul 07 15:22:55 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWalphaNloop.h"
#include "BWalphaNloop_emxutil.h"

/* Function Definitions */
void BWalphaNloop(emxArray_real_T *alpha, emxArray_real_T *eta, emxArray_real_T *
                  updater)
{
  int32_T i0;
  int32_T i;
  emxArray_real_T *C;
  emxArray_real_T *a;
  emxArray_real_T *b;
  emxArray_real_T *r0;
  emxArray_real_T *b_alpha;
  emxArray_real_T *b_updater;
  int32_T loop_ub;
  int32_T i1;
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
  int32_T exitg1;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);

  /* [alpha,eta]=BWALPHANLOOP(alpha,eta,updater) normalised forward variables for Baum-Welch algorithm */
  /*    alpha    = normalised forward variables */
  /*    eta      = normalisation factor */
  /*    updater  = updater matrices (M*outProj*eta) */
  i0 = eta->size[0];
  i = 1;
  emxInit_real_T(&C, 2, TRUE);
  emxInit_real_T(&a, 2, TRUE);
  emxInit_real_T(&b, 2, TRUE);
  emxInit_real_T(&r0, 2, TRUE);
  emxInit_real_T(&b_alpha, 2, TRUE);
  emxInit_real_T(&b_updater, 2, TRUE);
  while (i - 1 <= (int32_T)((real_T)i0 + -1.0) - 1) {
    loop_ub = alpha->size[1];
    i1 = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)a, i1, (int32_T)sizeof(real_T));
    for (i1 = 0; i1 < loop_ub; i1++) {
      a->data[a->size[0] * i1] = alpha->data[((int32_T)((2.0 + (real_T)(i - 1))
        - 1.0) + alpha->size[0] * i1) - 1];
    }

    loop_ub = updater->size[0];
    b_loop_ub = updater->size[1];
    i1 = b->size[0] * b->size[1];
    b->size[0] = loop_ub;
    b->size[1] = b_loop_ub;
    emxEnsureCapacity((emxArray__common *)b, i1, (int32_T)sizeof(real_T));
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        b->data[i2 + b->size[0] * i1] = updater->data[(i2 + updater->size[0] *
          i1) + updater->size[0] * updater->size[1] * ((int32_T)((2.0 + (real_T)
          (i - 1)) - 1.0) - 1)];
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
        i2 = r0->size[0] * r0->size[1];
        r0->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)r0, i2, (int32_T)sizeof(real_T));
        i2 = r0->size[0] * r0->size[1];
        r0->size[1] = i1;
        emxEnsureCapacity((emxArray__common *)r0, i2, (int32_T)sizeof(real_T));
        loop_ub = i1;
        for (i2 = 0; i2 < loop_ub; i2++) {
          r0->data[i2] = 0.0;
        }

        i2 = C->size[0] * C->size[1];
        C->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)C, i2, (int32_T)sizeof(real_T));
        i2 = C->size[0] * C->size[1];
        C->size[1] = i1;
        emxEnsureCapacity((emxArray__common *)C, i2, (int32_T)sizeof(real_T));
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
      emxEnsureCapacity((emxArray__common *)C, i1, (int32_T)sizeof(real_T));
      loop_ub = b->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        C->data[C->size[0] * i1] = 0.0;
        b_loop_ub = a->size[1];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          C->data[C->size[0] * i1] += a->data[a->size[0] * i2] * b->data[i2 +
            b->size[0] * i1];
        }
      }
    }

    loop_ub = C->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      alpha->data[i + alpha->size[0] * i1] = C->data[C->size[0] * i1];
    }

    i1 = alpha->size[1];
    if (i1 == 0) {
      alpha1 = 0.0;
    } else {
      alpha1 = alpha->data[i];
      loop_ub = 2;
      do {
        exitg1 = 0;
        i1 = alpha->size[1];
        if (loop_ub <= i1) {
          alpha1 += alpha->data[i + alpha->size[0] * (loop_ub - 1)];
          loop_ub++;
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
    }

    eta->data[i] = 1.0 / alpha1;
    loop_ub = alpha->size[1];
    alpha1 = eta->data[i];
    i1 = b_alpha->size[0] * b_alpha->size[1];
    b_alpha->size[0] = 1;
    b_alpha->size[1] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b_alpha, i1, (int32_T)sizeof(real_T));
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_alpha->data[b_alpha->size[0] * i1] = alpha->data[i + alpha->size[0] * i1]
        * alpha1;
    }

    loop_ub = b_alpha->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      alpha->data[i + alpha->size[0] * i1] = b_alpha->data[b_alpha->size[0] * i1];
    }

    loop_ub = updater->size[0];
    b_loop_ub = updater->size[1];
    alpha1 = eta->data[i];
    i1 = b_updater->size[0] * b_updater->size[1];
    b_updater->size[0] = loop_ub;
    b_updater->size[1] = b_loop_ub;
    emxEnsureCapacity((emxArray__common *)b_updater, i1, (int32_T)sizeof(real_T));
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_updater->data[i2 + b_updater->size[0] * i1] = updater->data[(i2 +
          updater->size[0] * i1) + updater->size[0] * updater->size[1] *
          ((int32_T)((2.0 + (real_T)(i - 1)) - 1.0) - 1)] * alpha1;
      }
    }

    loop_ub = b_updater->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_loop_ub = b_updater->size[0];
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        updater->data[(i2 + updater->size[0] * i1) + updater->size[0] *
          updater->size[1] * ((int32_T)((2.0 + (real_T)(i - 1)) - 1.0) - 1)] =
          b_updater->data[i2 + b_updater->size[0] * i1];
      }
    }

    i++;
  }

  emxFree_real_T(&b_updater);
  emxFree_real_T(&b_alpha);
  emxFree_real_T(&r0);
  emxFree_real_T(&b);
  emxFree_real_T(&a);
  emxFree_real_T(&C);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (BWalphaNloop.c) */
