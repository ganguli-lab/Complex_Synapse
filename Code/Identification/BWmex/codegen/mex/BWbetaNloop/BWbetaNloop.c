/*
 * BWbetaNloop.c
 *
 * Code generation for function 'BWbetaNloop'
 *
 * C source code generated on: Wed Jun 25 13:24:49 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWbetaNloop.h"
#include "BWbetaNloop_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "eml_error.h"
#include "BWbetaNloop_mexutil.h"
#include "BWbetaNloop_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 7, "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m"
};

static emlrtRSInfo b_emlrtRSI = { 10, "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m"
};

static emlrtRSInfo c_emlrtRSI = { 22, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo d_emlrtRSI = { 58, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo e_emlrtRSI = { 61, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo f_emlrtRSI = { 66, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo g_emlrtRSI = { 68, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo h_emlrtRSI = { 41, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtRSInfo i_emlrtRSI = { 56, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtRSInfo k_emlrtRSI = { 24, "indexIntRelop",
  "C:/Program Files/MATLAB/R2013b/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"
};

static emlrtRSInfo o_emlrtRSI = { 16, "max",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/datafun/max.m" };

static emlrtRSInfo y_emlrtRSI = { 9, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtRSInfo ab_emlrtRSI = { 12, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtRSInfo bb_emlrtRSI = { 64, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo cb_emlrtRSI = { 21, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo db_emlrtRSI = { 54, "eml_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m" };

static emlrtMCInfo emlrtMCI = { 67, 5, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtMCInfo b_emlrtMCI = { 66, 15, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtMCInfo e_emlrtMCI = { 57, 5, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtMCInfo f_emlrtMCI = { 56, 15, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtMCInfo i_emlrtMCI = { 94, 13, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo j_emlrtMCI = { 93, 23, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo k_emlrtMCI = { 99, 13, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo l_emlrtMCI = { 98, 23, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRTEInfo emlrtRTEI = { 1, 21, "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m"
};

static emlrtRTEInfo c_emlrtRTEI = { 65, 1, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtECInfo emlrtECI = { -1, 10, 5, "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m"
};

static emlrtBCInfo emlrtBCI = { -1, -1, 10, 12, "beta", "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m",
  0 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 10, 58, "beta", "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m",
  0 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 10, 39, "updater", "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m",
  0 };

static emlrtRTEInfo e_emlrtRTEI = { 9, 1, "BWbetaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWbetaNloop.m"
};

static emlrtRSInfo tb_emlrtRSI = { 93, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo ub_emlrtRSI = { 98, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo wb_emlrtRSI = { 94, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo xb_emlrtRSI = { 99, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo yb_emlrtRSI = { 67, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo ac_emlrtRSI = { 57, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

/* Function Declarations */
static int32_T div_nzp_s32_floor(int32_T numerator, int32_T denominator);
static const mxArray *message(const emlrtStack *sp, const mxArray *b,
  emlrtMCInfo *location);

/* Function Definitions */
static int32_T div_nzp_s32_floor(int32_T numerator, int32_T denominator)
{
  int32_T quotient;
  uint32_T absNumerator;
  uint32_T absDenominator;
  int32_T quotientNeedsNegation;
  uint32_T tempAbsQuotient;
  if (numerator >= 0) {
    absNumerator = (uint32_T)numerator;
  } else {
    absNumerator = (uint32_T)-numerator;
  }

  if (denominator >= 0) {
    absDenominator = (uint32_T)denominator;
  } else {
    absDenominator = (uint32_T)-denominator;
  }

  quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
  tempAbsQuotient = absNumerator / absDenominator;
  if ((uint32_T)quotientNeedsNegation) {
    absNumerator %= absDenominator;
    if (absNumerator > (uint32_T)0) {
      tempAbsQuotient++;
    }
  }

  if ((uint32_T)quotientNeedsNegation) {
    quotient = -(int32_T)tempAbsQuotient;
  } else {
    quotient = (int32_T)tempAbsQuotient;
  }

  return quotient;
}

static const mxArray *message(const emlrtStack *sp, const mxArray *b,
  emlrtMCInfo *location)
{
  const mxArray *pArray;
  const mxArray *m1;
  pArray = b;
  return emlrtCallMATLABR2012b(sp, 1, &m1, 1, &pArray, "message", TRUE, location);
}

void BWbetaNloop(const emlrtStack *sp, emxArray_real_T *beta, const
                 emxArray_real_T *updater)
{
  int32_T i;
  int32_T siz[2];
  int32_T i0;
  int32_T b_i;
  emxArray_int32_T *r0;
  emxArray_real_T *C;
  emxArray_real_T *a;
  emxArray_real_T *b;
  emxArray_real_T *r1;
  emxArray_real_T *b_updater;
  emxArray_int32_T *r2;
  real_T c_i;
  int32_T loop_ub;
  int32_T i1;
  int32_T nx;
  int32_T k;
  real_T b_a;
  const mxArray *y;
  static const int32_T iv6[2] = { 1, 21 };

  const mxArray *m4;
  char_T cv4[21];
  static const char_T cv5[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'p', 'm', 'a', 'x', 's', 'i', 'z', 'e' };

  int32_T sz[2];
  int32_T b_loop_ub;
  uint32_T varargin_1[2];
  const mxArray *b_y;
  static const int32_T iv7[2] = { 1, 40 };

  char_T cv6[40];
  static const char_T cv7[40] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'g', 'e', 't', 'R', 'e', 's', 'h', 'a', 'p', 'e', 'D',
    'i', 'm', 's', '_', 'n', 'o', 't', 'S', 'a', 'm', 'e', 'N', 'u', 'm', 'e',
    'l' };

  boolean_T b0;
  boolean_T guard2 = FALSE;
  const mxArray *c_y;
  static const int32_T iv8[2] = { 1, 21 };

  static const char_T cv8[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

  const mxArray *d_y;
  static const int32_T iv9[2] = { 1, 45 };

  char_T cv9[45];
  static const char_T cv10[45] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'm', 't', 'i', 'm', 'e', 's', '_', 'n', 'o', 'D',
    'y', 'n', 'a', 'm', 'i', 'c', 'S', 'c', 'a', 'l', 'a', 'r', 'E', 'x', 'p',
    'a', 'n', 's', 'i', 'o', 'n' };

  boolean_T guard1 = FALSE;
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
  int32_T iv10[1];
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &b_st;
  d_st.tls = b_st.tls;
  e_st.prev = &st;
  e_st.tls = st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /* beta=BWBETANLOOP(eta,updater) normalised backward variables for Baum-Welch algorithm */
  /*    beta     = normalised backward variables */
  /*    updater  = updater matrices (M*outProj*eta) */
  st.site = &emlrtRSI;
  i = beta->size[0];
  for (i0 = 0; i0 < 2; i0++) {
    siz[i0] = i;
  }

  i0 = beta->size[1];
  emlrtForLoopVectorCheckR2012b(beta->size[1], -1.0, 2.0, mxDOUBLE_CLASS,
    (int32_T)-(2.0 + (-1.0 - (real_T)beta->size[1])), &e_emlrtRTEI, sp);
  b_i = 0;
  emxInit_int32_T(sp, &r0, 1, &emlrtRTEI, TRUE);
  c_emxInit_real_T(sp, &C, 1, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &a, 2, &emlrtRTEI, TRUE);
  c_emxInit_real_T(sp, &b, 1, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &r1, 2, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &b_updater, 2, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r2, 1, &emlrtRTEI, TRUE);
  while (b_i <= (int32_T)-(2.0 + (-1.0 - (real_T)i0)) - 1) {
    c_i = (real_T)i0 + -(real_T)b_i;
    loop_ub = beta->size[0];
    i1 = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r0, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      r0->data[i1] = i1;
    }

    i1 = beta->size[1];
    i = (int32_T)(c_i - 1.0);
    emlrtDynamicBoundsCheckFastR2012b(i, 1, i1, &emlrtBCI, sp);
    st.site = &b_emlrtRSI;
    i1 = updater->size[2];
    i = (int32_T)(c_i - 1.0);
    emlrtDynamicBoundsCheckFastR2012b(i, 1, i1, &c_emlrtBCI, &st);
    i1 = updater->size[0];
    i = updater->size[1];
    nx = i1 * i;
    b_st.site = &c_emlrtRSI;
    c_st.site = &h_emlrtRSI;
    c_st.site = &h_emlrtRSI;
    c_st.site = &i_emlrtRSI;
    b_a = 1.0;
    for (k = 0; k < 2; k++) {
      if (siz[k] <= 0) {
        b_a = 0.0;
      } else {
        b_a *= (real_T)siz[k];
      }
    }

    c_st.site = &i_emlrtRSI;
    if (2.147483647E+9 >= b_a) {
    } else {
      y = NULL;
      m4 = mxCreateCharArray(2, iv6);
      for (i = 0; i < 21; i++) {
        cv4[i] = cv5[i];
      }

      emlrtInitCharArrayR2013a(&b_st, 21, m4, cv4);
      emlrtAssign(&y, m4);
      c_st.site = &i_emlrtRSI;
      d_st.site = &ac_emlrtRSI;
      error(&c_st, message(&d_st, y, &e_emlrtMCI), &f_emlrtMCI);
    }

    for (i1 = 0; i1 < 2; i1++) {
      sz[i1] = siz[i1];
    }

    b_st.site = &d_emlrtRSI;
    loop_ub = updater->size[0];
    b_loop_ub = updater->size[1];
    i1 = b_updater->size[0] * b_updater->size[1];
    b_updater->size[0] = loop_ub;
    b_updater->size[1] = b_loop_ub;
    emxEnsureCapacity(&b_st, (emxArray__common *)b_updater, i1, (int32_T)sizeof
                      (real_T), &emlrtRTEI);
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      for (i = 0; i < loop_ub; i++) {
        b_updater->data[i + b_updater->size[0] * i1] = updater->data[(i +
          updater->size[0] * i1) + updater->size[0] * updater->size[1] *
          ((int32_T)(c_i - 1.0) - 1)];
      }
    }

    for (i1 = 0; i1 < 2; i1++) {
      varargin_1[i1] = (uint32_T)b_updater->size[i1];
    }

    c_st.site = &o_emlrtRSI;
    i = (int32_T)varargin_1[0];
    if ((int32_T)varargin_1[1] > (int32_T)varargin_1[0]) {
      i = (int32_T)varargin_1[1];
    }

    b_st.site = &d_emlrtRSI;
    c_st.site = &o_emlrtRSI;
    if (nx < i) {
    } else {
      i = nx;
    }

    if (sz[0] > i) {
      b_st.site = &e_emlrtRSI;
      eml_error(&b_st);
    }

    if (sz[1] > i) {
      b_st.site = &e_emlrtRSI;
      eml_error(&b_st);
    }

    i1 = a->size[0] * a->size[1];
    a->size[0] = sz[0];
    a->size[1] = sz[1];
    emxEnsureCapacity(&st, (emxArray__common *)a, i1, (int32_T)sizeof(real_T),
                      &c_emlrtRTEI);
    b_st.site = &f_emlrtRSI;
    c_st.site = &k_emlrtRSI;
    if (nx == sz[0] * sz[1]) {
    } else {
      b_y = NULL;
      m4 = mxCreateCharArray(2, iv7);
      for (i = 0; i < 40; i++) {
        cv6[i] = cv7[i];
      }

      emlrtInitCharArrayR2013a(&st, 40, m4, cv6);
      emlrtAssign(&b_y, m4);
      b_st.site = &f_emlrtRSI;
      e_st.site = &yb_emlrtRSI;
      error(&b_st, message(&e_st, b_y, &emlrtMCI), &b_emlrtMCI);
    }

    b_st.site = &g_emlrtRSI;
    c_st.site = &y_emlrtRSI;
    if (1 > nx) {
      b0 = FALSE;
    } else {
      b0 = (nx > 2147483646);
    }

    if (b0) {
      c_st.site = &ab_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }

    for (k = 0; k + 1 <= nx; k++) {
      loop_ub = updater->size[0];
      b_loop_ub = updater->size[1];
      i1 = r1->size[0] * r1->size[1];
      r1->size[0] = loop_ub;
      r1->size[1] = b_loop_ub;
      emxEnsureCapacity(&st, (emxArray__common *)r1, i1, (int32_T)sizeof(real_T),
                        &emlrtRTEI);
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        for (i = 0; i < loop_ub; i++) {
          r1->data[i + r1->size[0] * i1] = updater->data[(i + updater->size[0] *
            i1) + updater->size[0] * updater->size[1] * ((int32_T)(c_i - 1.0) -
            1)];
        }
      }

      i1 = r1->size[0];
      a->data[k] = r1->data[k % r1->size[0] + r1->size[0] * div_nzp_s32_floor(k,
        i1)];
    }

    loop_ub = beta->size[0];
    i1 = beta->size[1];
    i = (int32_T)c_i;
    i = emlrtDynamicBoundsCheckFastR2012b(i, 1, i1, &b_emlrtBCI, sp);
    i1 = b->size[0];
    b->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)b, i1, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      b->data[i1] = beta->data[i1 + beta->size[0] * (i - 1)];
    }

    st.site = &b_emlrtRSI;
    b_st.site = &cb_emlrtRSI;
    i1 = beta->size[0];
    if (!(a->size[1] == i1)) {
      guard2 = FALSE;
      if ((a->size[0] == 1) && (a->size[1] == 1)) {
        guard2 = TRUE;
      } else {
        i1 = beta->size[0];
        if (i1 == 1) {
          guard2 = TRUE;
        } else {
          c_y = NULL;
          m4 = mxCreateCharArray(2, iv8);
          for (i = 0; i < 21; i++) {
            cv4[i] = cv8[i];
          }

          emlrtInitCharArrayR2013a(&b_st, 21, m4, cv4);
          emlrtAssign(&c_y, m4);
          c_st.site = &ub_emlrtRSI;
          d_st.site = &xb_emlrtRSI;
          error(&c_st, message(&d_st, c_y, &k_emlrtMCI), &l_emlrtMCI);
        }
      }

      if (guard2 == TRUE) {
        d_y = NULL;
        m4 = mxCreateCharArray(2, iv9);
        for (i = 0; i < 45; i++) {
          cv9[i] = cv10[i];
        }

        emlrtInitCharArrayR2013a(&b_st, 45, m4, cv9);
        emlrtAssign(&d_y, m4);
        c_st.site = &tb_emlrtRSI;
        d_st.site = &wb_emlrtRSI;
        error(&c_st, message(&d_st, d_y, &i_emlrtMCI), &j_emlrtMCI);
      }
    }

    guard1 = FALSE;
    if (a->size[1] == 1) {
      guard1 = TRUE;
    } else {
      i1 = beta->size[0];
      if (i1 == 1) {
        guard1 = TRUE;
      } else {
        varargin_1[0] = (uint32_T)a->size[0];
        b_st.site = &bb_emlrtRSI;
        c_st.site = &db_emlrtRSI;
        i1 = C->size[0];
        C->size[0] = (int32_T)varargin_1[0];
        emxEnsureCapacity(&c_st, (emxArray__common *)C, i1, (int32_T)sizeof
                          (real_T), &emlrtRTEI);
        loop_ub = (int32_T)varargin_1[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          C->data[i1] = 0.0;
        }

        if ((a->size[0] < 1) || (a->size[1] < 1)) {
        } else {
          b_a = 1.0;
          beta1 = 0.0;
          TRANSB = 'N';
          TRANSA = 'N';
          m_t = (ptrdiff_t)(a->size[0]);
          n_t = (ptrdiff_t)(1);
          k_t = (ptrdiff_t)(a->size[1]);
          lda_t = (ptrdiff_t)(a->size[0]);
          ldb_t = (ptrdiff_t)(a->size[1]);
          ldc_t = (ptrdiff_t)(a->size[0]);
          alpha1_t = (double *)(&b_a);
          Aia0_t = (double *)(&a->data[0]);
          Bib0_t = (double *)(&b->data[0]);
          beta1_t = (double *)(&beta1);
          Cic0_t = (double *)(&C->data[0]);
          dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, alpha1_t, Aia0_t, &lda_t,
                Bib0_t, &ldb_t, beta1_t, Cic0_t, &ldc_t);
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
        for (i = 0; i < b_loop_ub; i++) {
          C->data[i1] += a->data[i1 + a->size[0] * i] * b->data[i];
        }
      }
    }

    iv10[0] = r0->size[0];
    emlrtSubAssignSizeCheckR2012b(iv10, 1, *(int32_T (*)[1])C->size, 1,
      &emlrtECI, sp);
    i1 = r2->size[0];
    r2->size[0] = r0->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)r2, i1, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = r0->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r2->data[i1] = r0->data[i1];
    }

    loop_ub = C->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      beta->data[r2->data[i1] + beta->size[0] * ((int32_T)(c_i - 1.0) - 1)] =
        C->data[i1];
    }

    b_i++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, sp);
  }

  emxFree_int32_T(&r2);
  emxFree_real_T(&b_updater);
  emxFree_real_T(&r1);
  emxFree_real_T(&b);
  emxFree_real_T(&a);
  emxFree_real_T(&C);
  emxFree_int32_T(&r0);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (BWbetaNloop.c) */
