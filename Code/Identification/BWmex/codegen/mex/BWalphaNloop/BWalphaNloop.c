/*
 * BWalphaNloop.c
 *
 * Code generation for function 'BWalphaNloop'
 *
 * C source code generated on: Wed Jun 25 13:48:35 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWalphaNloop.h"
#include "BWalphaNloop_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "eml_error.h"
#include "BWalphaNloop_mexutil.h"
#include "BWalphaNloop_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 7, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRSInfo b_emlrtRSI = { 9, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRSInfo c_emlrtRSI = { 10, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRSInfo d_emlrtRSI = { 11, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRSInfo e_emlrtRSI = { 12, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRSInfo f_emlrtRSI = { 13, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRSInfo g_emlrtRSI = { 22, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo h_emlrtRSI = { 58, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo i_emlrtRSI = { 61, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo j_emlrtRSI = { 66, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo k_emlrtRSI = { 68, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo l_emlrtRSI = { 41, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtRSInfo m_emlrtRSI = { 56, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtRSInfo n_emlrtRSI = { 86, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtRSInfo o_emlrtRSI = { 24, "indexIntRelop",
  "C:/Program Files/MATLAB/R2013b/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"
};

static emlrtRSInfo q_emlrtRSI = { 117, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtRSInfo r_emlrtRSI = { 233, "indexIntRelop",
  "C:/Program Files/MATLAB/R2013b/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"
};

static emlrtRSInfo s_emlrtRSI = { 16, "max",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/datafun/max.m" };

static emlrtRSInfo t_emlrtRSI = { 18, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtRSInfo u_emlrtRSI = { 88, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtRSInfo w_emlrtRSI = { 59, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

static emlrtRSInfo cb_emlrtRSI = { 215, "indexIntRelop",
  "C:/Program Files/MATLAB/R2013b/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"
};

static emlrtRSInfo db_emlrtRSI = { 9, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtRSInfo eb_emlrtRSI = { 12, "eml_int_forloop_overflow_check",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"
};

static emlrtRSInfo fb_emlrtRSI = { 64, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo gb_emlrtRSI = { 21, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo hb_emlrtRSI = { 54, "eml_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m" };

static emlrtRSInfo ib_emlrtRSI = { 15, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo kb_emlrtRSI = { 32, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo tb_emlrtRSI = { 110, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo ub_emlrtRSI = { 111, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo vb_emlrtRSI = { 112, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo wb_emlrtRSI = { 113, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo xb_emlrtRSI = { 114, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo yb_emlrtRSI = { 115, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo ac_emlrtRSI = { 17, "sum",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtRSInfo bc_emlrtRSI = { 61, "sum",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtRSInfo cc_emlrtRSI = { 8, "isequal",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/isequal.m" };

static emlrtRSInfo dc_emlrtRSI = { 30, "eml_isequal_core",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_isequal_core.m"
};

static emlrtRSInfo ec_emlrtRSI = { 56, "eml_isequal_core",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_isequal_core.m"
};

static emlrtRSInfo fc_emlrtRSI = { 1, "mrdivide",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mrdivide.p" };

static emlrtRSInfo gc_emlrtRSI = { 15, "rdivide",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/rdivide.m" };

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

static emlrtMCInfo m_emlrtMCI = { 18, 9, "sum",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtMCInfo n_emlrtMCI = { 17, 19, "sum",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtMCInfo o_emlrtMCI = { 23, 9, "sum",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtMCInfo p_emlrtMCI = { 20, 19, "sum",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtRTEInfo emlrtRTEI = { 1, 34, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtRTEInfo c_emlrtRTEI = { 65, 1, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtECInfo emlrtECI = { -1, 13, 5, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtBCInfo emlrtBCI = { -1, -1, 13, 17, "updater", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 13, 39, "eta", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 13, 34, "updater", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtECInfo b_emlrtECI = { -1, 12, 5, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtBCInfo d_emlrtBCI = { -1, -1, 12, 11, "alpha", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 12, 31, "eta", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo f_emlrtBCI = { -1, -1, 12, 24, "alpha", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 11, 28, "alpha", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtECInfo c_emlrtECI = { -1, 10, 5, "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m"
};

static emlrtBCInfo h_emlrtBCI = { -1, -1, 10, 11, "alpha", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo i_emlrtBCI = { -1, -1, 10, 53, "updater", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo j_emlrtBCI = { -1, -1, 10, 24, "alpha", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtBCInfo k_emlrtBCI = { -1, -1, 11, 5, "eta", "BWalphaNloop",
  "D:/Users/Subhy/Documents/GitHub/Complex_Synapse/Code/Identification/BWmex/BWalphaNloop.m",
  0 };

static emlrtRSInfo hc_emlrtRSI = { 20, "sum",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtRSInfo ic_emlrtRSI = { 93, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo jc_emlrtRSI = { 98, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo mc_emlrtRSI = { 23, "sum",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtRSInfo nc_emlrtRSI = { 18, "sum",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/datafun/sum.m" };

static emlrtRSInfo oc_emlrtRSI = { 94, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo pc_emlrtRSI = { 99, "mtimes",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo qc_emlrtRSI = { 67, "reshape",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/elmat/reshape.m" };

static emlrtRSInfo rc_emlrtRSI = { 57, "eml_assert_valid_size_arg",
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

void BWalphaNloop(const emlrtStack *sp, emxArray_real_T *alpha, emxArray_real_T *
                  eta, emxArray_real_T *updater)
{
  int32_T i;
  int32_T siz[2];
  int32_T i0;
  int32_T b_i;
  emxArray_int32_T *r0;
  emxArray_real_T *C;
  emxArray_int32_T *r1;
  emxArray_real_T *r2;
  emxArray_real_T *a;
  emxArray_real_T *b;
  emxArray_real_T *r3;
  emxArray_real_T *b_updater;
  emxArray_int32_T *r4;
  emxArray_real_T *b_alpha;
  emxArray_int32_T *r5;
  emxArray_int32_T *r6;
  emxArray_int32_T *r7;
  int32_T i1;
  int32_T loop_ub;
  int32_T nx;
  int32_T k;
  real_T b_a;
  const mxArray *y;
  static const int32_T iv8[2] = { 1, 21 };

  const mxArray *m4;
  char_T cv4[21];
  static const char_T cv5[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'p', 'm', 'a', 'x', 's', 'i', 'z', 'e' };

  int32_T sz[2];
  int32_T b_loop_ub;
  uint32_T varargin_1[2];
  const mxArray *b_y;
  static const int32_T iv9[2] = { 1, 40 };

  char_T cv6[40];
  static const char_T cv7[40] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'g', 'e', 't', 'R', 'e', 's', 'h', 'a', 'p', 'e', 'D',
    'i', 'm', 's', '_', 'n', 'o', 't', 'S', 'a', 'm', 'e', 'N', 'u', 'm', 'e',
    'l' };

  boolean_T b0;
  const mxArray *c_y;
  static const int32_T iv10[2] = { 1, 45 };

  char_T cv8[45];
  static const char_T cv9[45] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'm', 't', 'i', 'm', 'e', 's', '_', 'n', 'o', 'D',
    'y', 'n', 'a', 'm', 'i', 'c', 'S', 'c', 'a', 'l', 'a', 'r', 'E', 'x', 'p',
    'a', 'n', 's', 'i', 'o', 'n' };

  const mxArray *d_y;
  static const int32_T iv11[2] = { 1, 21 };

  static const char_T cv10[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

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
  int32_T iv12[2];
  boolean_T overflow;
  boolean_T p;
  int32_T exitg2;
  const mxArray *e_y;
  static const int32_T iv13[2] = { 1, 30 };

  char_T cv11[30];
  static const char_T cv12[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 's', 'u', 'm', '_', 's', 'p', 'e', 'c', 'i', 'a',
    'l', 'E', 'm', 'p', 't', 'y' };

  boolean_T guard1 = FALSE;
  const mxArray *f_y;
  static const int32_T iv14[2] = { 1, 36 };

  char_T cv13[36];
  static const char_T cv14[36] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'a', 'u', 't', 'o', 'D', 'i', 'm', 'I', 'n', 'c',
    'o', 'm', 'p', 'a', 't', 'i', 'b', 'i', 'l', 'i', 't', 'y' };

  int32_T exitg1;
  int32_T iv15[2];
  int32_T iv16[2];
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &b_st;
  f_st.tls = b_st.tls;
  g_st.prev = &st;
  g_st.tls = st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /* [alpha,eta]=BWALPHANLOOP(alpha,eta,updater) normalised forward variables for Baum-Welch algorithm */
  /*    alpha    = normalised forward variables */
  /*    eta      = normalisation factor */
  /*    updater  = updater matrices (M*outProj*eta) */
  st.site = &emlrtRSI;
  i = alpha->size[1];
  for (i0 = 0; i0 < 2; i0++) {
    siz[i0] = i;
  }

  st.site = &b_emlrtRSI;
  i0 = eta->size[0];
  b_i = 0;
  emxInit_int32_T(sp, &r0, 1, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &C, 2, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r1, 1, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &r2, 2, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &a, 2, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &b, 2, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &r3, 2, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &b_updater, 2, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r4, 1, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &b_alpha, 2, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r5, 1, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r6, 1, &emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r7, 1, &emlrtRTEI, TRUE);
  while (b_i <= (int32_T)((real_T)i0 + -1.0) - 1) {
    i = alpha->size[0];
    i1 = (int32_T)(2.0 + (real_T)b_i);
    emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &h_emlrtBCI, sp);
    loop_ub = alpha->size[1];
    i = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r0, i, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i = 0; i < loop_ub; i++) {
      r0->data[i] = i;
    }

    loop_ub = alpha->size[1];
    i = alpha->size[0];
    i1 = (int32_T)((2.0 + (real_T)b_i) - 1.0);
    i = emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &j_emlrtBCI, sp);
    i1 = a->size[0] * a->size[1];
    a->size[0] = 1;
    a->size[1] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)a, i1, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i1 = 0; i1 < loop_ub; i1++) {
      a->data[a->size[0] * i1] = alpha->data[(i + alpha->size[0] * i1) - 1];
    }

    st.site = &c_emlrtRSI;
    i = updater->size[2];
    i1 = (int32_T)((2.0 + (real_T)b_i) - 1.0);
    emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &i_emlrtBCI, &st);
    i = updater->size[0];
    i1 = updater->size[1];
    nx = i * i1;
    b_st.site = &g_emlrtRSI;
    c_st.site = &l_emlrtRSI;
    c_st.site = &l_emlrtRSI;
    for (k = 0; k < 2; k++) {
      d_st.site = &n_emlrtRSI;
      e_st.site = &o_emlrtRSI;
    }

    c_st.site = &m_emlrtRSI;
    b_a = 1.0;
    for (k = 0; k < 2; k++) {
      if (siz[k] <= 0) {
        b_a = 0.0;
      } else {
        d_st.site = &q_emlrtRSI;
        b_a *= (real_T)siz[k];
      }
    }

    c_st.site = &m_emlrtRSI;
    d_st.site = &o_emlrtRSI;
    e_st.site = &r_emlrtRSI;
    if (2.147483647E+9 >= b_a) {
    } else {
      y = NULL;
      m4 = mxCreateCharArray(2, iv8);
      for (i = 0; i < 21; i++) {
        cv4[i] = cv5[i];
      }

      emlrtInitCharArrayR2013a(&b_st, 21, m4, cv4);
      emlrtAssign(&y, m4);
      c_st.site = &m_emlrtRSI;
      f_st.site = &rc_emlrtRSI;
      error(&c_st, message(&f_st, y, &e_emlrtMCI), &f_emlrtMCI);
    }

    for (i = 0; i < 2; i++) {
      sz[i] = siz[i];
    }

    b_st.site = &h_emlrtRSI;
    loop_ub = updater->size[0];
    b_loop_ub = updater->size[1];
    i = b_updater->size[0] * b_updater->size[1];
    b_updater->size[0] = loop_ub;
    b_updater->size[1] = b_loop_ub;
    emxEnsureCapacity(&b_st, (emxArray__common *)b_updater, i, (int32_T)sizeof
                      (real_T), &emlrtRTEI);
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_updater->data[i1 + b_updater->size[0] * i] = updater->data[(i1 +
          updater->size[0] * i) + updater->size[0] * updater->size[1] *
          ((int32_T)((2.0 + (real_T)b_i) - 1.0) - 1)];
      }
    }

    for (i = 0; i < 2; i++) {
      varargin_1[i] = (uint32_T)b_updater->size[i];
    }

    c_st.site = &s_emlrtRSI;
    d_st.site = &t_emlrtRSI;
    e_st.site = &u_emlrtRSI;
    i = (int32_T)varargin_1[0];
    if ((int32_T)varargin_1[1] > (int32_T)varargin_1[0]) {
      i = (int32_T)varargin_1[1];
    }

    b_st.site = &h_emlrtRSI;
    c_st.site = &s_emlrtRSI;
    d_st.site = &t_emlrtRSI;
    e_st.site = &w_emlrtRSI;
    if (nx < i) {
    } else {
      i = nx;
    }

    if (sz[0] > i) {
      b_st.site = &i_emlrtRSI;
      eml_error(&b_st);
    }

    if (sz[1] > i) {
      b_st.site = &i_emlrtRSI;
      eml_error(&b_st);
    }

    i = b->size[0] * b->size[1];
    b->size[0] = sz[0];
    b->size[1] = sz[1];
    emxEnsureCapacity(&st, (emxArray__common *)b, i, (int32_T)sizeof(real_T),
                      &c_emlrtRTEI);
    b_st.site = &j_emlrtRSI;
    c_st.site = &o_emlrtRSI;
    d_st.site = &cb_emlrtRSI;
    if (nx == sz[0] * sz[1]) {
    } else {
      b_y = NULL;
      m4 = mxCreateCharArray(2, iv9);
      for (i = 0; i < 40; i++) {
        cv6[i] = cv7[i];
      }

      emlrtInitCharArrayR2013a(&st, 40, m4, cv6);
      emlrtAssign(&b_y, m4);
      b_st.site = &j_emlrtRSI;
      g_st.site = &qc_emlrtRSI;
      error(&b_st, message(&g_st, b_y, &emlrtMCI), &b_emlrtMCI);
    }

    b_st.site = &k_emlrtRSI;
    c_st.site = &db_emlrtRSI;
    if (1 > nx) {
      b0 = FALSE;
    } else {
      b0 = (nx > 2147483646);
    }

    if (b0) {
      c_st.site = &eb_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }

    for (k = 0; k + 1 <= nx; k++) {
      loop_ub = updater->size[0];
      b_loop_ub = updater->size[1];
      i = r2->size[0] * r2->size[1];
      r2->size[0] = loop_ub;
      r2->size[1] = b_loop_ub;
      emxEnsureCapacity(&st, (emxArray__common *)r2, i, (int32_T)sizeof(real_T),
                        &emlrtRTEI);
      for (i = 0; i < b_loop_ub; i++) {
        for (i1 = 0; i1 < loop_ub; i1++) {
          r2->data[i1 + r2->size[0] * i] = updater->data[(i1 + updater->size[0] *
            i) + updater->size[0] * updater->size[1] * ((int32_T)((2.0 + (real_T)
            b_i) - 1.0) - 1)];
        }
      }

      i = r2->size[0];
      b->data[k] = r2->data[k % r2->size[0] + r2->size[0] * div_nzp_s32_floor(k,
        i)];
    }

    st.site = &c_emlrtRSI;
    b_st.site = &gb_emlrtRSI;
    i = alpha->size[1];
    if (!(i == b->size[0])) {
      i = alpha->size[1];
      if ((i == 1) || ((b->size[0] == 1) && (b->size[1] == 1))) {
        c_y = NULL;
        m4 = mxCreateCharArray(2, iv10);
        for (i = 0; i < 45; i++) {
          cv8[i] = cv9[i];
        }

        emlrtInitCharArrayR2013a(&b_st, 45, m4, cv8);
        emlrtAssign(&c_y, m4);
        c_st.site = &ic_emlrtRSI;
        f_st.site = &oc_emlrtRSI;
        error(&c_st, message(&f_st, c_y, &i_emlrtMCI), &j_emlrtMCI);
      } else {
        d_y = NULL;
        m4 = mxCreateCharArray(2, iv11);
        for (i = 0; i < 21; i++) {
          cv4[i] = cv10[i];
        }

        emlrtInitCharArrayR2013a(&b_st, 21, m4, cv4);
        emlrtAssign(&d_y, m4);
        c_st.site = &jc_emlrtRSI;
        f_st.site = &pc_emlrtRSI;
        error(&c_st, message(&f_st, d_y, &k_emlrtMCI), &l_emlrtMCI);
      }
    }

    i = alpha->size[1];
    if ((i == 1) || (b->size[0] == 1)) {
      i = C->size[0] * C->size[1];
      C->size[0] = 1;
      C->size[1] = b->size[1];
      emxEnsureCapacity(&st, (emxArray__common *)C, i, (int32_T)sizeof(real_T),
                        &emlrtRTEI);
      loop_ub = b->size[1];
      for (i = 0; i < loop_ub; i++) {
        C->data[C->size[0] * i] = 0.0;
        b_loop_ub = a->size[1];
        for (i1 = 0; i1 < b_loop_ub; i1++) {
          C->data[C->size[0] * i] += a->data[a->size[0] * i1] * b->data[i1 +
            b->size[0] * i];
        }
      }
    } else {
      varargin_1[1] = (uint32_T)b->size[1];
      i = r3->size[0] * r3->size[1];
      r3->size[0] = 1;
      emxEnsureCapacity(&st, (emxArray__common *)r3, i, (int32_T)sizeof(real_T),
                        &emlrtRTEI);
      i = r3->size[0] * r3->size[1];
      r3->size[1] = (int32_T)varargin_1[1];
      emxEnsureCapacity(&st, (emxArray__common *)r3, i, (int32_T)sizeof(real_T),
                        &emlrtRTEI);
      loop_ub = (int32_T)varargin_1[1];
      for (i = 0; i < loop_ub; i++) {
        r3->data[i] = 0.0;
      }

      b_st.site = &fb_emlrtRSI;
      c_st.site = &hb_emlrtRSI;
      i = C->size[0] * C->size[1];
      C->size[0] = 1;
      emxEnsureCapacity(&c_st, (emxArray__common *)C, i, (int32_T)sizeof(real_T),
                        &emlrtRTEI);
      i = C->size[0] * C->size[1];
      C->size[1] = (int32_T)varargin_1[1];
      emxEnsureCapacity(&c_st, (emxArray__common *)C, i, (int32_T)sizeof(real_T),
                        &emlrtRTEI);
      loop_ub = (int32_T)varargin_1[1];
      for (i = 0; i < loop_ub; i++) {
        C->data[i] = 0.0;
      }

      d_st.site = &ib_emlrtRSI;
      if (b->size[1] < 1) {
      } else {
        i = alpha->size[1];
        if (i < 1) {
        } else {
          d_st.site = &kb_emlrtRSI;
          b_a = 1.0;
          beta1 = 0.0;
          TRANSB = 'N';
          TRANSA = 'N';
          e_st.site = &tb_emlrtRSI;
          m_t = (ptrdiff_t)(1);
          e_st.site = &ub_emlrtRSI;
          n_t = (ptrdiff_t)(b->size[1]);
          e_st.site = &vb_emlrtRSI;
          i = alpha->size[1];
          k_t = (ptrdiff_t)(i);
          e_st.site = &wb_emlrtRSI;
          lda_t = (ptrdiff_t)(1);
          e_st.site = &xb_emlrtRSI;
          i = alpha->size[1];
          ldb_t = (ptrdiff_t)(i);
          e_st.site = &yb_emlrtRSI;
          ldc_t = (ptrdiff_t)(1);
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

    iv12[0] = 1;
    iv12[1] = r0->size[0];
    emlrtSubAssignSizeCheckR2012b(iv12, 2, *(int32_T (*)[2])C->size, 2,
      &c_emlrtECI, sp);
    i = r4->size[0];
    r4->size[0] = r0->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)r4, i, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = r0->size[0];
    for (i = 0; i < loop_ub; i++) {
      r4->data[i] = r0->data[i];
    }

    loop_ub = C->size[1];
    for (i = 0; i < loop_ub; i++) {
      alpha->data[(b_i + alpha->size[0] * r4->data[i]) + 1] = C->data[C->size[0]
        * i];
    }

    st.site = &d_emlrtRSI;
    i = alpha->size[0];
    i1 = (int32_T)(2.0 + (real_T)b_i);
    emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &g_emlrtBCI, &st);
    b_st.site = &ac_emlrtRSI;
    c_st.site = &cc_emlrtRSI;
    overflow = FALSE;
    d_st.site = &dc_emlrtRSI;
    e_st.site = &ec_emlrtRSI;
    p = FALSE;
    k = 0;
    do {
      exitg2 = 0;
      if (k < 2) {
        loop_ub = alpha->size[1];
        i = b_alpha->size[0] * b_alpha->size[1];
        b_alpha->size[0] = 1;
        b_alpha->size[1] = loop_ub;
        emxEnsureCapacity(&e_st, (emxArray__common *)b_alpha, i, (int32_T)sizeof
                          (real_T), &emlrtRTEI);
        for (i = 0; i < loop_ub; i++) {
          b_alpha->data[b_alpha->size[0] * i] = alpha->data[(b_i + alpha->size[0]
            * i) + 1];
        }

        if (b_alpha->size[k] != 0) {
          exitg2 = 1;
        } else {
          k++;
        }
      } else {
        p = TRUE;
        exitg2 = 1;
      }
    } while (exitg2 == 0);

    if (!p) {
    } else {
      overflow = TRUE;
    }

    if (!overflow) {
    } else {
      e_y = NULL;
      m4 = mxCreateCharArray(2, iv13);
      for (i = 0; i < 30; i++) {
        cv11[i] = cv12[i];
      }

      emlrtInitCharArrayR2013a(&st, 30, m4, cv11);
      emlrtAssign(&e_y, m4);
      b_st.site = &ac_emlrtRSI;
      g_st.site = &nc_emlrtRSI;
      error(&b_st, message(&g_st, e_y, &m_emlrtMCI), &n_emlrtMCI);
    }

    i = alpha->size[1];
    guard1 = FALSE;
    if (i == 1) {
      guard1 = TRUE;
    } else {
      i = alpha->size[1];
      if (i != 1) {
        guard1 = TRUE;
      } else {
        overflow = FALSE;
      }
    }

    if (guard1 == TRUE) {
      overflow = TRUE;
    }

    if (overflow) {
    } else {
      f_y = NULL;
      m4 = mxCreateCharArray(2, iv14);
      for (i = 0; i < 36; i++) {
        cv13[i] = cv14[i];
      }

      emlrtInitCharArrayR2013a(&st, 36, m4, cv13);
      emlrtAssign(&f_y, m4);
      b_st.site = &hc_emlrtRSI;
      g_st.site = &mc_emlrtRSI;
      error(&b_st, message(&g_st, f_y, &o_emlrtMCI), &p_emlrtMCI);
    }

    i = alpha->size[1];
    if (i == 0) {
      b_a = 0.0;
    } else {
      b_a = alpha->data[b_i + 1];
      b_st.site = &bc_emlrtRSI;
      c_st.site = &db_emlrtRSI;
      i = alpha->size[1];
      if (2 > i) {
        overflow = FALSE;
      } else {
        i = alpha->size[1];
        overflow = (i > 2147483646);
      }

      if (overflow) {
        c_st.site = &eb_emlrtRSI;
        check_forloop_overflow_error(&c_st);
      }

      k = 2;
      do {
        exitg1 = 0;
        i = alpha->size[1];
        if (k <= i) {
          b_a += alpha->data[(b_i + alpha->size[0] * (k - 1)) + 1];
          k++;
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
    }

    st.site = &d_emlrtRSI;
    b_st.site = &fc_emlrtRSI;
    c_st.site = &gc_emlrtRSI;
    i = eta->size[0];
    i1 = (int32_T)(2.0 + (real_T)b_i);
    eta->data[emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &k_emlrtBCI, sp) - 1] =
      1.0 / b_a;
    i = alpha->size[0];
    i1 = (int32_T)(2.0 + (real_T)b_i);
    emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &d_emlrtBCI, sp);
    loop_ub = alpha->size[1];
    i = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r0, i, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i = 0; i < loop_ub; i++) {
      r0->data[i] = i;
    }

    st.site = &e_emlrtRSI;
    i = alpha->size[0];
    i1 = (int32_T)(2.0 + (real_T)b_i);
    emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &f_emlrtBCI, &st);
    i = eta->size[0];
    i1 = (int32_T)(2.0 + (real_T)b_i);
    emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &e_emlrtBCI, &st);
    loop_ub = alpha->size[1];
    b_a = eta->data[b_i + 1];
    i = r3->size[0] * r3->size[1];
    r3->size[0] = 1;
    r3->size[1] = loop_ub;
    emxEnsureCapacity(&st, (emxArray__common *)r3, i, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i = 0; i < loop_ub; i++) {
      r3->data[r3->size[0] * i] = alpha->data[(b_i + alpha->size[0] * i) + 1] *
        b_a;
    }

    iv15[0] = 1;
    iv15[1] = r0->size[0];
    emlrtSubAssignSizeCheckR2012b(iv15, 2, *(int32_T (*)[2])r3->size, 2,
      &b_emlrtECI, sp);
    i = r5->size[0];
    r5->size[0] = r0->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)r5, i, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = r0->size[0];
    for (i = 0; i < loop_ub; i++) {
      r5->data[i] = r0->data[i];
    }

    loop_ub = r3->size[1];
    for (i = 0; i < loop_ub; i++) {
      alpha->data[(b_i + alpha->size[0] * r5->data[i]) + 1] = r3->data[r3->size
        [0] * i];
    }

    loop_ub = updater->size[0];
    i = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r0, i, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i = 0; i < loop_ub; i++) {
      r0->data[i] = i;
    }

    loop_ub = updater->size[1];
    i = r1->size[0];
    r1->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r1, i, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i = 0; i < loop_ub; i++) {
      r1->data[i] = i;
    }

    i = updater->size[2];
    i1 = (int32_T)((2.0 + (real_T)b_i) - 1.0);
    emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &emlrtBCI, sp);
    st.site = &f_emlrtRSI;
    i = updater->size[2];
    i1 = (int32_T)((2.0 + (real_T)b_i) - 1.0);
    emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &c_emlrtBCI, &st);
    i = eta->size[0];
    i1 = (int32_T)(2.0 + (real_T)b_i);
    emlrtDynamicBoundsCheckFastR2012b(i1, 1, i, &b_emlrtBCI, &st);
    loop_ub = updater->size[0];
    b_loop_ub = updater->size[1];
    b_a = eta->data[b_i + 1];
    i = r2->size[0] * r2->size[1];
    r2->size[0] = loop_ub;
    r2->size[1] = b_loop_ub;
    emxEnsureCapacity(&st, (emxArray__common *)r2, i, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        r2->data[i1 + r2->size[0] * i] = updater->data[(i1 + updater->size[0] *
          i) + updater->size[0] * updater->size[1] * ((int32_T)((2.0 + (real_T)
          b_i) - 1.0) - 1)] * b_a;
      }
    }

    iv16[0] = r0->size[0];
    iv16[1] = r1->size[0];
    emlrtSubAssignSizeCheckR2012b(iv16, 2, *(int32_T (*)[2])r2->size, 2,
      &emlrtECI, sp);
    i = r6->size[0];
    r6->size[0] = r1->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)r6, i, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = r1->size[0];
    for (i = 0; i < loop_ub; i++) {
      r6->data[i] = r1->data[i];
    }

    i = r7->size[0];
    r7->size[0] = r0->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)r7, i, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    loop_ub = r0->size[0];
    for (i = 0; i < loop_ub; i++) {
      r7->data[i] = r0->data[i];
    }

    loop_ub = r2->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_loop_ub = r2->size[0];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        updater->data[(r7->data[i1] + updater->size[0] * r6->data[i]) +
          updater->size[0] * updater->size[1] * ((int32_T)((2.0 + (real_T)b_i) -
          1.0) - 1)] = r2->data[i1 + r2->size[0] * i];
      }
    }

    b_i++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, sp);
  }

  emxFree_int32_T(&r7);
  emxFree_int32_T(&r6);
  emxFree_int32_T(&r5);
  emxFree_real_T(&b_alpha);
  emxFree_int32_T(&r4);
  emxFree_real_T(&b_updater);
  emxFree_real_T(&r3);
  emxFree_real_T(&b);
  emxFree_real_T(&a);
  emxFree_real_T(&r2);
  emxFree_int32_T(&r1);
  emxFree_real_T(&C);
  emxFree_int32_T(&r0);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (BWalphaNloop.c) */
