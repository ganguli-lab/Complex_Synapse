/*
 * BWbetaNloop_data.c
 *
 * Code generation for function 'BWbetaNloop_data'
 *
 * C source code generated on: Wed Jun 25 13:24:49 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "BWbetaNloop.h"
#include "BWbetaNloop_data.h"

/* Variable Definitions */
const volatile char_T *emlrtBreakCheckR2012bFlagVar;
emlrtRSInfo j_emlrtRSI = { 86, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

emlrtRSInfo l_emlrtRSI = { 239, "indexIntRelop",
  "C:/Program Files/MATLAB/R2013b/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"
};

emlrtRSInfo m_emlrtRSI = { 117, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

emlrtRSInfo n_emlrtRSI = { 233, "indexIntRelop",
  "C:/Program Files/MATLAB/R2013b/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"
};

emlrtRSInfo p_emlrtRSI = { 18, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

emlrtRSInfo q_emlrtRSI = { 88, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

emlrtRSInfo r_emlrtRSI = { 225, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

emlrtRSInfo s_emlrtRSI = { 59, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

emlrtRSInfo t_emlrtRSI = { 124, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

emlrtRSInfo u_emlrtRSI = { 153, "eml_min_or_max",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_min_or_max.m" };

emlrtRSInfo v_emlrtRSI = { 23, "eml_relop",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_relop.m" };

emlrtRSInfo w_emlrtRSI = { 225, "indexIntRelop",
  "C:/Program Files/MATLAB/R2013b/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"
};

emlrtRSInfo x_emlrtRSI = { 215, "indexIntRelop",
  "C:/Program Files/MATLAB/R2013b/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"
};

emlrtRSInfo eb_emlrtRSI = { 15, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

emlrtRSInfo fb_emlrtRSI = { 20, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

emlrtRSInfo gb_emlrtRSI = { 32, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

emlrtRSInfo hb_emlrtRSI = { 63, "eml_refblas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"
};

emlrtRSInfo ib_emlrtRSI = { 85, "eml_refblas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"
};

emlrtRSInfo jb_emlrtRSI = { 89, "eml_refblas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"
};

emlrtRSInfo kb_emlrtRSI = { 90, "eml_refblas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"
};

emlrtRSInfo lb_emlrtRSI = { 94, "eml_refblas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"
};

emlrtRSInfo mb_emlrtRSI = { 96, "eml_refblas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"
};

emlrtRSInfo nb_emlrtRSI = { 110, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

emlrtRSInfo ob_emlrtRSI = { 111, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

emlrtRSInfo pb_emlrtRSI = { 112, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

emlrtRSInfo qb_emlrtRSI = { 113, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

emlrtRSInfo rb_emlrtRSI = { 114, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

emlrtRSInfo sb_emlrtRSI = { 115, "eml_blas_xgemm",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

emlrtMCInfo c_emlrtMCI = { 42, 9, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

emlrtMCInfo d_emlrtMCI = { 41, 19, "eml_assert_valid_size_arg",
  "C:/Program Files/MATLAB/R2013b/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

/* End of code generation (BWbetaNloop_data.c) */
