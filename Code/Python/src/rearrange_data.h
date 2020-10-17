/* -*- Mode: C -*- */
/* Common code for creating GUFuncs with BLAS/Lapack

These functions copy matrices back anf forth between numpy anf fortran forms
Many BLAS/Lapack routines require inputs in (semi)contiguous fortran form
and modify the inputs, so data needs to be copied.
*/
#ifndef GUF_REARRANGE
#define GUF_REARRANGE
/*
*****************************************************************************
**                            Includes                                     **
*****************************************************************************
*/
#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#include "gufunc_common_f.h"
/*
*****************************************************************************
**                            Factories                                    **
*****************************************************************************
Macros for declaring data rearrangement functions
*/

/* Copying from numpy to contiguous fortran */
#define DECLARE_FUNC_LINEARIZE(SHAPE) \
    void linearize_DOUBLE_## SHAPE(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data);

/* Copying from contiguous fortran to numpy */
#define DECLARE_FUNC_DELINEARIZE(SHAPE) \
    void delinearize_DOUBLE_## SHAPE(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data);

/*
*****************************************************************************
**                           Declarations                                  **
*****************************************************************************
*/
DECLARE_FUNC_LINEARIZE(matrix)
DECLARE_FUNC_DELINEARIZE(matrix)
DECLARE_FUNC_LINEARIZE(vec)
DECLARE_FUNC_DELINEARIZE(vec)

#endif
