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
#define DECLARE_FUNC_LINEARIZE(SHAPE)                                                               \
    void linearize_FLOAT_## SHAPE(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data);  \
    void linearize_DOUBLE_## SHAPE(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data); \
    void linearize_CFLOAT_## SHAPE(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data); \
    void linearize_CDOUBLE_## SHAPE(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data);

/* Copying from contiguous fortran to numpy */
#define DECLARE_FUNC_DELINEARIZE(SHAPE)                                                               \
    void delinearize_FLOAT_## SHAPE(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data);  \
    void delinearize_DOUBLE_## SHAPE(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data); \
    void delinearize_CFLOAT_## SHAPE(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data); \
    void delinearize_CDOUBLE_## SHAPE(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data);

/* Fill a numpyy matrix with some value */
#define DECLARE_FUNC_FILL(NAME, SHAPE)                                        \
    void NAME ##_FLOAT_## SHAPE(void *dst_in, const LINEARIZE_DATA_t* data);  \
    void NAME ##_DOUBLE_## SHAPE(void *dst_in, const LINEARIZE_DATA_t* data); \
    void NAME ##_CFLOAT_## SHAPE(void *dst_in, const LINEARIZE_DATA_t* data); \
    void NAME ##_CDOUBLE_## SHAPE(void *dst_in, const LINEARIZE_DATA_t* data);
/*
*****************************************************************************
**                           Declarations                                  **
*****************************************************************************
*/
DECLARE_FUNC_LINEARIZE(matrix)
DECLARE_FUNC_DELINEARIZE(matrix)
DECLARE_FUNC_DELINEARIZE(triu)
DECLARE_FUNC_DELINEARIZE(tril)
DECLARE_FUNC_DELINEARIZE(diag)
DECLARE_FUNC_FILL(nan, matrix)
DECLARE_FUNC_FILL(zero, matrix)
DECLARE_FUNC_FILL(zero, triu)
DECLARE_FUNC_FILL(zero, tril)
DECLARE_FUNC_FILL(eye, matrix)
DECLARE_FUNC_LINEARIZE(vec)
DECLARE_FUNC_DELINEARIZE(vec)
DECLARE_FUNC_FILL(nan, vec)

/* combination of delinearize_ and zero_, _triu and _tril */
void delinearize_FLOAT_trilu(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data, int lower);
void delinearize_DOUBLE_trilu(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data, int lower);
void delinearize_CFLOAT_trilu(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data, int lower);
void delinearize_CDOUBLE_trilu(void *dst_in, const void *src_in, const LINEARIZE_DATA_t* data, int lower);

/* These functions convert floating point variables to integers */
fortran_int FLOAT_real_int(fortran_real val);
fortran_int DOUBLE_real_int(fortran_doublereal val);
fortran_int CFLOAT_real_int(fortran_complex val);
fortran_int CDOUBLE_real_int(fortran_doublecomplex val);
/*
*****************************************************************************
**                          Integer versions                               **
*****************************************************************************
*/
/* Copying from numpy to contiguous fortran */
static NPY_INLINE void
linearize_INT_vec(void *dst_in,
                const void *src_in,
                const LINEARIZE_DATA_t *data)
{
    fortran_int *dst = (fortran_int *) dst_in;
    npy_int *src = (npy_int *) src_in;

    if (dst) {
        fortran_int len = (fortran_int)data->columns;
        fortran_int strides = (fortran_int)(data->column_strides/sizeof(npy_int));
        int j;
        for (j = 0; j < len; ++j) {
            *dst = (fortran_int)*src;
            src += strides;
            dst += 1;
        }
    }
}

/* Copying from contiguous fortran to numpy */
static NPY_INLINE void
delinearize_INT_vec(void *dst_in,
                    const void *src_in,
                    const LINEARIZE_DATA_t *data)
{
    fortran_int *src = (fortran_int *) src_in;
    npy_int *dst = (npy_int *) dst_in;

    if (dst) {
        fortran_int len = (fortran_int)data->columns;
        fortran_int strides = (fortran_int)(data->column_strides/sizeof(npy_int));
        int j;
        for (j = 0; j < len; ++j) {
            *dst = (npy_int)*src;
            src += 1;
            dst += strides;
        }
    }
}
/*
*****************************************************************************
**                          Singularity checks                             **
*****************************************************************************
*/

void FNAME(strchk)(int *n, float *a, int *lda, int *info);
void FNAME(dtrchk)(int *n, double *a, int *lda, int *info);
void FNAME(ctrchk)(int *n, f2c_complex *a, int *lda, int *info);
void FNAME(ztrchk)(int *n, f2c_doublecomplex *a, int *lda, int *info);


#endif
