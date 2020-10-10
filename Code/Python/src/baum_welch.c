/* -*- Mode: C -*- */
/* Basic GUFuncs with BLAS
*/
/*
Adapted from https://github.com/numpy/numpy/numpy/linalg/umath_linalg.c.src
Copyright/licence info for that file:
* Copyright (c) 2005-2017, NumPy Developers.
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*   - Redistributions of source code must retain the above
*     copyright notice, this list of conditions and the
*     following disclaimer.
*   - Redistributions in binary form must reproduce the above copyright
*     notice, this list of conditions and the following disclaimer
*     in the documentation and/or other materials provided with the
*     distribution.
*   - Neither the name of the author nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
*****************************************************************************
**                            Includes                                     **
*****************************************************************************
*/
#define NPY_NO_DEPRECATED_API NPY_API_VERSION

#include "rearrange_data.h"

static const char* gufuncs_blas_version = "0.2.0";

/*
*****************************************************************************
**                   Doc string for Python functions                       **
*****************************************************************************
*/

PyDoc_STRVAR(matmul__doc__,
/* "matmul(X: ndarray, Y: ndarray) -> (Z: ndarray)\n\n" */
"Matrix-matrix product.\n\n"
".. deprecated:: 0.2.0\n"
"    This `gufunc` is no longer needed as NumPy switched to a gufunc in v1.16.\n\n"
"Uses BLAS routine `_gemm` for acceleration.\n"
"Does matrix-matrix, matrix-vector, vector-matrix and vector-vector versions,\n"
"with vector versions used *only* when one-dimensional.\n"
"\nParameters\n-----------\n"
"X: ndarray (...,M,N) or (N,)\n"
"    Matrix multiplying from left.\n"
"Y: ndarray (...,N,P) or (N,)\n"
"    Matrix multiplying from right.\n"
"\nReturns\n-------\n"
"Z: ndarray (...,M,P), (...,M), (...,P) or ()\n"
"    Result of matrix multiplication.");

PyDoc_STRVAR(rmatmul__doc__,
/* "rmatmul(Y: ndarray, X: ndarray) -> (Z: ndarray)\n\n" */
"Reversed matrix-matrix product.\n\n"
"Uses BLAS routine `_gemm` for acceleration.\n"
"Does matrix-matrix, matrix-vector, vector-matrix and vector-vector versions,\n"
"with vector versions used *only* when one-dimensional.\n"
"\nParameters\n-----------\n"
"Y: ndarray (...,N,P) or (N,)\n"
"    Matrix multiplying from right.\n"
"X: ndarray (...,M,N) or (N,)\n"
"    Matrix multiplying from left.\n"
"\nReturns\n-------\n"
"Z: ndarray (...,M,P), (...,M), (...,P) or ()\n"
"    Result of matrix multiplication.");

PyDoc_STRVAR(norm__doc__,
/* "norm(X: ndarray) -> (Z: ndarray)\n\n" */
"Euclidean norm of a vector.\n\n"
"Unlike `numpy.linalg.norm`, it only computes the vector 2-norm.\n"
"\nParameters\n-----------\n"
"X: ndarray (...,N)\n"
"    Vector, or array of vectors.\n"
"\nReturns\n-------\n"
"Z: float\n"
"    Euclidean norm of X.");

/*
*****************************************************************************
**                   BLAS/Lapack declarations                              **
*****************************************************************************
*/

/* copy vector x into y */
extern void
FNAME(dcopy)(int *n,
            double *sx, int *incx,
            double *sy, int *incy);

/* x -> sqrt(x'*x) */
extern double
FNAME(dnrm2)(int *n, double sx[], int *inc_x);

/* z -> a x*y + b z */
extern void
FNAME(dgemm)(char tra[], char trb[], int *m, int *n, int *k,
    double *alpha, double a[], int *lda, double b[], int *ldb,
    double *beta, float c[], int *ldc);

/*
******************************************************************************
**                              NORM                                        **
******************************************************************************
*/
/* to hold arguments for BLAS _nrm2 */
typedef struct nrm2_params_struct
{
    /* X is (N,) of base type */
    void *X;

    fortran_int N;
    fortran_int INCX;
} NRM2_PARAMS_t;


/***************************************************
* Calling BLAS/Lapack function _nrm2               *
****************************************************/

static NPY_INLINE double
call_dnrm2(NRM2_PARAMS_t *params)
{
    return BLAS(dnrm2)(&params->N, params->X, &params->INCX);
}

/********************************************************************
* Initialize the parameters to use in for the lapack function _nrm2 *
* Handles buffer allocation
*********************************************************************/
/*initialise parameters for BLAS
 N = length of vector */
static NPY_INLINE int
init_DOUBLE_nrm2(NRM2_PARAMS_t *params, npy_intp N_in)
{
    npy_uint8 *mem_buff = NULL;
    size_t safe_N = N_in;
    fortran_int N = (fortran_int)N_in;

    mem_buff = malloc(safe_N * sizeof(fortran_doublereal));
    if (!mem_buff) {
        goto error;
    }

    params->X = mem_buff;
    params->N = N;
    params->INCX = 1;

    return 1;
 error:
    free(mem_buff);
    memset(params, 0, sizeof(*params));
    PyErr_NoMemory();

    return 0;
}

/************************
* Deallocate buffer     *
*************************/

static NPY_INLINE void
release_DOUBLE_nrm2(NRM2_PARAMS_t *params)
{
    /* memory block base is in X */
    free(params->X);
    memset(params, 0, sizeof(*params));
}

/****************************
* Inner GUfunc loop         *
*****************************/

static void
DOUBLE_norm(char **args, npy_intp *dimensions, npy_intp *steps,
              void *NPY_UNUSED(func))
{
INIT_OUTER_LOOP_2

    /* dimensions of vector */
    npy_intp len_n = *dimensions++;
    /* 1st arg */
    npy_intp stride_n = *steps++;
    NRM2_PARAMS_t params;
    LINEARIZE_DATA_t y_in;

    /* allocate buffer, init params */
    if(init_DOUBLE_nrm2(&params, len_n)) {
        /* init size data */
        init_linearize_vdata(&y_in, len_n, stride_n);

        BEGIN_OUTER_LOOP

            /* copy to buffer */
            linearize_DOUBLE_vec(params.X, args[0], &y_in);
            /* call BLAS */
            *(double *)args[1] = call_dnrm2(&params);

        END_OUTER_LOOP_2
        /* deallocate buffer */
        release_DOUBLE_nrm2(&params);
    }

}

/*
******************************************************************************
**                               MATMUL                                     **
******************************************************************************
*/
/* to hold arguments for BLAS _gemm */
typedef struct gemm_params_struct
{
  void *A; /* A is scalar of base type */
  void *B; /* B is scalar of base type */
  void *X; /* X is (M,K) of base type */
  void *Y; /* Y is (K,N) of base type */
  void *Z; /* Z is (M,N) of base type */

  fortran_int M;
  fortran_int N;
  fortran_int K;
  fortran_int LDX;
  fortran_int LDY;
  fortran_int LDZ;
  char TRANSX;
  char TRANSY;
} GEMM_PARAMS_t;

/****************************************
* Calling BLAS/Lapack function _gemm    *
*****************************************/
/* Z -> A*X.Y + B*Z, need A=1, B=0 */
static NPY_INLINE void
call_dgemm(GEMM_PARAMS_t *params)
{
    BLAS(dgemm)(&params->TRANSX, &params->TRANSY,
       &params->M, &params->N, &params->K,
       params->A,  params->X, &params->LDX, params->Y, &params->LDY,
       params->B, params->Z, &params->LDZ);
}

/********************************************************************
* Initialize the parameters to use in for the lapack function _gemm *
* Handles buffer allocation
*********************************************************************/
/* initialise parameters for BLAS
    M: left dimension
    N: right dimension
    K: inner dimension */
static NPY_INLINE int
init_DOUBLE_matm(GEMM_PARAMS_t *params,
                npy_intp M_in, npy_intp N_in, npy_intp K_in)
{
    npy_uint8 *mem_buff = NULL;
    npy_uint8 *a, *b, *c;
    fortran_int M, N, K, ldx, ldy, ldz;
    size_t safe_M, safe_N, safe_K;

    M = (fortran_int)M_in;
    N = (fortran_int)N_in;
    K = (fortran_int)K_in;
    safe_M = M_in;
    safe_N = N_in;
    safe_K = K_in;
    ldx = fortran_int_max(M, 1);
    ldy = fortran_int_max(K, 1);
    ldz = fortran_int_max(M, 1);

    mem_buff = malloc(safe_M * safe_K * sizeof(fortran_doublereal)
                      + safe_K * safe_N * sizeof(fortran_doublereal)
                      + safe_M * safe_N * sizeof(fortran_doublereal));
    if (!mem_buff) {
        goto error;
    }
    a = mem_buff;  /* X at start of buffer */
    b = a + safe_M * safe_K * sizeof(fortran_doublereal);  /* Y after space for X */
    c = b + safe_K * safe_N * sizeof(fortran_doublereal);  /* Z after space for Y */

    params->TRANSX = 'N';
    params->TRANSY = 'N';
    params->A = &d_one;  /* Z -> A*X.Y + B*Z, need A=1, B=0 */
    params->B = &d_zero;
    params->X = a;
    params->Y = b;
    params->Z = c;
    params->M = M;
    params->N = N;
    params->K = K;
    params->LDX = ldx;
    params->LDY = ldy;
    params->LDZ = ldz;

    return 1;
 error:
    free(mem_buff);
    memset(params, 0, sizeof(*params));
    PyErr_NoMemory();

    return 0;
}

/************************************
* Deallocate buffer                 *
*************************************/

static NPY_INLINE void
release_DOUBLE_matm(GEMM_PARAMS_t *params)
{
    /* memory block base is in X */
    free(params->X);
    memset(params, 0, sizeof(*params));
}

/*****************************
* BLAS Inner GUfunc loop     *
******************************/

/* matmul_signature = "(m,k),(k,n)->(m,n)"; */

static void
DOUBLE_matmul(char **args, npy_intp *dimensions, npy_intp *steps,
              void *NPY_UNUSED(func))
{
INIT_OUTER_LOOP_3

    npy_intp len_m = *dimensions++;  /* dimensions of left */
    npy_intp len_k = *dimensions++;  /* dimensions of inner */
    npy_intp len_n = *dimensions++;  /* dimensions of right */
    npy_intp stride_x_m = *steps++;  /* 1st arg */
    npy_intp stride_x_k = *steps++;
    npy_intp stride_y_k = *steps++;  /* 2nd arg */
    npy_intp stride_y_n = *steps++;
    npy_intp stride_z_m = *steps++;  /* output */
    npy_intp stride_z_n = *steps++;
    GEMM_PARAMS_t params;
    LINEARIZE_DATA_t x_in, y_in, z_out;

    /* allocate buffer */
    if(init_DOUBLE_matm(&params, len_m, len_n, len_k)) {
        /* initialise size parameters */
        init_linearize_data(&x_in, len_k, len_m, stride_x_k, stride_x_m);
        init_linearize_data(&y_in, len_n, len_k, stride_y_n, stride_y_k);
        init_linearize_data(&z_out, len_n, len_m, stride_z_n, stride_z_m);

        BEGIN_OUTER_LOOP
            /* copy inputs */
            linearize_DOUBLE_matrix(params.X, args[0], &x_in);
            linearize_DOUBLE_matrix(params.Y, args[1], &y_in);
            /* call BLAS */
            call_dgemm(&params);
            /* copy output */
            delinearize_DOUBLE_matrix(args[2], params.Z, &z_out);
        END_OUTER_LOOP_3
        /* deallocate buffer */
        release_DOUBLE_matm(&params);
    }
}


/************************************************
*               RMATMUL                         *
*************************************************/

/* rmatmul_signature = "(k,n),(m,k)->(m,n)" */

static void
DOUBLE_rmatmul(char **args, npy_intp *dimensions, npy_intp *steps,
              void *NPY_UNUSED(func))
{
    /*swap (X,Y) for (Y,X)
    x,y here are not the same as X,Y  in docstring: x=Y, y=X
    args = {y, x, z}
    rargs = {x, y, z} */
    char *rargs[] = {args[1], args[0], args[2]};
    /* dimensions = {N, len_k, len_n, len_m};
       rdimensions[] = {N, len_m, len_k, len_n}; */
    npy_intp rdimensions[] = {dimensions[0], dimensions[3], dimensions[1],
                            dimensions[2]};
    /* steps = {strides_y, strides_x, strides_z,
         strides_y_r, strides_y_c,
         strides_x_r, strides_x_c,
         strides_z_r, strides_z_c, };
       rsteps = {strides_x, strides_y, strides_z,
         strides_x_r, strides_x_c,
         strides_y_r, strides_y_c,
         strides_z_r, strides_z_c, }; */
    npy_intp rsteps[] = {steps[1], steps[0], steps[2],
        steps[5], steps[6],
        steps[3], steps[4],
        steps[7], steps[8]};
    /* now that we've swapped a,b and transposed, proceed as if we're in lstsq */
    DOUBLE_matmul(rargs, rdimensions, rsteps, NULL);
}



/*
*****************************************************************************
**                             Ufunc definition                            **
*****************************************************************************
*/

/* types argument for creating 'norm' ufunc */
static char ufn_types_1_3[] = {NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};
static char ufn_types_1_2[] = {NPY_DOUBLE, NPY_DOUBLE};

/* array of functions for each ufunc loop */
static PyUFuncGenericFunction
FUNC_ARRAY_NAME(matmul)[] = {
    DOUBLE_matmul
};

static PyUFuncGenericFunction
FUNC_ARRAY_NAME(rmatmul)[] = {
    DOUBLE_rmatmul
};

static PyUFuncGenericFunction
FUNC_ARRAY_NAME(norm)[] = {
    DOUBLE_norm
};

/* info for creating ufunc object */
GUFUNC_DESCRIPTOR_t gufunc_descriptors[] = {
    {"matmul", "(m?,k),(k,n?)->(m?,n?)", matmul__doc__, 1, 2, 1,
        FUNC_ARRAY_NAME(matmul), ufn_types_1_3},
    {"rmatmul", "(k,n?),(m?,k)->(m?,n?)", rmatmul__doc__, 1, 2, 1,
        FUNC_ARRAY_NAME(rmatmul), ufn_types_1_3},
    {"norm", "(n)->()", norm__doc__, 1, 1, 1,
        FUNC_ARRAY_NAME(norm), ufn_types_1_2}
};

/*
*****************************************************************************
**               Module initialization stuff                               **
*****************************************************************************
*/

/* Methods to add to module (none, we add ufuncs after creating them) */
static PyMethodDef GUfuncs_BLAS_Methods[] = {
    /* Sentinel */
    {NULL, NULL, 0, NULL}
};

/* arguments for module creation */
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_gufuncs_blas",
        NULL,
        -1,
        GUfuncs_BLAS_Methods,
        NULL,
        NULL,
        NULL,
        NULL
};

/* create module */
PyObject *PyInit__gufuncs_blas(void)
{
    PyObject *m;

    init_constants();
    m = PyModule_Create(&moduledef);
    if (m == NULL) {
        return NULL;
    }

    import_array();
    import_ufunc();

    /* Load the ufunc operators into the module's namespace */
    int failure = addUfuncs(m, gufunc_descriptors, 3, gufuncs_blas_version);

    if (PyErr_Occurred() || failure) {
        PyErr_SetString(PyExc_RuntimeError,
                        "cannot load _gufuncs_blas module.");
        return NULL;
    }

    return m;
}

