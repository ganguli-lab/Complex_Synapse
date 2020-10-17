/* -*- Mode: C -*- */
/* Common code for creating GUFuncs and calling BLAS/Lapack
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
/*         Table of Contents
46.  Includes
63.  Fortran compatibility tools
95.  Structs used for data rearrangement
*/

#ifndef GUF_INCLUDE
#define GUF_INCLUDE

#undef NO_FORTRAN
#define NPY_NO_DEPRECATED_API NPY_API_VERSION
/*
*****************************************************************************
**                              Includes                                   **
*****************************************************************************
Needs to come before typedefs
*/
#include "Python.h"
#include "numpy/ndarraytypes.h"

/*
*****************************************************************************
**                         Fortran types                                   **
*****************************************************************************
*/

/* Aliases for fortran data types */
typedef struct { float r, i; } f2c_complex;
typedef struct { double r, i; } f2c_doublecomplex;

typedef int               fortran_int;
typedef float             fortran_real;
typedef double            fortran_doublereal;
typedef f2c_complex       fortran_complex;
typedef f2c_doublecomplex fortran_doublecomplex;

/* complex number types with choice of interface */
typedef union {
    fortran_complex f;
    npy_cfloat npy;
    float array[2];
} COMPLEX_t;

typedef union {
    fortran_doublecomplex f;
    npy_cdouble npy;
    double array[2];
} DOUBLECOMPLEX_t;

/*
*****************************************************************************
**                              Includes                                   **
*****************************************************************************
Needs to come after typedefs
*/

#define FORTRAN_TYPES 1
#include "gufunc_common.h"

/*
*****************************************************************************
**                         To use BLAS/Lapack                              **
*****************************************************************************
*/

static NPY_INLINE fortran_int
fortran_int_min(fortran_int x, fortran_int y)
{
    return x < y ? x : y;
}

static NPY_INLINE fortran_int
fortran_int_max(fortran_int x, fortran_int y)
{
    return x > y ? x : y;
}

#ifdef NO_APPEND_FORTRAN
# define FNAME(x) x
#else
# define FNAME(x) x##_
#endif

#define BLAS(FUNC)    \
    FNAME(FUNC)

#define LAPACK(FUNC)  \
    FNAME(FUNC)

/*
*****************************************************************************
**               Structs used for data rearrangement                       **
*****************************************************************************
*/

/*
 * this struct contains information about how to linearize a matrix in a local
 * buffer so that it can be used by blas functions.  All strides are specified
 * in bytes and are converted to elements later in type specific functions.
 *
 * rows: number of rows in the matrix
 * columns: number of columns in the matrix
 * row_strides: the number bytes between consecutive rows.
 * column_strides: the number of bytes between consecutive columns.
 * output_lead_dim: BLAS/LAPACK-side leading dimension, in elements
 * conj: should the matrix be complex conjugated?
 *
 * For vectors, we set row/col quantites to the same value.
**/
typedef struct linearize_data_struct
{
    npy_intp rows;
    npy_intp columns;
    npy_intp row_strides;
    npy_intp column_strides;
    npy_intp output_lead_dim;
    npy_intp conj;
} LINEARIZE_DATA_t;

/* Set all parameters */
static NPY_INLINE void
init_linearize_data_exc(LINEARIZE_DATA_t *lin_data,
                        npy_intp rows,
                        npy_intp columns,
                        npy_intp row_strides,
                        npy_intp column_strides,
                        npy_intp output_lead_dim,
                        npy_intp conj)
{
    lin_data->rows = rows;
    lin_data->columns = columns;
    lin_data->row_strides = row_strides;
    lin_data->column_strides = column_strides;
    lin_data->output_lead_dim = output_lead_dim;
    lin_data->conj = conj;
}

/* Set all parameters
assuming no conjugate needed */
static NPY_INLINE void
init_linearize_data_ex(LINEARIZE_DATA_t *lin_data,
                        npy_intp rows,
                        npy_intp columns,
                        npy_intp row_strides,
                        npy_intp column_strides,
                        npy_intp output_lead_dim)
{
    init_linearize_data_exc(lin_data, rows, columns,
                            row_strides, column_strides, output_lead_dim, 0);
}

/* Set all parameters
assuming we use the whole buffer to store matrix on BLAS/Lapack side */
static NPY_INLINE void
init_linearize_datac(LINEARIZE_DATA_t *lin_data,
                    npy_intp rows,
                    npy_intp columns,
                    npy_intp row_strides,
                    npy_intp column_strides,
                    npy_intp conj)
{
    init_linearize_data_exc(lin_data, rows, columns,
                            row_strides, column_strides, columns, conj);
}

/* Set all parameters
assuming no conjugate needed
assuming we use the whole buffer to store matrix on BLAS/Lapack side */
static NPY_INLINE void
init_linearize_data(LINEARIZE_DATA_t *lin_data,
    npy_intp rows,
    npy_intp columns,
    npy_intp row_strides,
    npy_intp column_strides)
{
    init_linearize_datac(lin_data, rows, columns, row_strides, column_strides, 0);
}


/* Set all parameters
for vectors */
static NPY_INLINE void
init_linearize_vdatac(LINEARIZE_DATA_t *lin_data,
                    npy_intp len,
                    npy_intp strides,
                    npy_intp conj)
{
    init_linearize_datac(lin_data, len, len, strides, strides, conj);
}

/* Set all parameters
for vectors
assuming no conjugate needed */
static NPY_INLINE void
init_linearize_vdata(LINEARIZE_DATA_t *lin_data,
                    npy_intp len,
                    npy_intp strides)
{
    init_linearize_vdatac(lin_data, len, strides, 0);
}

#endif
