/* -*- Mode: C -*- */
/* Common code for creating GUFuncs with BLAS/Lapack

These functions copy matrices back anf forth between numpy and fortran forms
Many BLAS/Lapack routines require inputs in (semi)contiguous fortran form
and modify the inputs, so data needs to be copied.
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

/*
*****************************************************************************
**                   BLAS/Lapack calling macros                            **
*****************************************************************************
*/

/* copy vector x into y */
extern void
FNAME(dcopy)(int *n, double *sx, int *incx, double *sy, int *incy);

/*
*****************************************************************************
**                    Data rearrangement functions                         **
*****************************************************************************
*/

              /* rearranging of 2D matrices using blas */

/*********************
*  Copying matrices  *
**********************/
/* Copying full matrix from numpy to contiguous fortran */
void
linearize_DOUBLE_matrix(void *dst_in,
                        const void *src_in,
                        const LINEARIZE_DATA_t* data)
{
    double *src = (double *) src_in;
    double *dst = (double *) dst_in;

    if (dst) {
        int i, j;
        fortran_int columns = (fortran_int)data->columns;
        fortran_int column_strides =
                (fortran_int)(data->column_strides/sizeof(double));
        fortran_int one = 1;
        for (i = 0; i < data->rows; i++) {
            if (column_strides > 0) {
                FNAME(dcopy)(&columns,
                              (void *)src, &column_strides,
                              (void *)dst, &one);
            }
            else if (column_strides < 0) {
                /*
                * BLAS _copy assumes dst points to the first element in memory.
                * numpy points it at the first element of the array.
                * We have to compensate.
                */
                FNAME(dcopy)(  &columns,
                                (void *)(src + (columns-1)*column_strides),
                                &column_strides,
                                (void *)dst, &one);
            }
            else {
                /*
                * Zero stride has undefined behavior in some BLAS
                * implementations (e.g. OSX Accelerate), so do it
                * manually
                */
                for (j = 0; j < columns; ++j) {
                    memcpy((void *)(dst + j), (void *)src, sizeof(double));
                }
            }
            src += data->row_strides/sizeof(double);
            dst += data->output_lead_dim;
        }
    }
}

/* Copying full matrix from contiguous fortran to numpy */
void
delinearize_DOUBLE_matrix(void *dst_in,
                        const void *src_in,
                        const LINEARIZE_DATA_t* data)
{
    double *src = (double *) src_in;
    double *dst = (double *) dst_in;

    if (src) {
        int i;
        fortran_int columns = (fortran_int)data->columns;
        fortran_int column_strides =
            (fortran_int)(data->column_strides/sizeof(double));
        fortran_int one = 1;
        for (i = 0; i < data->rows; i++) {
            if (column_strides > 0) {
                FNAME(dcopy)(  &columns,
                                (void *)src, &one,
                                (void *)dst, &column_strides);
            }
            else if (column_strides < 0) {
                /*
                * BLAS _copy/lacgv assume dst points to first element in memory
                * numpy points it at the first element of the array.
                * We have to compensate.
                */
                double *dst_end = dst + (columns-1)*column_strides;
                FNAME(dcopy)(  &columns,
                                (void *)src, &one,
                                (void *)dst_end, &column_strides);
            }
            else {
                /*
                * Zero stride has undefined behavior in some BLAS
                * implementations (e.g. OSX Accelerate), so do it
                * manually
                */
                if (columns > 0) {
                    memcpy((void *)dst, (void *)(src + (columns-1)), sizeof(double));
                }
            }
            src += data->output_lead_dim;
            dst += data->row_strides/sizeof(double);
        }
    }
}

/*********************
*      Vectors       *
**********************/

/* Copying vector from numpy to contiguous fortran */
void
linearize_DOUBLE_vec(void *dst_in,
                    const void *src_in,
                    const LINEARIZE_DATA_t *data)
{
    double *src = (double *) src_in;
    double *dst = (double *) dst_in;

    if (dst) {
        fortran_int len = (fortran_int)data->columns;
        fortran_int strides = (fortran_int)(data->column_strides/sizeof(double));
        fortran_int one = 1;
        if (strides > 0) {
            FNAME(dcopy)(  &len,
                            (void *)src, &strides,
                            (void *)dst, &one);
        }
        else if (strides < 0) {
            /*
            * Lapack _copy assumes dst points to first element in memory
            * instead of first element of array & tries to compensate.
            * We have to undo that
            */
            FNAME(dcopy)(  &len,
                            (void *)(src + (len-1)*strides), &strides,
                            (void *)dst, &one);
        }
        else {
            /*
            * Zero stride has undefined behavior in some BLAS
            * implementations (e.g. OSX Accelerate), so do it
            * manually
            */
            int j;
            for (j = 0; j < len; ++j) {
                memcpy((void *)(dst + j), (void *)src, sizeof(double));
            }
        }
    }
}

/* Copying vector from contiguous fortran to numpy */
void
delinearize_DOUBLE_vec(void *dst_in,
                     const void *src_in,
                     const LINEARIZE_DATA_t *data)
{
    double *src = (double *) src_in;
    double *dst = (double *) dst_in;

    if (dst) {
        fortran_int len = (fortran_int)data->columns;
        fortran_int strides = (fortran_int)(data->column_strides/sizeof(double));
        fortran_int one = 1;
        if (strides > 0) {
            FNAME(dcopy)(  &len,
                            (void *)src, &one,
                            (void *)dst, &strides);
        }
        else if (strides < 0) {
            /*
            * BLAS _copy assumes dst points to the first element in memory.
            * numpy points it at the first element of the array.
            * We have to compensate.
            */
            FNAME(dcopy)(  &len,
                            (void *)src, &one,
                            (void *)(dst + (len-1)*strides), &strides);
        }
        else {
            /*
            * Zero stride has undefined behavior in some BLAS
            * implementations (e.g. OSX Accelerate), so do it
            * manually
            */
            if (len > 0) {
                int j;
                for (j = 0; j < len; ++j) {
                    memcpy((void *)dst, (void *)(src + j), sizeof(double));
                }
            }
        }
    }
}
