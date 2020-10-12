/* -*- Mode: C -*- */
/* Common code for creating GUFuncs
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

#ifndef GUC_INCLUDE
#define GUC_INCLUDE

/*
*****************************************************************************
**                             Includes                                    **
*****************************************************************************
*/

#define NPY_NO_DEPRECATED_API NPY_API_VERSION

#include "Python.h"
#include "numpy/ndarraytypes.h"
#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"
#include "numpy/npy_math.h"
#include "numpy/npy_3kcompat.h"
/* #include "npy_config.h" */

/*
*****************************************************************************
**                         Outer loop macros                               **
*****************************************************************************
Macros for looping over non-core dimensions provided by numpy
*/

/* Declare and initialise loop variables:
    dN = length of non-core dimension
    N_ = loop variable
    s0 = stride along non-core dimension of first argument
When number of ufunc inputs + outputs == 1 */
#define INIT_OUTER_LOOP_1         \
    npy_intp dN = *dimensions++;  \
    npy_intp N_;                  \
    npy_intp s0 = *steps++;

/* Declare and initialise loop variables:
    s1 = stride along non-core dimension of second argument
When number of ufunc inputs + outputs == 2 */
#define INIT_OUTER_LOOP_2   \
    INIT_OUTER_LOOP_1       \
    npy_intp s1 = *steps++;

/* Declare and initialise loop variables:
    s2 = stride along non-core dimension of third argument
When number of ufunc inputs + outputs == 3 */
#define INIT_OUTER_LOOP_3   \
    INIT_OUTER_LOOP_2       \
    npy_intp s2 = *steps++;

/* Declare and initialise loop variables:
    s3 = stride along non-core dimension of fourth argument
When number of ufunc inputs + outputs == 4 */
#define INIT_OUTER_LOOP_4   \
    INIT_OUTER_LOOP_3       \
    npy_intp s3 = *steps++;

/* Declare and initialise loop variables:
    s4 = stride along non-core dimension of fifth argument
When number of ufunc inputs + outputs == 5 */
#define INIT_OUTER_LOOP_5   \
    INIT_OUTER_LOOP_4       \
    npy_intp s4 = *steps++;

/* Declare and initialise loop variables:
    s5 = stride along non-core dimension of sixth argument
When number of ufunc inputs + outputs == 6 */
#define INIT_OUTER_LOOP_6   \
    INIT_OUTER_LOOP_5       \
    npy_intp s5 = *steps++;

/* Declare and initialise loop variables:
    s6 = stride along non-core dimension of sixth argument
When number of ufunc inputs + outputs == 7 */
#define INIT_OUTER_LOOP_7   \
    INIT_OUTER_LOOP_6       \
    npy_intp s6 = *steps++;

#define BEGIN_OUTER_LOOP        \
    for (N_ = 0; N_ < dN; N_++) {

/* Start of loop over non-core dimension */
#define BEGIN_OUTER_LOOP_2  BEGIN_OUTER_LOOP
#define BEGIN_OUTER_LOOP_3  BEGIN_OUTER_LOOP
#define BEGIN_OUTER_LOOP_4  BEGIN_OUTER_LOOP
#define BEGIN_OUTER_LOOP_5  BEGIN_OUTER_LOOP
#define BEGIN_OUTER_LOOP_6  BEGIN_OUTER_LOOP
#define BEGIN_OUTER_LOOP_7  BEGIN_OUTER_LOOP

/* End of loop over non-core dimension */
#define END_OUTER_LOOP  }

/* End of loop over non-core dimension
When number of ufunc inputs + outputs == 1 */
#define END_OUTER_LOOP_1  \
    args[0] += s0;        \
    END_OUTER_LOOP

/* End of loop over non-core dimension
When number of ufunc inputs + outputs == 2 */
#define END_OUTER_LOOP_2  \
    args[1] += s1;        \
    END_OUTER_LOOP_1

/* End of loop over non-core dimension
When number of ufunc inputs + outputs == 3 */
#define END_OUTER_LOOP_3  \
    args[2] += s2;        \
    END_OUTER_LOOP_2

/* End of loop over non-core dimension
When number of ufunc inputs + outputs == 4 */
#define END_OUTER_LOOP_4  \
    args[3] += s3;        \
    END_OUTER_LOOP_3

/* End of loop over non-core dimension
When number of ufunc inputs + outputs == 5 */
#define END_OUTER_LOOP_5  \
    args[4] += s4;        \
    END_OUTER_LOOP_4

/* End of loop over non-core dimension
When number of ufunc inputs + outputs == 6 */
#define END_OUTER_LOOP_6  \
    args[5] += s5;        \
    END_OUTER_LOOP_5

/* End of loop over non-core dimension
When number of ufunc inputs + outputs == 7 */
#define END_OUTER_LOOP_7  \
    args[6] += s6;        \
    END_OUTER_LOOP_6

/*
*****************************************************************************
**                      Error signaling functions                          **
*****************************************************************************
*/

/* Get current floating-point error status (to be reset if no new errors)
and set invalid value to false */
static NPY_INLINE int
get_fp_invalid_and_clear(void)
{
    int status;
    status = npy_clear_floatstatus_barrier((char*)&status);
    return !!(status & NPY_FPE_INVALID);
}

/* Reset floating-point error status if no new errors
else set invalid value to true */
static NPY_INLINE void
set_fp_invalid_or_clear(int error_occurred)
{
    if (error_occurred) {
        npy_set_floatstatus_invalid();
    }
    else {
        npy_clear_floatstatus_barrier((char*)&error_occurred);
    }
}

/*
*****************************************************************************
**                      Some handy constants                               **
*****************************************************************************
*/

static NPY_INLINE npy_intp
npy_int_min(npy_intp x, npy_intp y)
{
    return x < y ? x : y;
}

static NPY_INLINE npy_intp
npy_int_max(npy_intp x, npy_intp y)
{
    return x > y ? x : y;
}

/* comples number types with choice of interface */
#ifndef FORTRAN_TYPES
typedef union {
    npy_cfloat npy;
    float array[2];
} COMPLEX_t;

typedef union {
    npy_cdouble npy;
    double array[2];
} DOUBLECOMPLEX_t;
#endif

/* Constants for use as template parameters (numpy .src processing) */
static float s_one;
static float s_zero;
static float s_minus_one;
static float s_inf;
static float s_ninf;
static float s_nan;
static float s_eps;
static double d_one;
static double d_zero;
static double d_minus_one;
static double d_inf;
static double d_ninf;
static double d_nan;
static double d_eps;
static COMPLEX_t c_one;
static COMPLEX_t c_zero;
static COMPLEX_t c_minus_one;
static COMPLEX_t c_inf;
static COMPLEX_t c_ninf;
static COMPLEX_t c_nan;
static DOUBLECOMPLEX_t z_one;
static DOUBLECOMPLEX_t z_zero;
static DOUBLECOMPLEX_t z_minus_one;
static DOUBLECOMPLEX_t z_inf;
static DOUBLECOMPLEX_t z_ninf;
static DOUBLECOMPLEX_t z_nan;
static char char_C;
static char char_F;
static char char_L;
static char char_N;
static char char_R;
static char char_T;
static char char_U;


static void init_constants(void)
{
    /*
    this is needed as NPY_INFINITY and NPY_NAN macros
    can't be used as initializers. I prefer to just set
    all the constants the same way.
    */
    s_one  = 1.0f;
    s_zero = 0.0f;
    s_minus_one = -1.0f;
    s_inf = NPY_INFINITYF;
    s_ninf = -NPY_INFINITYF;
    s_nan = NPY_NANF;
    s_eps = npy_spacingf(s_one);

    d_one  = 1.0;
    d_zero = 0.0;
    d_minus_one = -1.0;
    d_inf = NPY_INFINITY;
    d_ninf = -NPY_INFINITY;
    d_nan = NPY_NAN;
    d_eps = npy_spacing(d_one);

    c_one.array[0]  = 1.0f;
    c_one.array[1]  = 0.0f;
    c_zero.array[0] = 0.0f;
    c_zero.array[1] = 0.0f;
    c_minus_one.array[0] = -1.0f;
    c_minus_one.array[1] = 0.0f;
    c_inf.array[0] = NPY_INFINITYF;
    c_inf.array[1] = 0.0f;
    c_ninf.array[0] = -NPY_INFINITYF;
    c_ninf.array[1] = 0.0f;
    c_nan.array[0] = NPY_NANF;
    c_nan.array[1] = NPY_NANF;

    z_one.array[0]  = 1.0;
    z_one.array[1]  = 0.0;
    z_zero.array[0] = 0.0;
    z_zero.array[1] = 0.0;
    z_minus_one.array[0] = -1.0;
    z_minus_one.array[1] = 0.0;
    z_inf.array[0] = NPY_INFINITY;
    z_inf.array[1] = 0.0;
    z_ninf.array[0] = -NPY_INFINITY;
    z_ninf.array[1] = 0.0;
    z_nan.array[0] = NPY_NAN;
    z_nan.array[1] = NPY_NAN;

    char_C = 'C';
    char_F = 'F';
    char_L = 'L';
    char_N = 'N';
    char_R = 'R';
    char_T = 'T';
    char_U = 'U';
}

/*
*****************************************************************************
**                             Ufunc definition                            **
*****************************************************************************
*/
/* For the 'data' argument for ufunc creation */
static void *null_data_array[] = {
    (void *)NULL, (void *)NULL, (void *)NULL, (void *)NULL, (void *)NULL,
    (void *)NULL, (void *)NULL, (void *)NULL, (void *)NULL, (void *)NULL,
    (void *)NULL, (void *)NULL, (void *)NULL, (void *)NULL, (void *)NULL,
    (void *)NULL, (void *)NULL, (void *)NULL, (void *)NULL, (void *)NULL,
};

/* For the 'types' argument for ufunc creation */
static char ufn_types_2_2[] = {NPY_FLOAT, NPY_FLOAT,
                                NPY_DOUBLE, NPY_DOUBLE};
static char ufn_types_2_3[] = {NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
                                NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};
static char ufn_types_2_4[] = {NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
                            NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};
static char ufn_types_2_5[] = {
                NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
                NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};
static char ufn_types_2_6[] = {
        NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
        NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};

static char ufn_types_3_3[] = { NPY_LONG, NPY_LONG, NPY_LONG,
                                NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
                                NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE };

static char ufn_types_4_2[] = {NPY_FLOAT, NPY_FLOAT,
                                NPY_DOUBLE, NPY_DOUBLE,
                                NPY_CFLOAT, NPY_CFLOAT,
                                NPY_CDOUBLE, NPY_CDOUBLE };
static char ufn_types_4_3[] = {NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
                                NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
                                NPY_CFLOAT, NPY_CFLOAT, NPY_CFLOAT,
                                NPY_CDOUBLE, NPY_CDOUBLE, NPY_CDOUBLE };
static char ufn_types_4_4[] = {NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
                            NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
                            NPY_CFLOAT, NPY_CFLOAT, NPY_CFLOAT, NPY_CFLOAT,
                            NPY_CDOUBLE, NPY_CDOUBLE, NPY_CDOUBLE, NPY_CDOUBLE};
static char ufn_types_4_5[] = {
            NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
            NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
            NPY_CFLOAT, NPY_CFLOAT, NPY_CFLOAT, NPY_CFLOAT, NPY_CFLOAT,
            NPY_CDOUBLE, NPY_CDOUBLE, NPY_CDOUBLE, NPY_CDOUBLE, NPY_CDOUBLE};

static char ufn_types_5_3[] = { NPY_INT, NPY_INT, NPY_INT,
                                NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
                                NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
                                NPY_CFLOAT, NPY_CFLOAT, NPY_CFLOAT,
                                NPY_CDOUBLE, NPY_CDOUBLE, NPY_CDOUBLE };

/* For the 'func' argument for ufunc creation */
#define FUNC_ARRAY_NAME(NAME) NAME ## _funcs

#define GUFUNC_FUNC_ARRAY_REAL(NAME)                    \
    static PyUFuncGenericFunction                       \
    FUNC_ARRAY_NAME(NAME)[] = {                         \
        FLOAT_ ## NAME,                                 \
        DOUBLE_ ## NAME                                 \
    }

#define GUFUNC_FUNC_ARRAY_REAL_INT(NAME)                \
    static PyUFuncGenericFunction                       \
    FUNC_ARRAY_NAME(NAME)[] = {                         \
        INT_ ## NAME,                                   \
        FLOAT_ ## NAME,                                 \
        DOUBLE_ ## NAME                                 \
    }

#define GUFUNC_FUNC_ARRAY_REAL_COMPLEX(NAME)            \
    static PyUFuncGenericFunction                       \
    FUNC_ARRAY_NAME(NAME)[] = {                         \
        FLOAT_ ## NAME,                                 \
        DOUBLE_ ## NAME,                                \
        CFLOAT_ ## NAME,                                \
        CDOUBLE_ ## NAME                                \
    }

#define GUFUNC_FUNC_ARRAY_REAL_COMPLEX_INT(NAME)        \
    static PyUFuncGenericFunction                       \
    FUNC_ARRAY_NAME(NAME)[] = {                         \
        INT_ ## NAME,                                   \
        FLOAT_ ## NAME,                                 \
        DOUBLE_ ## NAME,                                \
        CFLOAT_ ## NAME,                                \
        CDOUBLE_ ## NAME                                \
    }

/* for collecting arguments for ufunc creation */
typedef struct gufunc_descriptor_struct {
    char *name;
    char *signature;
    char *doc;
    int ntypes;
    int nin;
    int nout;
    PyUFuncGenericFunction *funcs;
    char *types;
} GUFUNC_DESCRIPTOR_t;

/* Creating ufuncs and adding them to the module dict .
To be called in PyInit__gufuncs_<module name>.
    module: module object from PyModule_Create
    guf_descriptors[]: array of GUFUNC_DESCRIPTOR_t above
    gufunc_count: length of guf_descriptors[]
    version_string: char array containing version number.
*/
static int
addUfuncs(PyObject *module, const GUFUNC_DESCRIPTOR_t guf_descriptors[],
            int gufunc_num, const char *version_string)
{
    PyObject *dictionary;
    PyObject *version;

    /* module's dictionary (__dir__), borrowed reference */
    dictionary = PyModule_GetDict(module);

    /* new reference */
    version = PyString_FromString(version_string);
    PyDict_SetItemString(dictionary, "__version__", version);
    Py_DECREF(version);

    PyObject *f;
    int i;
    /* int gufunc_num = sizeof(*guf_descriptors) / sizeof(guf_descriptors[0]); */
    for (i = 0; i < gufunc_num; i++) {
        /* current gufunc descriptor */
        const GUFUNC_DESCRIPTOR_t* d = &guf_descriptors[i];
        /* create gufunc object (new reference) */
        f = PyUFunc_FromFuncAndDataAndSignature(d->funcs, null_data_array,
                            d->types, d->ntypes, d->nin, d->nout,
                            PyUFunc_None, d->name, d->doc, 0, d->signature);
        if (f == NULL) {
            return -1;
        }
        /* add gufunc object to module's dictionary */
        PyDict_SetItemString(dictionary, d->name, f);
        Py_DECREF(f);
    }
    return 0;
}

#endif
