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

static const char* baum_welch_version = "0.1.0";

/*
*****************************************************************************
**                   Doc string for Python functions                       **
*****************************************************************************
*/

/* alpha_beta_signature = "(r,p,m,m),(r,m),(tm),(t)->(t,m),(t,m),(t)"; */

PyDoc_STRVAR(alpha_beta__doc__,
/* "alpha_beta(M, S, P, R: ndarray) -> (A, B, E: ndarray)\n\n" */
"Calculate BW forward/backward variables\n"
"\nParameters\n-----------\n"
"updaters : la.lnarray, (R,P,M,M)\n"
"    Plasticity matrices multiplied by readout indicators of 'to' state,\n"
"    `Prob(i(t+1)=j, w(t+1)=r|i(t)=i, mu(t)=p)`.\n"
"initial : la.lnarray, (R,M)\n"
"    Initial state distribution multiplied by readout indicators of state,\n"
"    `Prob(w(0)=r, i(0)=i)`.\n"
"plast_type : ArrayLike, (T-1,E), int[0:P]\n"
"    id of plasticity type after each time-step.\n"
"readouts : ArrayLike, (T,E), int[0:R]\n"
"    id of readout from synapse at each time-step.\n"
"\nReturns\n-------\n"
"alpha : la.lnarray, (T,M)\n"
"    Normalised Baum-Welch forward variable.\n"
"beta : la.lnarray, (T,M)\n"
"    Scaled Baum-Welch backward variable.\n"
"eta : la.lnarray, (T,)\n"
"    Norm of Baum-Welch forward variable.\n"
);


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

/* set vector x to constant */
extern void
FNAME(dlaset)(char trans[], int *m, int *n,
            double *alpha, double *beta,
            double a[], int *lda);

/* x -> x / a */
extern void
FNAME(drscl)(int *n, double *alpha, double x[], int *incx);

/* x -> ||x||_1 */
extern double
FNAME(dasum)(int *n, double x[], int *incx);

/* y -> a A*x + b y */
extern void
FNAME(dgemv)(char trans[], int *m, int *n,
    double *alpha, double a[], int *lda, double x[], int *incx,
    double *beta, double y[], int *incy);

/* Need:
    dasum, (dscal), dgemv, drscl
*/

/*
******************************************************************************
**                             ALPHA_BETA                                   **
******************************************************************************
*/
/* to hold arguments for BLAS _gemm */
typedef struct gemv_params_struct
{
  void *A; /* A is scalar of base type */
  void *B; /* B is scalar of base type */
  void *E; /* E is scalar of base type */
  void *X; /* X is (M,) of base type */
  void *Y; /* Y is (M) of base type */
  void *Z;  /* Z is (M,M) of base type */
  void *ZZ; /* ZZ is (R,P,M,M) of base type */
  void *S;  /* S is (M) of base type */
  void *SS; /* SS is (R,M) of base type */

  /* strides for choosing */
  size_t IP;
  size_t IR;
  size_t ISR;
  /* for fortran functions */
  fortran_int M;
  fortran_int P;
  fortran_int R;
  fortran_int ISM;
  fortran_int IX;
  fortran_int IY;
  fortran_int LDZ;
  char TRANSZ;
} GEMV_PARAMS_t;

/***************************************************
* Calling BLAS functions drscl, dasum, dgemv       *
****************************************************/

/* Copy vector: X -> S */
static NPY_INLINE void
call_dcopy_x_i(GEMV_PARAMS_t *params)
{
    BLAS(dcopy)(&params->M, params->S, &params->ISM, params->X, &params->IX);
}

/* Set vector: X -> [1,1,...] */
static NPY_INLINE void
call_dlaset_x(GEMV_PARAMS_t *params)
{
    LAPACK(dlaset)(&char_N, &params->M, &params->IX,
                    params->A, params->A, params->X, &params->M);
}

/* Scale vector: X -> X/E */
static NPY_INLINE void
call_drscl_x_e(GEMV_PARAMS_t *params)
{
    if (*(npy_double *)params->E > 0.0)
    {
        BLAS(drscl)(&params->M, params->E, params->X, &params->IX);
    }
}

/* L1 norm: E -> ||X||_1 */
static NPY_INLINE void
call_dasum_e_x(GEMV_PARAMS_t *params)
{
    npy_double norm = BLAS(dasum)(&params->M, params->X, &params->IX);
    *(npy_double *)(params->E) = norm;
}

/* Matrix vector product: X -> Z.X */
static NPY_INLINE void
call_dgemv_x_z_x(GEMV_PARAMS_t *params)
{
    BLAS(dgemv)(&params->TRANSZ, &params->M, &params->M,
        params->A,  params->Z, &params->LDZ, params->X, &params->IX,
        params->B, params->Y, &params->IY);
    npy_uint8 *a;
    a = params->X;
    params->X = params->Y;
    params->Y = a;
}

/********************************************************************
* Initialize the parameters to use in for the lapack function _gemv *
* Handles buffer allocation                                         *
*********************************************************************/
/* initialise parameters for BLAS
    M: matrix dimension
    P: number of matrices
    R: number of readouts */
static NPY_INLINE int
init_DOUBLE_alpha_beta(GEMV_PARAMS_t *params,
                npy_intp M_in, npy_intp P_in, npy_intp R_in)
{
    npy_uint8 *mem_buff = NULL;
    npy_uint8 *a, *b, *c, *d;
    fortran_int M, P, R, ldz;
    size_t safe_M, safe_P, safe_R;

    M = (fortran_int)M_in;
    P = (fortran_int)P_in;
    R = (fortran_int)R_in;
    safe_M = M_in;
    safe_P = P_in;
    safe_R = R_in;
    ldz = fortran_int_max(M, 1);

    mem_buff = malloc(safe_M * sizeof(fortran_doublereal)
                      + safe_M * sizeof(fortran_doublereal)
                      + safe_R * safe_M * sizeof(fortran_doublereal)
                      + safe_R * safe_P * safe_M * safe_M * sizeof(fortran_doublereal)
                      );
    if (!mem_buff) {
        goto error;
    }
    a = mem_buff;  /* X at start of buffer */
    b = a + safe_M * sizeof(fortran_doublereal);  /* Y after space for X */
    c = b + safe_M * sizeof(fortran_doublereal);  /* S after space for Y */
    d = c + safe_R * safe_M * sizeof(fortran_doublereal);  /* Z after space for S */

    params->TRANSZ = 'N';
    params->A = &d_one;  /* Y -> A*Z.X + B*Y, need A=1, B=0 */
    params->B = &d_zero;
    params->X = a;
    params->Y = b;
    params->S = c;
    params->SS = c;
    params->Z = d;
    params->ZZ = d;
    params->IP = safe_M * safe_M * sizeof(fortran_doublereal);
    params->IR = safe_P * (params->IP);
    params->ISR = safe_M * sizeof(fortran_doublereal);
    params->M = M;
    params->P = P;
    params->R = R;
    params->ISM = 1;
    params->IX = 1;
    params->IY = 1;
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
release_DOUBLE_alpha_beta(GEMV_PARAMS_t *params)
{
    /* memory block base is in X */
    free(params->X);
    memset(params, 0, sizeof(*params));
}

/************************************
* Move parameters                   *
*************************************/

static NPY_INLINE int
choose_updater(GEMV_PARAMS_t *params, const char *plasttype, const char *readout)
{
    npy_int P_in = *(npy_int *)plasttype;
    npy_int R_in = *(npy_int *)readout;
    if (P_in < params->P && R_in < params->R)
    {
        npy_uint8 *a;
        a = params->ZZ;
        params->Z = a + (size_t)P_in * params->IP + (size_t)R_in * params->IR;
        return 0;
    } else {
        // printf("[P_in: %d, P_max: %d. %d\n", P_in, params->P, P_in < params->P);
        // printf(" R_in: %d, R_max: %d. %d]\n", R_in, params->R, R_in < params->R);
        return 1;
    }
}

static NPY_INLINE int
choose_init(GEMV_PARAMS_t *params, const char *readout)
{
    npy_int R_in = *(npy_int *)readout;
    if (R_in < params->R)
    {
        npy_uint8 *a;
        a = params->SS;
        params->S = a + (size_t)R_in * params->ISR;
        return 0;
    } else {
        // printf("[R_in: %d, R_max: %d. %d]\n", R_in, params->R, R_in < params->R);
        return 1;
    }
}

/*****************************
* BLAS Inner GUfunc loop     *
******************************/

/* copy initial, updaters into struct for fortran */
static void
linearize_updaters(const char **args, npy_intp *dimensions, npy_intp *steps,
                   GEMV_PARAMS_t *params, const LINEARIZE_DATA_t *lin_data)
{
    LINEARIZE_DATA_t *u_in = lin_data++;
    LINEARIZE_DATA_t *s_in = lin_data++;

    const char *ip_updater = args[0];  //  in-ptr: updater readouts
    const char *ip_updatep = args[0];  //  in-ptr: updater plast types
    const char *ip_initial = args[1];  //  in-ptr: init readouts
    npy_intp cnt_r, cnt_p;

    /* Copy updaters, initials */
    for (cnt_r = 0; cnt_r < dimensions[0]; cnt_r++)
    {
        choose_init(params, cnt_r);
        linearize_DOUBLE_vec(params->S, ip_initial, s_in);
        ip_updatep = ip_updater;

        for (cnt_p = 0; cnt_p < dimensions[1]; cnt_p++)
        {
            choose_updater(params, cnt_p, cnt_r);
            linearize_DOUBLE_matrix(params->Z, ip_updatep, u_in);
            ip_updatep += steps[0];
        }
        ip_updater += steps[1];
        ip_initial += steps[2];
    }
}

/* normalise, copy out alpha, eta */
static void
norm_alpha(char *alpha, char *eta, GEMV_PARAMS_t *params, const LINEARIZE_DATA_t *a_out)
{
    /* normalise */
    call_dasum_e_x(params);
    call_drscl_x_e(params);
    /* copy alpha out */
    delinearize_DOUBLE_vec(alpha, params->X, a_out);
    /* eta(0) is 1/eta from here */
    *(npy_double *)eta = *(npy_double *)(params->E);
}

/* scale, copy out beta, invert eta */
static void
scale_beta(char *beta, char *eta, GEMV_PARAMS_t *params, const LINEARIZE_DATA_t *b_out)
{
    /* scale by eta */
    *(npy_double *)(params->E) = *(npy_double *)eta;
    call_drscl_x_e(params);
    /* copy beta out */
    delinearize_DOUBLE_vec(beta, params->X, b_out);
    /* eta(T) back to eta from here */
    if (*(npy_double *)params->E > 0.0)
    {
        *(npy_double *)eta = 1. / *(npy_double *)eta;
    }
}

/* alpha_beta_signature = "(r,p,m,m),(r,m),(tm),(t)->(t,m),(t,m),(t)"; */

static void
DOUBLE_alpha_beta(char **args, npy_intp *dimensions, npy_intp *steps,
              void *NPY_UNUSED(func))
{
INIT_OUTER_LOOP_7

    npy_intp len_r = dimensions[0];  /* number of readouts */
    npy_intp len_p = dimensions[1];  /* number of matrices */
    npy_intp len_m = dimensions[2];  /* number of states */
    npy_intp len_tm = dimensions[3]; /* number of time-points - 1 */
    npy_intp len_t = dimensions[4];  /* number of time-points */
    npy_intp stride_u_r = *steps++;  /* 1st arg */
    npy_intp stride_u_p = *steps++;
    npy_intp stride_u_mr = *steps++;
    npy_intp stride_u_mc = *steps++;
    npy_intp stride_s_r = *steps++;  /* 2nd arg */
    npy_intp stride_s_m = *steps++;
    npy_intp stride_p_t = *steps++;  /* 3rd arg */
    npy_intp stride_r_t = *steps++;  /* 4th arg */
    npy_intp stride_a_t = *steps++;  /* 1st output */
    npy_intp stride_a_m = *steps++;
    npy_intp stride_b_t = *steps++;  /* 2nd output */
    npy_intp stride_b_m = *steps++;
    npy_intp stride_e_t = *steps++;  /* 3rd output */
    int error_occurred = get_fp_invalid_and_clear();
    GEMV_PARAMS_t params;
    LINEARIZE_DATA_t u_in, s_in, a_out, b_out;
    npy_intp cnt_t;

    if (len_t != len_tm + 1)
    {
        error_occurred = 1;
    }


   /* allocate buffer */
    if (error_occurred == 0 &&
        init_DOUBLE_alpha_beta(&params, len_m, len_p, len_r))
    {
        /* initialise size parameters */
        init_linearize_data(&u_in, len_m, len_m, stride_u_mc, stride_u_mr);
        init_linearize_vdata(&s_in, len_m, stride_s_m);
        init_linearize_vdata(&a_out, len_m, stride_a_m);
        init_linearize_vdata(&b_out, len_m, stride_b_m);

        /* Copy updaters, initials, assume no broadcasting */
        npy_intp rsteps[] = {stride_u_p, stride_u_r, stride_s_r};
        LINEARIZE_DATA_t lin_data[] = {u_in, s_in};
        linearize_updaters(args, dimensions, rsteps, &params, lin_data);

        BEGIN_OUTER_LOOP

            const char *ip_plast = args[2];  //  in-ptr: plast_type
            const char *ip_reado = args[3];  //  in-ptr: readouts
            char *op_alpha = args[4];  //  out-ptr: alpha
            char *op_beta = args[5];   //  out-ptr: beta
            char *op_eta = args[6];    //  out-ptr: eta
            params.TRANSZ = 'T';

            /* first alpha */
            error_occurred = choose_init(&params, ip_reado);
            if (error_occurred)
                { break; }

            /* set first alpha to init */
            call_dcopy_x_i(&params);
            /* normalise, copy out alpha, eta */
            norm_alpha(op_alpha, op_eta, &params, &a_out);
            /* eta(0) is 1/eta from here */

            for (cnt_t = 0; cnt_t < len_tm; cnt_t++)
            {
                ip_reado += stride_r_t;
                op_alpha += stride_a_t;
                op_eta += stride_e_t;

                error_occurred = choose_updater(&params, ip_plast, ip_reado);
                if (error_occurred)
                    { break; }

                /* forward step */
                call_dgemv_x_z_x(&params);
                /* normalise, copy out alpha, eta */
                norm_alpha(op_alpha, op_eta, &params, &a_out);
                /* eta(t) is 1/eta from here */

                ip_plast += stride_p_t;
            }
            if (error_occurred)
                { break; }

            params.TRANSZ = 'N';
            ip_plast -= stride_p_t;
            op_beta += len_tm * stride_b_t;
            /* ip_plast, ip_reado, op_beta, op_eta all point to last element */

            /* set last beta to one */
            call_dlaset_x(&params);
            /* scale, copy out beta, invert eta */
            scale_beta(op_beta, op_eta, &params, &b_out);
            /* eta(T) back to eta from here */

            for (cnt_t = 0; cnt_t < len_tm; cnt_t++)
            {
                op_beta -= stride_b_t;
                op_eta -= stride_e_t;

                choose_updater(&params, ip_plast, ip_reado);
                /* backward step */
                call_dgemv_x_z_x(&params);
                /* scale, copy out beta, invert eta */
                scale_beta(op_beta, op_eta, &params, &b_out);
                /* eta(T) back to eta from here */

                ip_plast -= stride_p_t;
                ip_reado -= stride_r_t;
            }
            // ip_plast += stride_p_t;

        END_OUTER_LOOP_7
        /* deallocate buffer */
        release_DOUBLE_alpha_beta(&params);
    }
    set_fp_invalid_or_clear(error_occurred);
}

/*
*****************************************************************************
**                             Ufunc definition                            **
*****************************************************************************
*/

/* types argument for creating 'norm' ufunc */
static char ufn_types_1_7[] = {NPY_DOUBLE, NPY_DOUBLE, NPY_INT, NPY_INT,
                               NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};

/* array of functions for each ufunc loop */
static PyUFuncGenericFunction
FUNC_ARRAY_NAME(alpha_beta)[] = {
    DOUBLE_alpha_beta
};

/* info for creating ufunc object */
GUFUNC_DESCRIPTOR_t gufunc_descriptors[] = {
    {"alpha_beta", "(r,p,m,m),(r,m),(tm),(t)->(t,m),(t,m),(t)", alpha_beta__doc__,
        1, 4, 3, FUNC_ARRAY_NAME(alpha_beta), ufn_types_1_7}
};

/*
*****************************************************************************
**               Module initialization stuff                               **
*****************************************************************************
*/

/* Methods to add to module (none, we add ufuncs after creating them) */
static PyMethodDef Baum_Welch_Methods[] = {
    /* Sentinel */
    {NULL, NULL, 0, NULL}
};

/* arguments for module creation */
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_baum_welch",
        NULL,
        -1,
        Baum_Welch_Methods,
        NULL,
        NULL,
        NULL,
        NULL
};

/* create module */
PyObject *PyInit__baum_welch(void)
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
    int failure = addUfuncs(m, gufunc_descriptors, 1, baum_welch_version);

    if (PyErr_Occurred() || failure) {
        PyErr_SetString(PyExc_RuntimeError,
                        "cannot load _baum_welch module.");
        return NULL;
    }

    return m;
}

