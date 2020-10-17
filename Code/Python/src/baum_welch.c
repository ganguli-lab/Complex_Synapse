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

/* m_init_signature = "(r,p,m,m),(tm),(t),(t,m),(t,m),(t)->(p,m,m),(m)"; */

PyDoc_STRVAR(plast_init__doc__,
/* "m_init(M, P, R, A, B, E: ndarray) -> (UM, US: ndarray)\n\n" */
"One Baum-Welch/Rabiner-Juang update of the model\n"
"\nParameters\n-----------\n"
"updaters : la.lnarray, (R,P,M,M) float[0:1]\n"
"    Plasticity matrices multiplied by readout probability given 'to' state.\n"
"plast_type : la.lnarray, (T-1,E), int[0:P]\n"
"    id of plasticity type after each time-step\n"
"alpha : la.lnarray, (T,M) float[0:1]\n"
"    Normalised Baum-Welch forward variable\n"
"beta : la.lnarray, (T,M) float\n"
"    Scaled Baum-Welch backward variable\n"
"eta : la.lnarray, (T,) float[1:]\n"
"    Norm of Baum-Welch forward variable\n"
"\nReturns\n-------\n"
"plast : array_like, (P,M,M), float[0:1]\n"
"    new estimate of transition probability matrix.\n"
"initial : array_like, (M,) float[0:1]\n"
"    new estimate of distribution of initial state,\n"
"    not assumed to be the steady-state distribution.\n"
);

PyDoc_STRVAR(plast_steady__doc__,
/* "m_init(M, P, R, A, B, E: ndarray) -> (UM, US: ndarray)\n\n" */
"One Baum-Welch/Rabiner-Juang update of the model\n"
"\nParameters\n-----------\n"
"updaters : la.lnarray, (R,P,M,M) float[0:1]\n"
"    Plasticity matrices multiplied by readout probability given 'to' state.\n"
"plast_type : la.lnarray, (T-1,E), int[0:P]\n"
"    id of plasticity type after each time-step\n"
"alpha : la.lnarray, (T,M) float[0:1]\n"
"    Normalised Baum-Welch forward variable\n"
"beta : la.lnarray, (T,M) float\n"
"    Scaled Baum-Welch backward variable\n"
"eta : la.lnarray, (T,) float[1:]\n"
"    Norm of Baum-Welch forward variable\n"
"\nReturns\n-------\n"
"plast : array_like, (P,M,M), float[0:1]\n"
"    new estimate of transition probability matrix.\n"
"initial : array_like, (M,) float[0:1]\n"
"    new estimate of distribution of initial state,\n"
"    assumed to be the steady-state distribution.\n"
);

/*
*****************************************************************************
**                   BLAS/Lapack declarations                              **
*****************************************************************************
*/

/* copy vector x into y */
extern void
FNAME(dcopy)(int *n, double *x, int *incx, double *y, int *incy);

/* set vector x to constant */
extern void
FNAME(dlaset)(char uplo[], int *m, int *n,
            double *alpha, double *beta,
            double a[], int *lda);

/* y -> a x + y */
// extern void
// FNAME(daxpy)(int *n, double *a, double x[], int *incx,
//             double y[], int *incy);

/* x -> x / a */
extern void
FNAME(drscl)(int *n, double *alpha, double x[], int *incx);

/* _ -> ||x||_1 */
extern double
FNAME(dasum)(int *n, double x[], int *incx);

// /* _ -> x^T y */
// extern double
// FNAME(ddot)(int *n, double x[], int *incx, double y[], int *incy);

/* A -> A + a x y^T */
extern void
FNAME(dger)(int *m, int *n,
            double *alpha, double x[], int *incx,
            double y[], int *incy, double a[], int *lda);

/* y -> a A*x + b y */
extern void
FNAME(dgemv)(char trans[], int *m, int *n,
            double *alpha, double a[], int *lda, double x[], int *incx,
            double *beta, double y[], int *incy);

/* y -> a A*x + b y */
extern void
FNAME(dgbmv)(char trans[], int *m, int *n, int *kl, int *ku,
            double *alpha, double a[], int *lda, double x[], int *incx,
            double *beta, double y[], int *incy);

/* Need:
    dlaset, dasum, (dscal), dgemv, drscl, dger, dgbmv, daxpy
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
call_dcopy_x_s(GEMV_PARAMS_t *params)
{
    BLAS(dcopy)(&params->M, params->S, &params->ISM, params->X, &params->IX);
}

/* Set vector: X -> [1,1,...] */
static NPY_INLINE void
call_dlaset_x_one(GEMV_PARAMS_t *params)
{
    LAPACK(dlaset)(&char_N, &params->M, &params->ISM,  /* N=1 */
                    params->A, params->A, params->X, &params->M);  /* alpha=beta=1 */
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
    *(fortran_doublereal *)(params->E) = norm;
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

    size_t Espace, Xspace, Yspace, ZZspace, SSspace, Zspace, Sspace;

    Espace = sizeof(fortran_doublereal);
    Xspace = safe_M * Espace;
    Yspace = Xspace;
    Sspace = Xspace;
    Zspace = safe_M * Xspace;
    ZZspace = safe_R * safe_P * Zspace;
    SSspace = safe_R * Sspace;

    mem_buff = malloc(Espace + Xspace + Yspace + ZZspace + SSspace);
    if (!mem_buff) {
        goto error;
    }
    /* E at start of buffer */
    a = mem_buff + Espace;  /* X after space for E */
    b = a + Xspace;  /* Y after space for X */
    c = b + Yspace;  /* ZZ after space for Y */
    d = c + ZZspace;  /* SS after space for ZZ */

    params->TRANSZ = 'N';
    params->A = &d_one;  /* Y -> A*Z.X + B*Y, need A=1, B=0 */
    params->B = &d_zero;
    params->E = mem_buff;
    params->X = a;
    params->Y = b;
    params->ZZ = c;
    params->SS = d;
    params->Z = c;
    params->S = d;
    params->IP = Zspace;
    params->IR = safe_P * Zspace;
    params->ISR = Sspace;
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
    /* memory block base is in E */
    free(params->E);
    memset(params, 0, sizeof(*params));
}

/************************************
* Move parameters                   *
*************************************/

static NPY_INLINE int
choose_updater_in_n(GEMV_PARAMS_t *params, npy_int P_in, npy_int R_in)
{
    if (P_in < (npy_int)(params->P) && R_in < (npy_int)(params->R))
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
choose_updater_in(GEMV_PARAMS_t *params, const char *plasttype, const char *readout)
{
    npy_int P_in = *(npy_int *)plasttype;
    npy_int R_in = *(npy_int *)readout;
    return choose_updater_in_n(params, P_in, R_in);
}

static NPY_INLINE int
choose_init_n(GEMV_PARAMS_t *params, npy_int R_in)
{
    if (R_in < (npy_int)(params->R))
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

static NPY_INLINE int
choose_init(GEMV_PARAMS_t *params, const char *readout)
{
    npy_int R_in = *(npy_int *)readout;
    return choose_init_n(params, R_in);
}

/*****************************
* BLAS Inner GUfunc loop     *
******************************/

/* copy initial, updaters into struct for fortran */
static void
linearize_updater_init(const char **args, npy_intp *dimensions, npy_intp *steps,
                   GEMV_PARAMS_t *params, const LINEARIZE_DATA_t *lin_data)
{
    const LINEARIZE_DATA_t *u_in = lin_data++;
    const LINEARIZE_DATA_t *s_in = lin_data++;

    const char *ip_updater = args[0];  //  in-ptr: updater, readouts
    const char *ip_updatep = args[0];  //  in-ptr: updater, plast types
    const char *ip_initial = args[1];  //  in-ptr: init, readouts
    npy_intp cnt_r, cnt_p;

    /* Copy updaters, initials */
    for (cnt_r = 0; cnt_r < dimensions[0]; cnt_r++)
    {
        choose_init_n(params, cnt_r);
        linearize_DOUBLE_vec(params->S, ip_initial, s_in);
        ip_updatep = ip_updater;

        for (cnt_p = 0; cnt_p < dimensions[1]; cnt_p++)
        {
            choose_updater_in_n(params, cnt_p, cnt_r);
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
    /* copy beta out */
    delinearize_DOUBLE_vec(beta, params->X, b_out);
    /* scale beta(t) by eta(t) (for next backward step) */
    *(npy_double *)(params->E) = *(npy_double *)eta;
    call_drscl_x_e(params);
    /* eta(t) back to eta from here */
    if (*(npy_double *)eta > 0.0)
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
    // printf("allocate memory\n");
    if (error_occurred == 0 &&
        init_DOUBLE_alpha_beta(&params, len_m, len_p, len_r))
    {
        /* initialise size parameters */
        init_linearize_data(&u_in, len_m, len_m, stride_u_mc, stride_u_mr);
        init_linearize_vdata(&s_in, len_m, stride_s_m);
        init_linearize_vdata(&a_out, len_m, stride_a_m);
        init_linearize_vdata(&b_out, len_m, stride_b_m);
        // printf("linearize updaters\n");

        /* Copy updaters, initials, assume no broadcasting */
        npy_intp rsteps[] = {stride_u_p, stride_u_r, stride_s_r};
        LINEARIZE_DATA_t lin_data[] = {u_in, s_in};
        linearize_updater_init(args, dimensions, rsteps, &params, lin_data);
        // printf("updaters linearized\n");

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
            call_dcopy_x_s(&params);
            /* normalise, copy out alpha, eta */
            norm_alpha(op_alpha, op_eta, &params, &a_out);
            /* eta(0) is 1/eta from here */
            // printf("begin alpha loop\n");

            for (cnt_t = 0; cnt_t < len_tm; cnt_t++)
            {
                ip_reado += stride_r_t;
                op_alpha += stride_a_t;
                op_eta += stride_e_t;

                error_occurred = choose_updater_in(&params, ip_plast, ip_reado);
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
            // printf("end alpha loop\n");

            params.TRANSZ = 'N';
            ip_plast -= stride_p_t;
            op_beta += len_tm * stride_b_t;
            /* ip_plast, ip_reado, op_beta, op_eta all point to last element */

            /* set last beta to one */
            call_dlaset_x_one(&params);
            /* copy out, scale beta, invert eta */
            scale_beta(op_beta, op_eta, &params, &b_out);
            /* eta(T) back to eta from here */
            // printf("begin beta loop\n");

            for (cnt_t = 0; cnt_t < len_tm; cnt_t++)
            {
                op_beta -= stride_b_t;
                op_eta -= stride_e_t;

                choose_updater_in(&params, ip_plast, ip_reado);
                /* backward step */
                call_dgemv_x_z_x(&params);
                /* copy out, scale beta, invert eta */
                scale_beta(op_beta, op_eta, &params, &b_out);
                /* eta(T) back to eta from here */

                ip_plast -= stride_p_t;
                ip_reado -= stride_r_t;
            }
            // ip_plast += stride_p_t;
            // printf("end beta loop\n");

        END_OUTER_LOOP_7
        /* deallocate buffer */
        release_DOUBLE_alpha_beta(&params);
        // printf("memory freed\n");
    }
    set_fp_invalid_or_clear(error_occurred);
}

/*
******************************************************************************
**                          UPDATE_M_INIT                                   **
******************************************************************************
*/
/* to hold arguments for BLAS _gemm */
typedef struct gbmv_params_struct
{
  void *A; /* A is scalar of double */
  void *B; /* B is scalar of double */
  void *E; /* E is scalar of double */
  void *X; /* X is (M) of double */
  void *Y; /* Y is (M) of double */
  void *ZZ;  /* ZZ is (R,P,M,M) of double, (M,M,P,R) in fortran */
  void *Z; /* Z is (M,M) of double */
  void *ZTMP; /* ZTMP is (R,P,M,M) of double, (M,M,P,R) in fortran */
  void *UZZ; /* UZZ is (P,M,M) of double, (M,M,P) in fortran */
  void *UZ;  /* UZ is (M,M) of double */
  void *US;  /* US is (M) of double */

  /* strides for choosing */
  size_t IP;
  size_t IR;
  /* for fortran functions */
  fortran_int K;
  fortran_int M;
  fortran_int P;
  fortran_int R;
  fortran_int RPM;
  fortran_int PM;
  fortran_int PMM;
  fortran_int ISM;
  fortran_int IX;
  fortran_int IY;
  fortran_int LDZ;
  char TRANSZ;
} GBMV_PARAMS_t;

/***************************************************
* Calling BLAS functions drscl, dasum, dgemv       *
****************************************************/

/* Clear US, ZTMP, UZZ -> 0 */
static NPY_INLINE void
call_dlaset_u_zero(GBMV_PARAMS_t *params)
{
    LAPACK(dlaset)(&char_N, &params->M, &params->ISM,  /* N=1 */
                    params->B, params->B, params->US, &params->LDZ);  /* alpha=beta=0 */
    LAPACK(dlaset)(&char_N, &params->M, &params->RPM,  /* N=RPM */
                    params->B, params->B, params->ZTMP, &params->LDZ);  /* alpha=beta=0 */
    LAPACK(dlaset)(&char_N, &params->M, &params->PM,  /* N=PM */
                    params->B, params->B, params->UZZ, &params->LDZ);  /* alpha=beta=0 */
}

/* Ranx one update: UZ -> E X.Y^T */
static NPY_INLINE void
call_dger_uz_x_y(GBMV_PARAMS_t *params)
{
    LAPACK(dger)(&params->M, &params->M, params->E,
            params->X, &params->IX, params->Y, &params->IY,
            params->UZ, &params->LDZ);
}

/* Multiply and add: US -> US + X o Y

To achieve this, X->diag(X)
Then it is US -> US + diag(X).Y
*/
call_dgbmv_us_x_y(GBMV_PARAMS_t *params)
{
    LAPACK(dgbmv)(&char_N, &params->M, &params->M, &params->K, &params->K,
            params->A, params->X, &params->IX, params->Y, &params->IY,
            params->A, params->US, &params->ISM);
}

/* Multiply and add: UZZ -> UZZ + Z o UZ

To achieve this, vectorise UZZ,Z,UZ, and Z->diag(Z)
Then it is UZZ -> UZZ + diag(Z).UZ
*/
call_dgbmv_uuz_z_uz(GBMV_PARAMS_t *params)
{
    LAPACK(dgbmv)(&char_N, &params->PMM, &params->PMM, &params->K, &params->K,
            params->A, params->Z, &params->ISM, params->UZ, &params->IX,
            params->A, params->UZZ, &params->IY);
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
init_DOUBLE_m_init(GBMV_PARAMS_t *params,
                npy_intp M_in, npy_intp P_in, npy_intp R_in)
{
    npy_uint8 *mem_buff = NULL;
    npy_uint8 *a, *b, *c, *d, *e, *f;
    fortran_int M, P, R, ldz;
    size_t safe_M, safe_P, safe_R;

    M = (fortran_int)M_in;
    P = (fortran_int)P_in;
    R = (fortran_int)R_in;
    safe_M = M_in;
    safe_P = P_in;
    safe_R = R_in;
    ldz = fortran_int_max(M, 1);

    size_t Espace, Xspace, Yspace, ZZspace, ZTMPspace, UZZspace, USspace, Zspace;

    Espace = sizeof(fortran_doublereal);
    Xspace = safe_M * Espace;
    Yspace = Xspace;
    USspace = Xspace;
    Zspace = safe_M * Xspace;
    UZZspace = safe_P * Zspace;
    ZZspace = safe_R * UZZspace;
    ZTMPspace = ZZspace;

    mem_buff = malloc(Espace + Xspace + Yspace + ZZspace + ZTMPspace + UZZspace + USspace);
    if (!mem_buff) {
        goto error;
    }
    a = mem_buff + Espace;  /* E at start of buffer, X after space for E */
    b = a + Xspace;  /* Y after space for X */
    c = b + Yspace;  /* ZZ after space for Y */
    d = c + ZZspace;  /* ZTMP after space for ZZ */
    e = d + ZTMPspace;  /* UZZ after space for ZTMP */
    f = e + UZZspace;  /* US after space for UZZ */

    params->TRANSZ = 'N';
    params->A = &d_one;  /* Y -> A*Z.X + B*Y, need A=1, B=0 */
    params->B = &d_zero;
    params->E = mem_buff;
    params->X = a;
    params->Y = b;
    params->ZZ = c;
    params->ZTMP = d;
    params->UZZ = e;
    params->US = f;
    params->Z = c;
    params->UZ = d;
    params->IP = Zspace;
    params->IR = UZZspace;
    params->K = 0;
    params->M = M;
    params->P = P;
    params->R = R;
    params->RPM = R * P * M;
    params->PM = P * M;
    params->PMM = P * M * M;
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
release_DOUBLE_m_init(GBMV_PARAMS_t *params)
{
    /* memory block base is in E */
    free(params->E);
    memset(params, 0, sizeof(*params));
}

/************************************
* Move parameters                   *
*************************************/

static NPY_INLINE int
choose_updater_inout_n(GBMV_PARAMS_t *params, npy_int P_in, npy_int R_in)
{
    if (P_in < (npy_int)(params->P) && R_in < (npy_int)(params->R))
    {
        npy_uint8 *a;
        size_t shift = (size_t)P_in * params->IP + (size_t)R_in * params->IR;
        a = params->ZZ;
        params->Z = a + shift;
        a = params->ZTMP;
        params->UZ = a + shift;
        return 0;
    } else {
        // printf("[P_in: %d, P_max: %d. %d\n", P_in, params->P, P_in < params->P);
        // printf(" R_in: %d, R_max: %d. %d]\n", R_in, params->R, R_in < params->R);
        return 1;
    }
}

static NPY_INLINE int
choose_updater_inout(GBMV_PARAMS_t *params, const char *plasttype, const char *readout)
{
    npy_int P_in = *(npy_int *)plasttype;
    npy_int R_in = *(npy_int *)readout;
    return choose_updater_inout_n(params, P_in, R_in);
}

static NPY_INLINE int
choose_updater_out_n(GBMV_PARAMS_t *params, npy_int P_in)
{
    // npy_int P_in = *(npy_int *)plasttype;
    if (P_in < (npy_int)(params->P))
    {
        npy_uint8 *a;
        a = params->UZZ;
        params->UZ = a + (size_t)P_in * params->IP;
        return 0;
    } else {
        // printf("[P_in: %d, P_max: %d. %d]\n", P_in, params->P, P_in < params->P);
        return 1;
    }
}


/*****************************
* BLAS Inner GUfunc loop     *
******************************/

/* copy initial, updaters into struct for fortran */
static void
linearize_updater(const char *arg, npy_intp *dimensions, npy_intp *steps,
                   GBMV_PARAMS_t *params, const LINEARIZE_DATA_t *lin_data)
{
    const char *ip_updater = arg;  //  in-ptr: updater, readouts
    const char *ip_updatep = arg;  //  in-ptr: updater, plast types
    npy_int cnt_r, cnt_p;

    /* Copy updaters, initials */
    for (cnt_r = 0; cnt_r < dimensions[0]; cnt_r++)
    {
        ip_updatep = ip_updater;

        for (cnt_p = 0; cnt_p < dimensions[1]; cnt_p++)
        {
            choose_updater_inout_n(params, cnt_p, cnt_r);
            linearize_DOUBLE_matrix(params->Z, ip_updatep, lin_data);
            ip_updatep += steps[1];
        }
        ip_updater += steps[0];
    }
}

/* copy initial, updaters from struct for numpy */
static void
delinearize_plast(char *arg, npy_intp *dimensions, npy_intp step,
                   GBMV_PARAMS_t *params, const LINEARIZE_DATA_t *lin_data)
{
    char *op_trans = arg;  //  out-ptr: updater, plast types
    npy_int cnt_r, cnt_p;
    npy_uint8 *a;

    /* Gather updaters */
    // printf("Gather r\n");
    for (cnt_r = 0; cnt_r < dimensions[0]; cnt_r++)
    {
        choose_updater_inout_n(params, 0, cnt_r);
        call_dgbmv_uuz_z_uz(params);
    }

    /* Copy updaters, initials */
    // printf("copy m\n");
    for (cnt_p = 0; cnt_p < dimensions[1]; cnt_p++)
    {
        choose_updater_out_n(params, cnt_p);
        delinearize_DOUBLE_matrix(op_trans, params->UZ, lin_data);
        op_trans += step;
    }
    // printf("copy init\n");
}

/* update_m_init_signature = "(r,p,m,m),(tm),(t),(t,m),(t,m),(t)->(p,m,m),(m)"; */

static void
DOUBLE_m_init(char **args, npy_intp *dimensions, npy_intp *steps, npy_intp steady)
{
INIT_OUTER_LOOP_8

    npy_intp len_r = dimensions[0];  /* number of readouts */
    npy_intp len_p = dimensions[1];  /* number of matrices */
    npy_intp len_m = dimensions[2];  /* number of states */
    npy_intp len_tm = dimensions[3]; /* number of time-points - 1 */
    npy_intp len_t = dimensions[4];  /* number of time-points */
    npy_intp stride_u_r = *steps++;  /* 1st arg */
    npy_intp stride_u_p = *steps++;
    npy_intp stride_u_mr = *steps++;
    npy_intp stride_u_mc = *steps++;
    npy_intp stride_p_t = *steps++;  /* 2nd arg */
    npy_intp stride_r_t = *steps++;  /* 3rd arg */
    npy_intp stride_a_t = *steps++;  /* 4th arg */
    npy_intp stride_a_m = *steps++;
    npy_intp stride_b_t = *steps++;  /* 5th arg */
    npy_intp stride_b_m = *steps++;
    npy_intp stride_e_t = *steps++;  /* 6th arg */
    npy_intp stride_uu_p = *steps++;  /* 1st output */
    npy_intp stride_uu_mr = *steps++;
    npy_intp stride_uu_mc = *steps++;
    npy_intp stride_ss_m = *steps++;  /* 2nd output */
    int error_occurred = get_fp_invalid_and_clear();
    GBMV_PARAMS_t params;
    LINEARIZE_DATA_t u_in, a_in, b_in, u_out, s_out;
    npy_intp cnt_t;

    if (len_t != len_tm + 1)
    {
        error_occurred = 1;
    }

    /* allocate buffer */
    // printf("allocate memory\n");
    if (error_occurred == 0 &&
        init_DOUBLE_m_init(&params, len_m, len_p, len_r))
    {
        /* initialise size parameters */
        init_linearize_data(&u_in, len_m, len_m, stride_u_mc, stride_u_mr);
        init_linearize_vdata(&a_in, len_m, stride_a_m);
        init_linearize_vdata(&b_in, len_m, stride_b_m);
        init_linearize_data(&u_out, len_m, len_m, stride_uu_mc, stride_uu_mr);
        init_linearize_vdata(&s_out, len_m, stride_ss_m);
        // printf("linearize updaters\n");


        /* Copy updaters, assume no broadcasting */
        npy_intp rsteps[] = {stride_u_r, stride_u_p};
        linearize_updater(args[0], dimensions, rsteps, &params, &u_in);
        // printf("updaters linearized\n");

        BEGIN_OUTER_LOOP

            const char *ip_plast = args[1];  //  in-ptr: plast_type
            const char *ip_reado = args[2];  //  in-ptr: readouts
            const char *ip_alpha = args[3];  //  in-ptr: alpha
            const char *ip_beta = args[4];   //  in-ptr: beta
            const char *ip_eta = args[5];    //  in-ptr: eta
            char *op_plast = args[6];    //  out-ptr: M
            char *op_initial = args[7];   //  out-ptr: init

            call_dlaset_u_zero(&params);
            linearize_DOUBLE_vec(params.X, ip_alpha, &a_in);
            linearize_DOUBLE_vec(params.Y, ip_beta, &b_in);
            call_dgbmv_us_x_y(&params);

            // printf("begin t loop\n");
            for (cnt_t = 0; cnt_t < len_tm; cnt_t++)
            {
                ip_beta += stride_b_t;
                ip_eta += stride_e_t;
                ip_reado += stride_r_t;
                // printf("copy b\n");
                linearize_DOUBLE_vec(params.Y, ip_beta, &b_in);
                // printf("copy e\n");
                *(fortran_doublereal *)(params.E) = *(npy_double *)ip_eta;

                // printf("choose updaters\n");
                error_occurred = choose_updater_inout(&params, ip_plast, ip_reado);
                if (error_occurred)
                    { break; }
                // printf("calc coeffs\n");
                call_dger_uz_x_y(&params);

                ip_alpha += stride_a_t;
                ip_plast += stride_p_t;
                // printf("copy a\n");
                linearize_DOUBLE_vec(params.X, ip_alpha, &a_in);

                // printf("calc init\n");
                /* if steady-state */
                if (steady)
                {
                    call_dgbmv_us_x_y(&params);
                }

            }
            if (error_occurred)
                { break; }
            // printf("end t loop\n");

            delinearize_plast(op_plast, dimensions, stride_uu_p, &params, &u_out);
            delinearize_DOUBLE_vec(op_initial, params.US, &s_out);
            // printf("copied out\n");
        END_OUTER_LOOP_8
        // printf("end outer loop\n");
        /* deallocate buffer */
        release_DOUBLE_m_init(&params);
        // printf("memory freed\n");
    }
    set_fp_invalid_or_clear(error_occurred);
}

/* update_m_init_signature = "(r,p,m,m),(tm),(t),(t,m),(t,m),(t)->(p,m,m),(m)"; */

static void
DOUBLE_plast_init(char **args, npy_intp *dimensions, npy_intp *steps,
              void *NPY_UNUSED(func))
{
    DOUBLE_m_init(args, dimensions, steps, 0);
}


/* update_m_init_signature = "(r,p,m,m),(tm),(t),(t,m),(t,m),(t)->(p,m,m),(m)"; */

static void
DOUBLE_plast_steady(char **args, npy_intp *dimensions, npy_intp *steps,
              void *NPY_UNUSED(func))
{
    DOUBLE_m_init(args, dimensions, steps, 1);
}

/*
*****************************************************************************
**                             Ufunc definition                            **
*****************************************************************************
*/

/* types argument for creating 'norm' ufunc */
static char ufn_types_1_7[] = {NPY_DOUBLE, NPY_DOUBLE, NPY_INT, NPY_INT,
                                NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};
static char ufn_types_1_8[] = {NPY_DOUBLE, NPY_INT, NPY_INT,
                                NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
                                NPY_DOUBLE, NPY_DOUBLE};

/* array of functions for each ufunc loop */
static PyUFuncGenericFunction
alpha_beta_funcs[] = {
    DOUBLE_alpha_beta
};
static PyUFuncGenericFunction
plast_init_funcs[] = {
    DOUBLE_plast_init
};
static PyUFuncGenericFunction
plast_steady_funcs[] = {
    DOUBLE_plast_steady
};

/* info for creating ufunc object */
GUFUNC_DESCRIPTOR_t gufunc_descriptors[] = {
    {"alpha_beta", "(r,p,m,m),(r,m),(tm),(t)->(t,m),(t,m),(t)", alpha_beta__doc__,
        1, 4, 3, alpha_beta_funcs, ufn_types_1_7},
    {"plast_init", "(r,p,m,m),(tm),(t),(t,m),(t,m),(t)->(p,m,m),(m)", plast_init__doc__,
        1, 6, 2, plast_init_funcs, ufn_types_1_8},
    {"plast_steady", "(r,p,m,m),(tm),(t),(t,m),(t,m),(t)->(p,m,m),(m)", plast_steady__doc__,
        1, 6, 2, plast_steady_funcs, ufn_types_1_8},
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
    int failure = addUfuncs(m, gufunc_descriptors, 3, baum_welch_version);

    if (PyErr_Occurred() || failure) {
        PyErr_SetString(PyExc_RuntimeError,
                        "cannot load _baum_welch module.");
        return NULL;
    }

    return m;
}

