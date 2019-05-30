# -*- coding: utf-8 -*-
"""Utilities for Markov processes
"""
import numpy as np
import numpy_linalg as la
from .markov import stochastify_c


def offdiag_inds(nst: int) -> la.lnarray:
    """Ravel indices of independent parameters of transition matrix.

    Parameters
    ----------
    nst : int
        Number of states.

    Returns
    -------
    K : la.lnarray (n(n-1),)
        Vector of ravel indices of off-diagonal elements, in order:
        mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2,n-1.
    """
    # when put into groups of size nst+1:
    # M[0,0], M[0,1], M[0,2], ..., M[0,n-1], M[1,0],
    # M[1,1], M[1,2], ..., M[1,n-1], M[2,0], M[2,1],
    # ...
    # M[n-2,n-2], M[n-2,n-1], M[n-1,0], ..., M[n-1,n-2],
    # unwanted elements are 1st in each group
    k_1st = la.arange(0, nst, nst+1)  # ravel ind of 1st element in group
    k = la.arange(nst**2)
    return la.delete(k, k_1st)


def serial_inds(nst: int) -> la.lnarray:
    """Ravel indices of independent parameters of serial transition matrix.

    Parameters
    ----------
    nst : int
        Number of states.

    Returns
    -------
    K : la.lnarray (2(n-1),)
        Vector of ravel indices of off-diagonal elements, in order:
        mat_01, mat_12, ..., mat_n-2,n-1,
        mat_10, mat_21, ..., mat_n-2,n-1.
    """
    return la.hstack((np.arange(1, nst**2, nst+1),
                      np.arange(nst, nst**2, nst+1)))


def ring_inds(nst: int) -> la.lnarray:
    """Ravel indices of independent parameters of ring transition matrix.

    Parameters
    ----------
    nst : int
        Number of states.

    Returns
    -------
    K : la.lnarray (2n,)
        Vector of ravel indices of off-diagonal elements, in order:
        mat_01, mat_12, ..., mat_n-2,n-1, mat_n-1,0,
        mat_0,n-1, mat_10, mat_21, ..., mat_n-1,n-2.
    """
    return la.hstack((np.arange(1, nst**2, nst+1), [nst*(nst-1), nst-1],
                      np.arange(nst, nst**2, nst+1)))


def param_inds(nst: int, serial=False, ring=False) -> la.lnarray:
    """Ravel indices of independent parameters of transition matrix.

    Parameters
    ----------
    nst : int
        Number of states.
    serial : bool, optional, default: False
        Is the rate vector meant for `serial_params_to_mat` or
        `gen_params_to_mat`?
    ring : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `gen_params_to_mat`?
    uniform : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `uni_ring_params_to_mat`?

    Returns
    -------
    mat : la.lnarray (k,k), k in (n(n-1), 2(n-1), 2n, 2)
        Matrix of independent elements, each axis in order of `mat_to_params`.
    """
    if serial:
        return serial_inds(nst)
    if ring:
        return ring_inds(nst)
    return offdiag_inds(nst)


def num_param(states, serial=False, ring=False, uniform=False):
    """Number of independent rates

    Parameters
    ----------
    states : int or la.lnarray (n,...)
        Number of states, or array over states.
    serial : bool, optional, default: False
        Is the rate vector meant for `serial_params_to_mat` or
        `gen_params_to_mat`?
    ring : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `gen_params_to_mat`?
    uniform : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `uni_ring_params_to_mat`?

    Returns
    -------
    params : int
        Number of rate parameters.
    """
    if isinstance(states, np.ndarray):
        states = states.shape[0]
    if uniform:
        return 2
    if serial:
        return 2 * states - 2
    if ring:
        return 2 * states
    return states * (states - 1)


def num_state(params, serial=False, ring=False):
    """Number of states from rate vector

    Parameters
    ----------
    params : int or np.ndarray (n,)
        Number of rate parameters, or vector of rates.
    serial : bool, optional, default: False
        Is the rate vector meant for `serial_params_to_mat` or
        `gen_params_to_mat`?
    ring : bool, optional, default: True
        Is the rate vector meant for `ring_params_to_mat` or
        `gen_params_to_mat`?
    uniform : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `uni_ring_params_to_mat`?

    Returns
    -------
    states : int
        Number of states.
    """
    if isinstance(params, np.ndarray):
        params = params.shape[-1]
    if serial:
        return params // 2 + 1
    if ring:
        return params // 2
    return (0.5 + np.sqrt(0.25 + params)).astype(int)


def mat_type(params, states):
    """Is it a (uniform) ring

    Parameters
    ----------
    params : int or np.ndarray (np,)
        Number of rate parameters, or vector of rates.
    states : int or np.ndarray (n,...)
        Number of states, or array over states.

    Returns
    -------
    serial : bool
        Is the rate vector meant for `serial_params_to_mat` or
        `gen_params_to_mat`?
    ring : bool
        Is the rate vector meant for `ring_params_to_mat` or
        `gen_params_to_mat`?
    uniform : bool
        Is the rate vector meant for `ring_params_to_mat` or
        `uni_ring_params_to_mat`?
    """
    if isinstance(params, np.ndarray):
        params = params.shape[-1]
    if isinstance(states, np.ndarray):
        states = states.shape[-1]
    uniform = (params == 2)
    ring = uniform or (params == 2 * states)
    serial = uniform or (params == 2 * states - 2)
    return serial, ring, uniform


def gen_params_to_mat(params):
    """Transition matrix from independent parameters.

    Parameters
    ----------
    params : np.ndarray (n(n-1),)
        Vector of off-diagonal elements, in order:
        mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2,n-1.

    Returns
    -------
    mat : la.lnarray (n,n)
        Continuous time stochastic matrix.
    """
    nst = num_state(params)
    mat = la.empty(nst**2)
    mat[offdiag_inds(nst)] = params
    mat.resize((nst, nst))
    stochastify_c(mat)
    return mat


def serial_params_to_mat(params):
    """Serial transition matrix from independent parameters.

    Parameters
    ----------
    params : np.ndarray (2(n-1),)
        Vector of independent elements, in order:
        mat_01, mat_12, ..., mat_n-2,n-1,
        mat_10, mat_21, ..., mat_n-1,n-2.

    Returns
    -------
    mat : la.lnarray (n,n)
        Continuous time stochastic matrix.
    """
    njmp = num_state(params, serial=True) - 1
    mat = la.diagflat(params[:njmp], 1) + la.diagflat(params[njmp:], -1)
    stochastify_c(mat)
    return mat


def uni_serial_params_to_mat(params, num_st):
    """Uniform serial transition matrix from independent parameters.

    Parameters
    ----------
    params : np.ndarray (2,)
        Vector of independent elements, in order:
        mat_01 = mat_12 = ... = mat_n-2,n-1,
        mat_10 = mat_21 = ... = mat_n-1,n-2.
    num_st : int
        Number of states.

    Returns
    -------
    mat : la.lnarray (n,n)
        Continuous time stochastic matrix.
    """
    serial_params = la.hstack((np.full(num_st-1, params[0]),
                               np.full(num_st-1, params[1])))
    return serial_params_to_mat(serial_params)


def ring_params_to_mat(params):
    """Ring transition matrix from independent parameters.

    Parameters
    ----------
    params : np.ndarray (2n,)
        Vector of independent elements, in order:
        mat_01, mat_12, ..., mat_n-2,n-1, mat_n-1,0,
        mat_0,n-1, mat_10, mat_21, ..., mat_n-1,n-2.

    Returns
    -------
    mat : la.lnarray (n,n)
        Continuous time stochastic matrix.
    """
    nst = num_state(params.size, ring=True)
    mat = la.diagflat(params[:nst-1], 1) + la.diagflat(params[nst+1:], -1)
    mat[nst-1, 0] = params[nst-1]
    mat[0, nst-1] = params[nst]
    stochastify_c(mat)
    return mat


def uni_ring_params_to_mat(params, num_st):
    """Ring transition matrix from independent parameters.

    Parameters
    ----------
    params : np.ndarray (2,)
        Vector of independent elements, in order:
        mat_01 = mat_12 = ... = mat_n-2,n-1 = mat_n-1,0,
        mat_0,n-1 = mat_10 = mat_21 = ... = mat_n-1,n-2.
    num_st : int
        Number of states.

    Returns
    -------
    mat : la.lnarray (n,n)
        Continuous time stochastic matrix.
    """
    ring_params = la.hstack((np.full(num_st, params[0]),
                             np.full(num_st, params[1])))
    return ring_params_to_mat(ring_params)


def params_to_mat(params, serial=False, ring=False, uniform=False, nst=2):
    """Transition matrix from independent parameters.

    Parameters
    ----------
    params : np.ndarray (n(n-1),) or (2(n-1),) or (2n,) or (2,)
        Vector of independent elements, in order that depends on flags below.
    serial : bool, optional, default: False
        Is the rate vector meant for `serial_params_to_mat` or
        `gen_params_to_mat`?
    ring : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `gen_params_to_mat`?
    uniform : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `uni_ring_params_to_mat`?
    nst : int, optional, default: 2
        Number of states. Only needed when `uniform` is True

    Returns
    -------
    mat : la.lnarray (n,n)
        Continuous time stochastic matrix.
    """
    if serial:
        if uniform:
            return uni_serial_params_to_mat(params, nst)
        return serial_params_to_mat(params)
    if ring:
        if uniform:
            return uni_ring_params_to_mat(params, nst)
        return ring_params_to_mat(params)
    return gen_params_to_mat(params)


def gen_mat_to_params(mat):
    """Independent parameters of transition matrix.

    Parameters
    ----------
    mat : np.ndarray (n,n)
        Continuous time stochastic matrix.

    Returns
    -------
    params : la.lnarray (n(n-1),)
        Vector of off-diagonal elements, in order:
        mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2,n-1.
    """
    nst = mat.shape[0]
    # m_00, m_01, m_02, ..., m_0n-1, m_10,
    # mat_12, ..., mat_n-2,n
    param = mat.flatten()
    return param[offdiag_inds(nst)]


def serial_mat_to_params(mat):
    """Independent parameters of serial transition matrix.

    Parameters
    ----------
    mat : np.ndarray (n,n)
        Continuous time stochastic matrix.

    Returns
    -------
    params : la.lnarray (2(n-1),)
        Vector of independent elements, in order:
        mat_01, mat_12, ..., mat_n-2,n-1,
        mat_10, mat_21, ..., mat_n-2,n-1.
    """
    nst = mat.shape[0]
    param = mat.flatten()
    return param[serial_inds(nst)]


def uni_serial_mat_to_params(mat, grad=True):
    """Independent parameters of uniform serial transition matrix.

    Parameters
    ----------
    mat : np.ndarray (n,n)
        Continuous time stochastic matrix.
    grad : bool, optional, default: True
        Is the output for a gradient (True) or a transition matrix (False).
        If True, return sum of (anti)clockwise transitions.
        If False, return the mean.

    Returns
    -------
    params : la.lnarray (2,)
        Vector of independent elements, in order (grad=False):
            mat_01 = mat_12 = ... = mat_n-2,n-1,
            mat_10 = mat_21 = ... = mat_n-1,n-2.
        Or, in order (grad=True):
            mat_01 + mat_12 + ... + mat_n-2,n-1,
            mat_10 + mat_21 + ... + mat_n-1,n-2.
    """
    njmp = mat.shape[0] - 1
    serl_params = serial_mat_to_params(mat)
    params = la.hstack([serl_params[:njmp].sum(), serl_params[njmp:].sum()])
    if not grad:
        params /= njmp
    return params


def ring_mat_to_params(mat):
    """Independent parameters of ring transition matrix.

    Parameters
    ----------
    mat : np.ndarray (n,n)
        Continuous time stochastic matrix.

    Returns
    -------
    params : la.lnarray (2n,)
        Vector of independent elements, in order:
        mat_01, mat_12, ..., mat_n-2,n-1, mat_n-1,0,
        mat_0,n-1, mat_10, mat_21, ..., mat_n-1,n-2.
    """
    nst = mat.shape[0]
    param = mat.flatten()
    return param[ring_inds(nst)]


def uni_ring_mat_to_params(mat, grad=True):
    """Independent parameters of ring transition matrix.

    Parameters
    ----------
    mat : np.ndarray (n,n)
        Continuous time stochastic matrix.
    grad : bool, optional, default: True
        Is the output for a gradient (True) or a transition matrix (False).
        If True, return sum of (anti)clockwise transitions.
        If False, return the mean.

    Returns
    -------
    params : la.lnarray (2,)
        Vector of independent elements, in order (grad=False):
            mat_01 = mat_12 = ... = mat_n-2,n-1 = mat_n-10,
            mat_0n-1 = mat10 = mat_21 = ... = mat_n-1,n-2.
        Or, in order (grad=True):
            mat_01 + mat_12 + ... + mat_n-2,n-1 + mat_n-10,
            mat_0n-1 + mat10 + mat_21 + ... + mat_n-1,n-2.
    """
    nst = mat.shape[0]
    ring_params = ring_mat_to_params(mat)
    params = la.array([ring_params[:nst].sum(), ring_params[nst:].sum()])
    if not grad:
        params /= nst
    return params


def mat_to_params(mat, serial=False, ring=False, uniform=False, grad=True):
    """Independent parameters of transition matrix.

    Parameters
    ----------
    mat : np.ndarray (n,n)
        Continuous time stochastic matrix.
    serial : bool, optional, default: False
        Is the rate vector meant for `serial_params_to_mat` or
        `gen_params_to_mat`?
    ring : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `gen_params_to_mat`?
    uniform : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `uni_ring_params_to_mat`?

    Returns
    -------
    params : la.lnarray (n(n-1),) or (2(n-1),) or (2n,) or (2,)
        Vector of independent elements. For the order, see `*_mat_to_params`.
    """
    if serial:
        if uniform:
            return uni_serial_mat_to_params(mat, grad)
        return serial_mat_to_params(mat)
    if ring:
        if uniform:
            return uni_ring_mat_to_params(mat, grad)
        return ring_mat_to_params(mat)
    return gen_mat_to_params(mat)


def tens_to_mat(tens, serial=False, ring=False, uniform=False, grad=True):
    """Independent parameters of 4th rank tensor.

    Parameters
    ----------
    tens : np.ndarray (n,n,n,n)
        Continuous time stochastic matrix.
    serial : bool, optional, default: False
        Is the rate vector meant for `serial_params_to_mat` or
        `gen_params_to_mat`?
    ring : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `gen_params_to_mat`?
    uniform : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `uni_ring_params_to_mat`?

    Returns
    -------
    mat : la.lnarray (n**2,n**2)
        Matrix of independent elements, each axis in order `*_mat_to_params`.
    """
    nst = tens.shape[0]
    mat = tens.reshape((nst**2, nst**2))
    inds = param_inds(nst, serial=serial, ring=ring)
    mat = mat[np.ix_(inds, inds)]
    if uniform:
        nind = len(inds) // 2
        mat = la.block([[mat[:nind, :nind].sum(), mat[:nind, nind:].sum()],
                        [mat[nind:, :nind].sum(), mat[nind:, nind:].sum()]])
        if not grad:
            mat /= nind
    return mat
