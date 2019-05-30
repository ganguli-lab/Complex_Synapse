# -*- coding: utf-8 -*-
"""Utilities for Markov processes
"""
import numpy as np
import numpy_linalg as la

from sl_py_tools.numpy_tricks import allfinite, tri_low_rank


def stochastify_c(mat: la.lnarray):  # make cts time stochastic
    """
    Make a matrix the generator of a continuous time Markov process.
    Changes diagonal to make row sums zero.
    **Modifies** in place, **does not** return.

    Parameters
    ----------
    mat : la.lnarray (...,n,n)
        square matrix with non-negative off-diagonal elements.
        **Modified** in place.
    """
    mat -= np.apply_along_axis(np.diagflat, -1, mat.sum(axis=-1))


def stochastify_d(mat: la.lnarray):  # make dsc time stochastic
    """
    Make a matrix the generator of a discrete time Markov process.
    Scales rows to make row sums one.

    Parameters
    ----------
    mat : la.lnarray (...,n,n)
        square matrix with non-negative elements.
        **Modified** in place
    """
    mat /= mat.sum(axis=-1, keepdims=True)


def isstochastic_c(mat: la.lnarray, thresh: float = 1e-5) -> bool:
    """Are row sums zero?
    """
    nonneg = mat.flattish(-2) >= 0
    nonneg[..., ::mat.shape[-1]+1] = True
    return nonneg.all() and (np.fabs(mat.sum(axis=-1)) < thresh).all()


def isstochastic_d(mat: la.lnarray, thresh: float = 1e-5) -> bool:
    """Are row sums one?
    """
    return (np.fabs(mat.sum(axis=-1) - 1) < thresh).all() and (mat >= 0).all()


def calc_peq(rates):
    """Calculate steady state distribution.

    Parameters
    ----------
    rates : np.ndarray (n,n) or Tuple[np.ndarray (n,n), np.ndarray (n,)]
        Continuous time stochastic matrix or LU factors of inverse fundamental.

    Returns
    -------
    peq : la.lnarray (n,)
        Steady-state distribution.
    (z_lu, ipv) : np.ndarray (n,n),(n,)
        LU factors of inverse fundamental matrix.
    """
    if isinstance(rates, tuple):
        z_lu, ipv = rates
        evc = la.ones(z_lu[0].shape[0])
        peq = la.gufuncs.rlu_solve(evc.r, z_lu, ipv)
    else:
        evc = la.ones(rates.shape[0])
        fund_inv = la.ones_like(rates) - rates
        peq, z_lu, ipv = la.gufuncs.rsolve_lu(evc.r, fund_inv)
    # check for singular matrix
    if not allfinite(z_lu) or tri_low_rank(z_lu):
        return la.full_like(evc, np.nan)

    return peq, (z_lu, ipv)


def adjoint(tensor: la.lnarray, measure: la.lnarray) -> la.lnarray:
    """Adjoint with respect to L2 inner product with measure

    Parameters
    ----------
    tensor : la.lnarray (...,n,n) or (...,1,n) or (...,n,1)
        The matrix/row/column vector to be adjointed.
    measure : la.lnarray (...,n)
        The measure for the inner-product wrt which we adjoint

    Parameters
    ----------
    tensor : la.lnarray (...,n,n) or (...,n,1) or (...,1,n)
        The adjoint matrix/column/row vector.
    """
    adj = tensor.t
    if adj.shape[-1] == 1:  # row -> col
        adj /= measure.c
    elif adj.shape[-2] == 1:  # col -> row
        adj *= measure.r
    else:  # mat -> mat
        adj *= measure.r / measure.c
    return adj


def rand_trans(nst: int, npl: int = 1, sparsity: float = 1.) -> la.lnarray:
    """
    Make a random transition matrix (continuous time).

    Parameters
    ----------
    n : int
        total number of states
    npl : int
        number of plasticity types
    sparsity : float
        sparsity

    Returns
    -------
    mat : la.lnarray
        transition matrix
    """
    mat = la.random.rand(npl, nst, nst)
    ind = la.random.rand(npl, nst, nst)
    mat[ind > sparsity] = 0.
    stochastify_c(mat)
    return mat


def sim_markov(rates, peq, num_jmp=None, max_time=None):
    """Simulate Markov process trajectory.

    Parameters
    ----------
    rates : la.lnarray (n,n)
        Continuous time stochastic matrix.
    peq : la.lnarray (n,)
        Initial-state distribution.
    num_jmp : int, optional, default: None
        Stop after this many jumps.
    max_time : float, optional, default: None
        Stop after this much time.

    Returns
    -------
    states : la.lnarray (w,)
        Vector of states visited.
    dwells : la.lnarray (w,)
    """
    num_states = len(peq)
    state_inds = la.arange(num_states)
    dwell = -1. / np.diagonal(rates)
    jump = rates * dwell.c
    jump[np.diag_indices(num_states)] = 0.

    if num_jmp is None:
        if max_time is None:
            raise ValueError("Must specify either num_jmp or max_time")
        num_jmp = np.inf
        est_num = int(5 * max_time * (peq / dwell).sum())
    if max_time is None:
        max_time = np.inf
        est_num = num_jmp

    states_from = la.array([np.random.choice(state_inds, size=est_num-1, p=p)
                            for p in jump])
    dwells_from = - dwell.c * np.log(la.random.rand(est_num))
    states = [np.random.choice(state_inds, p=peq)]
    dwells = [dwells_from[states[-1], 0]]
    num = 1
    time = dwells[-1]
    while num < num_jmp and time < max_time:
        states.append(states_from[states[-1], num - 1])
        dwells.append(dwells_from[states[-1], num])
        num += 1
        time += dwells[-1]
        # if num >= est_num:
        #     msg = f"n/N/N^ {num}/{num_jump}/{est_num}, t/T {time}/{max_time}"
        #     raise IndexError(msg)
    if time > max_time:
        dwells[-1] -= time - max_time
    return la.array(states), la.array(dwells)
