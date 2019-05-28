"""Utilities for Markov processes
"""
import numpy as np
import numpy_linalg as la

from sl_py_tools.numpy_tricks import allfinite, tri_low_rank


@la.wrappers.wrap_one
def offdiag_inds(nst):
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
    k_1st = (nst+1) * la.arange(nst)  # ravel ind of 1st element in group
    k = la.arange(nst**2)
    return la.delete(k, k_1st)


@la.wrappers.wrap_several
def offdiag_subs(nst):
    """Indices of independent parameters of transition matrix.

    Parameters
    ----------
    nst : int
        Number of states.

    Returns
    -------
    ind0, ind1 : la.lnarray (n(n-1),)
        Vectors of indices of off-diagonal elements, in order:
        mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2,n-1.
    """
    rav_inds = offdiag_inds(nst)
    ind0, ind1 = np.indices((nst, nst), int)
    return ind0.ravel()[rav_inds], ind1.ravel()[rav_inds]


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
    return ((np.fabs(mat.sum(axis=-1)) < thresh).all()
            and (mat[(...,) + offdiag_subs(mat.shape[-1])] > 0).all())


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
    """
    if isinstance(rates, tuple):
        z_lu, ip = rates
        ev = la.ones(z_lu[0].shape[0])
        peq = la.gufuncs.rlu_solve(ev.r, z_lu, ip)
    else:
        ev = la.ones(rates.shape[0])
        fund_inv = la.ones_like(rates) - rates
        peq, z_lu, ip = la.gufuncs.rsolve_lu(ev.r, fund_inv)
    # check for singular matrix
    if not allfinite(z_lu) or tri_low_rank(z_lu):
        return la.full_like(ev, np.nan)

    return peq


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


@la.wrappers.wrap_one
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


def num_param(states, ring=False, uniform=False):
    """Number of independent rates

    Parameters
    ----------
    states : int or la.lnarray (n,...)
        Number of states, or array over states.
    ring : bool, optional, default: True
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
    if ring:
        if uniform:
            return 2
        return 2 * states
    return states * (states - 1)


def num_state(params, ring=False, uniform=False):
    """Number of states from rate vector

    Parameters
    ----------
    params : int or np.ndarray (n,)
        Number of rate parameters, or vector of rates.
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
    if ring:
        if uniform:
            return np.nan
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
    ring : bool, optional, default: True
        Is the rate vector meant for `ring_params_to_mat` or
        `gen_params_to_mat`?
    uniform : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `uni_ring_params_to_mat`?
    """
    if isinstance(params, np.ndarray):
        params = params.shape[-1]
    if isinstance(states, np.ndarray):
        states = states.shape[-1]
    uniform = (params == 2)
    ring = uniform or (params == 2 * states)
    return ring, uniform


def gen_params_to_mat(params, *args, **kwds):
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
    nst = num_state(params, ring=False)
    mat = la.empty(nst**2)
    mat[offdiag_inds(nst)] = params
    mat.resize((nst, nst))
    stochastify_c(mat)
    return mat


def gen_mat_to_params(mat, *args, **kwds):
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


def ring_params_to_mat(params, *args, **kwds):
    """Ring transition matrix from independent parameters.

    Parameters
    ----------
    params : np.ndarray (2n,)
        Vector of independent elements, in order:
        mat_01, mat_12, ..., mat_n-2,n-1, mat_n-10,
        mat_0n-1, mat10, mat_21, ..., mat_n-1,n-2.

    Returns
    -------
    mat : la.lnarray (n,n)
        Continuous time stochastic matrix.
    """
    nst = num_state(params.size, ring=True, uniform=False)
    above, bottom, top, below = la.split(params, [nst-1, nst, nst+1])
    mat = la.diagflat(above, 1) + la.diagflat(bottom, 1-nst)
    mat += la.diagflat(top, nst-1) + la.diagflat(below, -1)
    stochastify_c(mat)
    return mat


def ring_mat_to_params(mat, *args, **kwds):
    """Independent parameters of ring transition matrix.

    Parameters
    ----------
    mat : np.ndarray (n,n)
        Continuous time stochastic matrix.

    Returns
    -------
    params : la.lnarray (2n,)
        Vector of independent elements, in order:
        mat_01, mat_12, ..., mat_n-2n-1, mat_n-10,
        mat_0n-1, mat10, mat_21, ..., mat_n-2,n-1.
    """
    nst = mat.shape[0]
    params = (np.diagonal(mat, 1), np.diagonal(mat, 1-nst),
              np.diagonal(mat, nst-1), np.diagonal(mat, -1))
    return la.hstack(params)


def uni_ring_params_to_mat(params, num_st, *args, **kwds):
    """Ring transition matrix from independent parameters.

    Parameters
    ----------
    params : np.ndarray (2,)
        Vector of independent elements, in order:
        mat_01 = mat_12 = ... = mat_n-2,n-1 = mat_n-10,
        mat_0n-1 = mat10 = mat_21 = ... = mat_n-1,n-2.
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


def uni_ring_mat_to_params(mat, grad=True, *args, **kwds):
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
    params = la.array([np.sum(ring_params[:nst]), np.sum(ring_params[nst:])])
    if not grad:
        params /= nst
    return params


def params_to_mat(params, ring=False, uniform=False, num_st=2, **kwds):
    """Transition matrix from independent parameters.

    Parameters
    ----------
    params : np.ndarray (n(n-1),)
        Vector of independent elements, in order:
        mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2,n-1.
    ring : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `gen_params_to_mat`?
    uniform : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `uni_ring_params_to_mat`?

    Returns
    -------
    mat : la.lnarray (n,n)
        Continuous time stochastic matrix.
    """
    if ring:
        if uniform:
            return uni_ring_params_to_mat(params, num_st, **kwds)
        return ring_params_to_mat(params, **kwds)
    return gen_params_to_mat(params, **kwds)


def mat_to_params(mat, ring=False, uniform=False, grad=True, **kwds):
    """Independent parameters of transition matrix.

    Parameters
    ----------
    mat : np.ndarray (n,n)
        Continuous time stochastic matrix.
    ring : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `gen_params_to_mat`?
    uniform : bool, optional, default: False
        Is the rate vector meant for `ring_params_to_mat` or
        `uni_ring_params_to_mat`?

    Returns
    -------
    params : la.lnarray (n(n-1),)
        Vector of independent elements, in order:
        mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2,n-1.
    """
    if ring:
        if uniform:
            return uni_ring_mat_to_params(mat, grad, **kwds)
        return ring_mat_to_params(mat, **kwds)
    return gen_mat_to_params(mat, **kwds)


def sim_markov(rates: la.lnarray, peq: la.lnarray, num_jump=None,
               max_time=None, eps=0.001):
    """Simulate Markov process trajectory.

    Parameters
    ----------
    rates : la.lnarray (n,n)
        Continuous time stochastic matrix.
    peq : la.lnarray (n,)
        Initial-state distribution.
    num_jump : int, optional, default: None
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

    if num_jump is None:
        if max_time is None:
            raise ValueError("Must specify either num_jump or max_time")
        num_jump = np.inf
        est_num = int(5 * max_time * (peq / dwell).sum())
    if max_time is None:
        max_time = np.inf
        est_num = num_jump

    states_from = la.array([np.random.choice(state_inds, size=est_num-1, p=p)
                            for p in jump])
    dwells_from = - dwell.c * np.log(la.random.rand(est_num))
    states = [np.random.choice(state_inds, p=peq)]
    dwells = [dwells_from[states[-1], 0]]
    num = 1
    time = dwells[-1]
    while num < num_jump and time < max_time:
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
