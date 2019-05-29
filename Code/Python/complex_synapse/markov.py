"""Utilities for Markov processes
"""
import numpy as np
import numpy_linalg as la

from sl_py_tools.numpy_tricks import allfinite, tri_low_rank


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
    valid = (mat.flattish(-2, -1)[..., offdiag_inds(mat.shape[-1])] >= 0).all()
    return valid and (np.fabs(mat.sum(axis=-1)) < thresh).all()


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


def num_state(params, serial=False, ring=False, uniform=False):
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
    if uniform:
        return np.nan
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


def serial_params_to_mat(params, *args, **kwds):
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
    above, below = np.split(params, 2)
    mat = la.diagflat(above, 1) + la.diagflat(below, -1)
    stochastify_c(mat)
    return mat


def serial_mat_to_params(mat, *args, **kwds):
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
    return la.hstack((np.diagonal(mat, 1), np.diagonal(mat, -1)))


def uni_serial_params_to_mat(params, num_st, *args, **kwds):
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


def uni_serial_mat_to_params(mat, grad=True, *args, **kwds):
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
    n = mat.shape[0]
    serl_params = serial_mat_to_params(mat)
    params = la.array([serl_params[:n-1].sum(), serl_params[n-1:].sum()])
    if not grad:
        params /= (n-1)
    return params


def ring_params_to_mat(params, *args, **kwds):
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
        mat_01, mat_12, ..., mat_n-2,n-1, mat_n-1,0,
        mat_0,n-1, mat_10, mat_21, ..., mat_n-1,n-2.
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
    params = la.array([ring_params[:nst].sum(), ring_params[nst:].sum()])
    if not grad:
        params /= nst
    return params


def params_to_mat(params, serial=False, ring=False, uniform=False, nst=2,
                  **kwds):
    """Transition matrix from independent parameters.

    Parameters
    ----------
    params : np.ndarray (n(n-1),) or (2(n-1),) or (2n,) or (2,)
        Vector of independent elements, in order:
        mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2,n-1.
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
            return uni_serial_params_to_mat(params, nst, **kwds)
        return serial_params_to_mat(params, **kwds)
    if ring:
        if uniform:
            return uni_ring_params_to_mat(params, nst, **kwds)
        return ring_params_to_mat(params, **kwds)
    return gen_params_to_mat(params, **kwds)


def mat_to_params(mat, serial=False, ring=False, uniform=False, grad=True,
                  **kwds):
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
            return uni_serial_mat_to_params(mat, grad, **kwds)
        return serial_mat_to_params(mat, **kwds)
    if ring:
        if uniform:
            return uni_ring_mat_to_params(mat, grad, **kwds)
        return ring_mat_to_params(mat, **kwds)
    return gen_mat_to_params(mat, **kwds)


def serial_inds(nst):
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
    return serial_mat_to_params(la.arange(nst**2).reshape(nst, nst))


def ring_inds(nst):
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
    return serial_mat_to_params(la.arange(nst**2).reshape(nst, nst))


def param_inds(nst, serial=False, ring=False, **kwds):
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
        return serial_inds(nst, **kwds)
    if ring:
        return ring_inds(nst, **kwds)
    return offdiag_inds(nst, **kwds)


def tens_to_mat(tens, serial=False, ring=False, uniform=False, grad=True,
                **kwds):
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
    k = param_inds(nst, serial=serial, ring=ring, **kwds)
    mat = mat[np.ix_(k, k)]
    if uniform:
        p = len(k) // 2
        mat = la.block([[mat[:p, :p].sum(), mat[:p, p:].sum()],
                        [mat[p:, :p].sum(), mat[p:, p:].sum()]])
        if not grad:
            mat /= p
    return mat
