# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 14:57:38 2017

Utility functions for constructing Markov transition matrices
and SynapseMemoryModel objects.

Methods
-------
binary_weights(nst)
    vector of +/-1.
linear_weights(nst)
    vector of numbers between +/-1.
stochastifyC(W)
    turn into continuous time transition matrix.
stochastifyD(mat)
    turn into discrete time transition matrix.
rand_trans(nst, sp)
    random transition matrix.
build_zero(nst)
    empty model (use with SynapseMemoryModel.Build).
build_rand(nst, sp)
    random model (use with SynapseMemoryModel.Build).
build_serial(nst, q)
    serial model (use with SynapseMemoryModel.Build).
build_multistate(nst, q)
    multistate model (use with SynapseMemoryModel.Build).

@author: Subhy
"""

from typing import Dict
import numpy as np
import numpy_linalg as la


def binary_weights(nst: int) -> la.lnarray:  # binary synaptic weights
    """
    Make binary (+/-1) weight vector

    Parameters
    ----------
    n : int
        total number of states

    Returns
    -------
    w : la.lnarray
        vector of synaptic weights
    """
    weight = la.ones(nst // 2)
    return la.hstack((-weight, weight))


def linear_weights(nst: int) -> la.lnarray:  # linear synaptic weights
    """
    Make linear (-1,...,1) weight vector

    Parameters
    ----------
    n : int
        total number of states

    Returns
    -------
    w : la.lnarray
        vector of synaptic weights
    """
    return la.linspace(-1., 1., nst)


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
    return (np.fabs(mat.sum(axis=-1)) < thresh).all()


def isstochastic_d(mat: la.lnarray, thresh: float = 1e-5) -> bool:
    """Are row sums one?
    """
    return (np.fabs(mat.sum(axis=-1) - 1.) < thresh).all()


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
    mat = np.random.rand(npl, nst, nst)
    ind = np.random.rand(npl, nst, nst)
    mat[ind > sparsity] = 0.
    stochastify_c(mat)
    return mat


def serial_trans(nst: int, jmp: float = 1.) -> la.lnarray:
    """
    Make a random transition matrix (continuous time).

    Parameters
    ----------
    n : int
        total number of states
    sparsity : float
        sparsity

    Returns
    -------
    mat : la.lnarray
        transition matrix
    """
    mat = la.stack((la.diag(jmp * np.ones(nst-1), 1),
                    la.diag(jmp * np.ones(nst-1), -1)))
    stochastify_c(mat)
    return mat


def build_generic(func, nst: int, npl: int = 2,
                  binary: bool = False) -> Dict[str, la.lnarray]:
    """Make a model from a matrix creating function.

    Factory for model builderss.
    For use with `SynapseMemoryModel.Build`.

    Parameters
    ----------
    Func : function
        function with signature `func(nst: int, npl: int) -> la.lnarray` and
        func(nst).shape = (npl, nst, nst)
    nst : int
        total number of states
    npl : optional, int = 2
        number of plasticity matrices
    binary : optional, bool = True
        is the weight vector binary? Otherwise it's linear. Default: False

    Returns
    -------
    dictionary:
        plast : la.lnarray
            potentiation/depression transition matrices
        weight : la.lnarray
            synaptic weights (linear/binary)
        signal : la.lnarray
            desired signal contribution from each plasticity type
    """
    plast = func(nst, npl)
    if binary:
        weight = binary_weights(nst)
    else:
        weight = linear_weights(nst)
    signal = la.linspace(1, -1, npl)
    return {'plast': plast, 'weight': weight, 'signal': signal}


def build_zero(nst: int, npl: int = 2,
               binary: bool = False) -> Dict[str, la.lnarray]:
    """Make an empty model.

    Make an empty model, i.e. all transition rates are zero.
    For use with `SynapseMemoryModel.Build`.

    Parameters
    ----------
    n : int
        total number of states
    binary : bool
        is the weight vector binary? Otherwise it's linear. Default: False

    Returns
    -------
    dictionary:
        plast : la.lnarray
            potentiation/depression transition matrices
        weight : la.lnarray
            synaptic weights (linear/binary)
        signal : la.lnarray
            desired signal contribution from each plasticity type
    """
    return build_generic(lambda n, p: la.zeros((p, n, n)),
                         nst, npl, binary)


def build_empty(nst: int, npl: int = 2,
                binary: bool = False) -> Dict[str, la.lnarray]:
    """Make an empty model.

    Make an empty model, i.e. all transition rates uninitialised.
    For use with `SynapseMemoryModel.Build`.

    Parameters
    ----------
    n : int
        total number of states
    binary : bool
        is the weight vector binary? Otherwise it's linear. Default: False

    Returns
    -------
    dictionary:
        plast : la.lnarray
            potentiation/depression transition matrices
        weight : la.lnarray
            synaptic weights (linear/binary)
        signal : la.lnarray
            desired signal contribution from each plasticity type
    """
    return build_generic(lambda n, p: la.empty((p, n, n)),
                         nst, npl, binary)


def build_rand(nst: int, npl: int = 2, binary: bool = False,
               sparsity: float = 1.) -> Dict[str, la.lnarray]:
    """Make a random model.

    Make a random model, i.e. transition rates are random numbers.
    For use with `SynapseMemoryModel.Build`.

    Parameters
    ----------
    n : int
        total number of states
    binary : bool
        is the weight vector binary? Otherwise it's linear. Default: False
    optional arguments
        passed to rand_trans
    sp : float
        sparsity, default: 1.

    Returns
    -------
    dictionary:
        plast : la.lnarray
            potentiation/depression transition matrices
        weight : la.lnarray
            synaptic weights (linear/binary)
        signal : la.lnarray
            desired signal contribution from each plasticity type
    """
    return build_generic(lambda n, p: rand_trans(n, p, sparsity),
                         nst, npl, binary)


def build_serial(nst: int, jmp: float) -> Dict[str, la.lnarray]:
    """Make a serial model.

    Make a serial model, i.e. uniform nearest neighbour transition rates,
    and binary weights.
    For use with `SynapseMemoryModel.Build`.

    Parameters
    ----------
    n : int
        total number of states
    jmp : float
        transition rate

    Returns
    -------
    dictionary:
        plast : la.lnarray
            potentiation/depression transition matrices
        weight : la.lnarray
            synaptic weights (linear/binary)
        signal : la.lnarray
            desired signal contribution from each plasticity type
    """
    out = build_empty(nst, 2, True)
    out['plast'] = serial_trans(nst, jmp)
    return out


def build_multistate(nst: int, jmp: float) -> Dict[str, la.lnarray]:
    """Make a multistate model.

    Make a multistate model, i.e. uniform nearest neighbour transition rates,
    and linear weights.

    For use with `SynapseMemoryModel.Build`

    Parameters
    ----------
    n : int
        total number of states
    jmp : float
        transition rate

    Returns
    -------
    dictionary:
        plast : la.lnarray
            potentiation/depression transition matrices
        weight : la.lnarray
            synaptic weights (linear/binary)
        signal : la.lnarray
            desired signal contribution from each plasticity type
    """
    out = build_empty(nst, 2, False)
    out['plast'] = serial_trans(nst, jmp)
    return out


def build_cascade(nst: int, jmp: float) -> Dict[str, la.lnarray]:
    """Make a cascade model.

    Make a cascade model with geometric transition rates and binary weights.
    For use with `SynapseMemoryModel.Build`.

    Parameters
    ----------
    n : int
        total number of states
    jmp : float
        transition rate parameter

    Returns
    -------
    dictionary:
        plast : la.lnarray
            potentiation/depression transition matrices
        weight : la.lnarray
            synaptic weights (linear/binary)
        signal : la.lnarray
            desired signal contribution from each plasticity type
    """
    out = build_zero(nst, 2, True)
    plast = out['plast']
    n = nst // 2
    jmp_vec = np.logspace(n-2, 0, num=n-1, base=jmp)
    jmp /= (1-jmp)

    inds = np.arange(1, n)
    plast[0, inds, n] = jmp_vec
    plast[0, 0, n] = jmp_vec[0] * jmp
    plast[0, inds, inds-1] = jmp_vec * jmp

    inds = nst - 1 - inds
    plast[1, inds, n-1] = jmp_vec
    plast[1, -1, n-1] = jmp_vec[0] * jmp
    plast[1, inds, inds+1] = jmp_vec * jmp

    stochastify_c(plast)
#    out['plast'] = plast
    return out


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
