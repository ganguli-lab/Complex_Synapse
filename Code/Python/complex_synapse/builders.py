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


def binary_weights(nst: int) -> np.ndarray:  # binary synaptic weights
    """
    Make binary (+/-1) weight vector

    Parameters
    ----------
    n : int
        total number of states

    Returns
    -------
    w : np.ndarray
        vector of synaptic weights
    """
    weight = np.ones(nst // 2)
    return np.hstack((-weight, weight))


def linear_weights(nst: int) -> np.ndarray:  # linear synaptic weights
    """
    Make linear (-1,...,1) weight vector

    Parameters
    ----------
    n : int
        total number of states

    Returns
    -------
    w : np.ndarray
        vector of synaptic weights
    """
    return np.linspace(-1., 1, nst)


def stochastify_c(mat: np.ndarray):  # make cts time stochastic
    """
    Make a matrix the generator of a continuous time Markov process.
    Changes diagonal to make row sums zero.
    **Modifies** in place, **does not** return.

    Parameters
    ----------
    mat : np.ndarray
        square matrix with non-negative off-diagonal elements.
        **Modified** in place.
    """
    mat -= np.apply_along_axis(np.diag, -1, mat.sum(axis=-1))


def stochastify_d(mat: np.ndarray):  # make dsc time stochastic
    """
    Make a matrix the generator of a discrete time Markov process.
    Scales rows to make row sums one.

    Parameters
    ----------
    mat : np.ndarray
        square matrix with non-negative elements.
        **Modified** in place
    """
    mat /= mat.sum(axis=-1, keepdims=True)


def isstochastic_c(mat: np.ndarray, thresh: float = 1e-5) -> bool:
    """Are row sums zero?
    """
    return (np.fabs(mat.sum(axis=-1)) < thresh).all()


def isstochastic_d(mat: np.ndarray, thresh: float = 1e-5) -> bool:
    """Are row sums one?
    """
    return (np.fabs(mat.sum(axis=-1) - 1.) < thresh).all()


def rand_trans(nst: int, sparsity: float = 1.) -> np.ndarray:
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
    mat : np.ndarray
        transition matrix
    """
    mat = np.random.rand(nst, nst)
    ind = np.random.rand(nst, nst)
    mat[ind > sparsity] = 0.
    stochastify_c(mat)
    return mat


def build_generic(nst: int, func, npl: int = 2,
                  binary: bool = False) -> Dict[str, np.ndarray]:
    """Make a model from a matrix creating function.

    Factory for model builderss.
    For use with `SynapseMemoryModel.Build`.

    Parameters
    ----------
    n : int
        total number of states
    Func : function
        function with signature `func(nst: int) -> np.ndarray` and
        func(nst).shape = (nst, nst)
    npl : optional, int = 2
        number of plasticity matrices
    binary : optional, bool = True
        is the weight vector binary? Otherwise it's linear. Default: False

    Returns
    -------
    dictionary:
        plast : np.ndarray
            potentiation/depression transition matrices
        weight : np.ndarray
            synaptic weights (linear/binary)
        signal : np.ndarray
            desired signal contribution from each plasticity type
    """
    plast = [func(nst) for __ in range(npl)]
    if binary:
        weight = binary_weights(nst)
    else:
        weight = linear_weights(nst)
    signal = np.linspace(1, -1, npl)
    return {'plast': np.array(plast), 'weight': weight, 'signal': signal}


def build_zero(nst: int, npl: int = 2,
               binary: bool = False) -> Dict[str, np.ndarray]:
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
        plast : np.ndarray
            potentiation/depression transition matrices
        weight : np.ndarray
            synaptic weights (linear/binary)
        signal : np.ndarray
            desired signal contribution from each plasticity type
    """
    return build_generic(nst, lambda n: np.zeros((n, n)), npl, binary)


def build_empty(nst: int, npl: int = 2,
                binary: bool = False) -> Dict[str, np.ndarray]:
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
        plast : np.ndarray
            potentiation/depression transition matrices
        weight : np.ndarray
            synaptic weights (linear/binary)
        signal : np.ndarray
            desired signal contribution from each plasticity type
    """
    return build_generic(nst, lambda n: np.empty((n, n)), npl, binary)


def build_rand(nst: int, npl: int = 2, binary: bool = False,
               sparsity: float = 1.) -> Dict[str, np.ndarray]:
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
        plast : np.ndarray
            potentiation/depression transition matrices
        weight : np.ndarray
            synaptic weights (linear/binary)
        signal : np.ndarray
            desired signal contribution from each plasticity type
    """
    return build_generic(nst, lambda n: rand_trans(nst, sparsity), npl, binary)


def build_serial(nst: int, jmp: float) -> Dict[str, np.ndarray]:
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
        plast : np.ndarray
            potentiation/depression transition matrices
        weight : np.ndarray
            synaptic weights (linear/binary)
        signal : np.ndarray
            desired signal contribution from each plasticity type
    """
    out = build_empty(nst, 2, True)
    out['plast'] = np.stack((np.diag(jmp * np.ones(nst-1), 1),
                             np.diag(jmp * np.ones(nst-1), -1)))
    stochastify_c(out['plast'])
    return out


def build_multistate(nst: int, jmp: float) -> Dict[str, np.ndarray]:
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
        plast : np.ndarray
            potentiation/depression transition matrices
        weight : np.ndarray
            synaptic weights (linear/binary)
        signal : np.ndarray
            desired signal contribution from each plasticity type
    """
    out = build_empty(nst, 2, False)
    out['plast'] = np.stack((np.diag(jmp * np.ones(nst-1), 1),
                             np.diag(jmp * np.ones(nst-1), -1)))
    stochastify_c(out['plast'])
    return out


def build_cascade(nst: int, jmp: float) -> Dict[str, np.ndarray]:
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
        plast : np.ndarray
            potentiation/depression transition matrices
        weight : np.ndarray
            synaptic weights (linear/binary)
        signal : np.ndarray
            desired signal contribution from each plasticity type
    """
    out = build_zero(nst, 2, True)
    plast = out['plast']
    n = nst // 2
    jmp_vec = np.logspace(n-2, 0, num=n-1, base=jmp)
    jmps = jmp / (1-jmp)

    inds = np.arange(1, n)
    plast[0, inds, n] = jmp_vec
    plast[0, 0, n] = jmp_vec[0] * jmps
    plast[0, inds, inds-1] = jmp_vec * jmps

    inds = nst - 1 - inds
    plast[1, inds, n-1] = jmp_vec
    plast[1, -1, n-1] = jmp_vec[0] * jmps
    plast[1, inds, inds+1] = jmp_vec * jmps

    stochastify_c(plast)
    out['plast'] = plast
    return out
