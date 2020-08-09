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

from typing import Dict, Union

import numpy as np

import numpy_linalg as la
from sl_py_tools.numpy_tricks import markov as ma
from sl_py_tools.numpy_tricks.markov import params as mp

RNG = la.random.default_rng()

# =============================================================================


def insert_axes(arr: np.ndarray, how_many: int, where: int = -1) -> np.ndarray:
    """Add multiple sigleton axes at once.

    New dimensions added to `arr` at:
    * `where, where+1, ... where+how_many-1` if `where` is non-negative.
    * `where, where-1, ... where-how_many+1` if `where` is negative.
    """
    sgn = -1 if where < 0 else 1
    axes = tuple(range(where, where + sgn * how_many, sgn))
    return np.expand_dims(arr, axes)


def scalarise(arg: np.ndarray) -> Union[np.ndarray, np.generic]:
    """Replace array with scalar if ndim==0."""
    if arg.ndim == 0:
        return arg[()]
    return arg


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
    return np.hstack((-weight, weight))


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


def serial_trans(nst: int, npl: int = 2, jmp: float = 1.) -> la.lnarray:
    """
    Make a serial transition matrix (continuous time).

    Parameters
    ----------
    n : int
        total number of states
    jmp : float
        jump rate

    Returns
    -------
    mat : la.lnarray
        transition matrix
    """
    if npl > 2:
        raise ValueError("serial mdel has 2 plasticity types")
    if npl == 1:
        # drn scalar, so no broadcast axis. drn = 0, so shape[-1] == 2.
        return mp.uni_serial_params_to_mat([jmp, jmp], nst, drn=0)
    # broadcast drn with axis 0. drn != 0, so shape[-1] == 1.
    return mp.uni_serial_params_to_mat([[jmp], [jmp]], nst, drn=(1, -1))


def cascade_trans(nst: int, npl: int = 2, jmp: float = 0.5) -> la.lnarray:
    """
    Make a cascade transition matrix (continuous time).

    Parameters
    ----------
    n : int
        total number of states
    jmp : float
        jump rate parameter of cascade model

    Returns
    -------
    mat : la.lnarray
        transition matrix
    """
    if npl > 2:
        raise ValueError("cascade mdel has 2 plasticity types")
    if npl == 1:
        # drn is scalar, so no broadcast axis. drn = 0, so shape[-1] == 2.
        return mp.std_cascade_params_to_mat([jmp, jmp], nst, drn=0)
    # broadcast drn with axis 0. drn != 0, so shape[-1] == 1.
    return mp.std_cascade_params_to_mat([[jmp], [jmp]], nst, drn=(1, -1))


def build_generic(func, nst: int, npl: int = 2, binary: bool = False,
                  **kwds) -> Dict[str, la.lnarray]:
    """Make a model from a matrix creating function.

    Factory for model builderss.
    For use with `SynapseMemoryModel.build`.

    Parameters
    ----------
    func : function
        function with signature `func(nst: int, npl: int, **kwds) -> lnarray`
        and `func(nst, npl, **kwds).shape = (npl, nst, nst)`
    nst : int
        total number of states
    npl : optional, int = 2
        number of plasticity matrices
    binary : optional, bool = True
        is the weight vector binary? Otherwise it's linear. Default: False
    Other keywords passed to `func`.

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
    plast = la.asarray(func(nst, npl, **kwds))
    if binary:
        weight = binary_weights(nst)
    else:
        weight = linear_weights(nst)
    signal = la.linspace(1, -1, npl)
    return {'plast': plast, 'weight': weight, 'signal': signal}


def build_standard(func, nst: int, npl: int = 2, binary: bool = False,
                   **kwds) -> Dict[str, la.lnarray]:
    """Make a model from a standard array creating function.

    Factory for model builderss.
    For use with `SynapseMemoryModel.build`.

    Parameters
    ----------
    func : function
        function with signature `func(shape: Tuple[int], **kwds) -> lnarray`
        and `func(shape, **kwds).shape = shape
    nst : int
        total number of states
    npl : optional, int = 2
        number of plasticity matrices
    binary : optional, bool = True
        is the weight vector binary? Otherwise it's linear. Default: False
    Other keywords passed to `func`.

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
    plast = la.asarray(func((npl, nst, nst), **kwds))
    if binary:
        weight = binary_weights(nst)
    else:
        weight = linear_weights(nst)
    signal = la.linspace(1, -1, npl)
    return {'plast': plast, 'weight': weight, 'signal': signal}


def build_zero(nst: int, npl: int = 2, binary: bool = False
               ) -> Dict[str, la.lnarray]:
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
    return build_standard(la.zeros, nst, npl, binary)


def build_empty(nst: int, npl: int = 2, binary: bool = False
                ) -> Dict[str, la.lnarray]:
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
    return build_standard(la.empty, nst, npl, binary)


def build_rand(nst: int, npl: int = 2, binary: bool = False,
               **kwds) -> Dict[str, la.lnarray]:
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
        passed to `sl_py_tools.numpy_tricks.markov.rand_trans`.
    sparsity : float
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
    return build_generic(ma.rand_trans, nst, npl, binary, **kwds)


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
    return build_generic(serial_trans, nst, 2, True, jmp=jmp)


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
    return build_generic(serial_trans, nst, 2, False, jmp=jmp)


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
    return build_generic(cascade_trans, nst, 2, True, jmp=jmp)
