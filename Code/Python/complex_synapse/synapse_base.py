# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 15:49:42 2017

@author: Subhy
"""
from __future__ import annotations
from typing import Union, Sequence, Dict
from numbers import Number
import numpy as np
import numpy_linalg as la
from . import builders as bld
cvl = la.convert

# types that can multiply with/add to a matrix
ArrayLike = Union[Number, Sequence[Number], np.ndarray]


class SynapseBase(np.lib.mixins.NDArrayOperatorsMixin):
    """Base class for complex synapses.

    Contains methods that modify instance variables, those that do not
    perform calculations and overloads of arithmetic operators, str and repr.
    Subclasses should override: `__init__`, `dict_copy` if necessary,
    including calls to super().

    Parameters (and attributes)
    ---------------------------
    plast : la.lnarray
        potentiation/depression transition rate matrix.
    frac : array_like
        fraction of events that are potentiating/depressing.
    weight : array_like
        synaptic weight of each state.
    signal : array_like
        desired signal contribution from each plasticity type.

    Properties
    ----------
    nstates
        number of states.
    nplast
        number of plasticity types.
    """
    # Attributes

    # potentiation/depression transition rate matrices.
    plast: la.lnarray
    # fraction of events that are potentiating./depressing
    frac: la.lnarray

    def __init__(self, plast: la.lnarray,
                 frac: ArrayLike = 0.5):
        """Class for complex synapse models.

        Parameters (and attributes)
        ---------------------------
        plast : la.lnarray
            potentiation/depression transition rate matrix.
        frac : array_like
            fraction of events that are potentiating/depressing.
        weight : array_like
            synaptic weight of each state.
        """
        # store inputs
        self.plast = la.asarray(plast)
        self.frac = la.asarray(frac).ravel()
        if len(self.frac) == len(self.plast) - 1:
            self.frac = la.concatenate((self.frac, [1. - self.frac.sum()]))

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Handling ufuncs with SynapseBases
        """
        args, _ = cvl.conv_in_attr('plast', SynapseBase, inputs)

        conv = [True] + [False] * (ufunc.nout-1)
        outputs, conv = cvl.conv_in_attr('plast', SynapseBase, kwargs, conv)

        results = self.plast.__array_ufunc__(ufunc, method, *args, **kwargs)

        return cvl.conv_out_attr(self, 'plast', results, outputs, conv)

    # -------------------------------------------------------------------------
    # Housekeeping
    # -------------------------------------------------------------------------

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = type(self).__name__ + "(\n"
        rpr += "    plast = "
        rpr += repr(self.plast).replace("\n", "\n" + " " * 12) + ",\n"
        rpr += f"    frac = {self.frac!r},\n)"
        return rpr

    def __str__(self) -> str:
        """Short representation of object"""
        return f"{type(self).__name__} with M={self.nstates}, fpm={self.frac}"

    # -------------------------------------------------------------------------
    # Utility methods
    # -------------------------------------------------------------------------

    def view(self, typ):
        """View of plasticity matrices as typ
        """
        return self.plast.view(typ)

    def dict_copy(self, keys=(), order='C', **kwds) -> Dict[str, la.lnarray]:
        """Dictionary with copies of data attributes
        """
        # In subclass:
        # keys += ('new_attr', ...)
        # return super().dict_copy(keys=keys, order=order, **kwds)
        keys += ('plast', 'frac')
        for k in keys:
            if k not in kwds.keys():
                kwds[k] = getattr(self, k).copy(order)
        return kwds

    def copy(self, *args, **kwargs) -> SynapseBase:
        """Copy of object, with copies of attributes
        """
        return type(self)(**self.dict_copy(*args, **kwargs))

    @property
    def nplast(self) -> int:
        """number of plasticity types."""
        return self.frac.size

    @property
    def nstates(self) -> int:
        """Number of states."""
        return self.plast.shape[-1]

    # -------------------------------------------------------------------------
    # Factory methods
    # -------------------------------------------------------------------------

    @classmethod
    def build(cls, builder, nst: int, frac: ArrayLike = 0.5,
              extra_args=(), **kwargs) -> 'SynapseBase':
        """Build model from function.

        Parameters
        ----------
        builder : function
            function object with parameters (n,...) that returns dictionary
            {plast, weight, etc}, ignoring any extra elements.
        nst: int
            total number of states, passed to `builder`.
        frac : float
            fraction of events that are potentiating, default=0.5.
        extra_args, **kwargs
            *extra_args, **kwargs: passed to `builder`.

        Returns
        -------
        synobj
            SynapseBase instance
        """
        return cls(**builder(nst, *extra_args, **kwargs), frac=frac)

    @classmethod
    def zero(cls, *args, **kwargs) -> SynapseBase:
        """Zero model

        Synapse model with all transition rates set to zero

        Parameters
        ----------
            All passed to `cls.build` or `builders.build_zero`
        nst: int
            total number of states
        npl: int
            total number of plasticity types
        frac : float
            fraction of events that are potentiating, default=0.5.
        binary : bool
            is the weight vector binary? Otherwise it's linear. Default: False
        ...
            extra arguments passed to `cls.build` or `builders.build_zero`.

        Returns
        -------
        synobj
            SynapseBase instance
        """
        return cls.build(bld.build_zero, *args, **kwargs)

    @classmethod
    def empty(cls, *args, **kwargs) -> SynapseBase:
        """Empty model

        Synapse model with all transition rates uninitialised.

        Parameters
        ----------
            All passed to `cls.build` or `builders.build_empty`
        nst: int
            total number of states
        npl: int
            total number of plasticity types
        frac : float
            fraction of events that are potentiating, default=0.5.
        binary : bool
            is the weight vector binary? Otherwise it's linear. Default: False
        ...
            extra arguments passed to `cls.build` or `builders.build_empty`.

        Returns
        -------
        synobj
            SynapseBase instance
        """
        return cls.build(bld.build_empty, *args, **kwargs)

    @classmethod
    def rand(cls, nst, *args, **kwargs) -> SynapseBase:
        """Random model

        Synapse model with random transition matrices

        Parameters
        ----------
            All passed to `cls.build`, `builders.build_rand`
            or `builders.rand_trans`:
        nst: int
            total number of states
        npl: int
            total number of plasticity types
        frac : float
            fraction of events that are potentiating, default=0.5.
        binary : bool
            is the weight vector binary? Otherwise it's linear. Default: False
        sp : float
            sparsity, default: 1.
        ...
            extra arguments passed to `cls.build`, `builders.build_rand`
            or `Builder.rand_trans`.

        Returns
        -------
        synobj
            SynapseBase instance
        """
        return cls.build(bld.build_rand, nst, *args, **kwargs)
