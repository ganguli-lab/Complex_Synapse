# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 15:49:42 2017

@author: Subhy
"""
from __future__ import annotations
from typing import ClassVar, Union, Sequence, Dict
from numbers import Number
import numpy as np
from . import builders as bld
from .builders import la, ma
cvl = la.convert_loop

# types that can multiply with/add to a matrix
ArrayLike = Union[Number, Sequence[Number], np.ndarray]


class SynapseBase(la.gufuncs.LNArrayOperatorsMixin):
    """Base class for complex synapses.

    Contains methods that modify instance variables, those that do not
    perform calculations and overloads of arithmetic operators, str and repr.
    Subclasses should override: `__init__`, `dict_copy`, `fix`, `valid_shapes`,
    `valid_values}`, `normalise` if necessary, including calls to super().
    Beware of `fix`: it is called in super().__init__, so any reference to
    subclass attributes should be enclosed by `try ... except AttributeError`.

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

    # Common constatnts / parameters

    # largest row sum for valid plast & frac
    StochThresh: ClassVar[float] = 1e-5

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
        self.fix()

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Handling ufunce with SynapseBases
        """
        args = cvl.conv_loop_in_attr('plast', SynapseBase, inputs)[0]

        conv = [False] * ufunc.nout
        conv[0] = True
        outputs = kwargs.pop('out', None)
        if outputs:
            out_args, conv = cvl.conv_loop_in_attr('plast', SynapseBase,
                                                   outputs)
            kwargs['out'] = tuple(out_args)
        else:
            outputs = (None,) * ufunc.nout

        results = self.plast.__array_ufunc__(ufunc, method, *args, **kwargs)

        if results is NotImplemented:
            return NotImplemented

        if ufunc.nout == 1:
            results = (results,)

        results = cvl.conv_loop_out_attr(self, 'plast', results, outputs, conv)

        return results[0] if len(results) == 1 else results

    # -------------------------------------------------------------------------
    # %%* Housekeeping
    # -------------------------------------------------------------------------

    def normalise(self):
        """Ensure that all attributes are valid.
        """
        ma.stochastify_c(self.plast)
        scale = -np.diagonal(self.plast).min()
        if scale > 1:
            self.plast /= scale
        ma.stochastify_d(self.frac)

    def fix(self):
        """Complete frac vector
        """
        if len(self.frac) == len(self.plast) - 1:
            self.frac = la.concatenate((self.frac, [1. - self.frac.sum()]))

    def valid_shapes(self) -> bool:
        """Do attributes (plast, weight, frac) have correct shapes?"""
        vld = self.plast.shape[-2] == self.plast.shape[-1]
        vld &= len(self.plast) == len(self.frac)
        return vld

    def valid_values(self) -> bool:
        """Do attributes (plast, frac) have valid values?"""
        vld = ma.isstochastic_c(self.plast, self.StochThresh)
        vld &= ma.isstochastic_d(self.frac, self.StochThresh)
        return vld

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = type(self).__name__ + "(\n"
        rpr += "    plast = "
        rpr += repr(self.plast).replace("\n", "\n" + " " * 12) + ",\n"
        rpr += f"    frac = {self.frac!r},\n)"
        return rpr

    def __str__(self) -> str:
        """Short representation of object"""
        return f"{type(self).__name_} with M={self.nstates}, fpm={self.frac}"

    # -------------------------------------------------------------------------
    # %%* Utility methods
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
    # %%* Factory methods
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
    def rand(cls, *args, **kwargs) -> SynapseBase:
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
        return cls.build(bld.build_rand, *args, **kwargs)
