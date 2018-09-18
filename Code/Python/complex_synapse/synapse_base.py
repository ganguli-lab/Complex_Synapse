# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 15:49:42 2017

@author: Subhy
"""
from typing import ClassVar, Union, Sequence
from numbers import Number
import numpy as np
from . import builders as bld
from .builders import la

# types that can multiply with/add to a matrix
ArrayLike = Union[Number, Sequence[Number], np.ndarray]


class SynapseBase(np.lib.NDArrayOperatorsMixin):
    """Base class for complex synapses.

    Contains methods that modify instance variables, those that do not
    perform calculations and overloads of arithmetic operators, str and repr.
    Subclasses should override `init`, `copy`, `fix`, `valid_{shapes,values}`.

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
                 frac: ArrayLike = 0.5,
                 **kwargs):
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

    (__matmul__,
     __rmatmul__,
     __imatmul__) = np.lib.mixins._numeric_methods(la.matmul, 'matmul')

    def copy(self) -> 'SynapseBase':
        """Copy of object, with copies of attributes
        """
        return type(self)(plast=self.plast.copy(),
                          frac=self.frac.copy())

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Handling ufunce with SynapseBases
        """
        args = []
        # For most inputs, we swap multiplication & division instead of inverse
        for input_ in inputs:
            if isinstance(input_, SynapseBase):
                args.append(input_.plast)
            else:
                args.append(input_)
        args = tuple(args)

        conv_out = [False] * ufunc.nout
        conv_out[0] = True
        outputs = kwargs.pop('out', None)
        if outputs:
            out_args = []
            for j, output in enumerate(outputs):
                if isinstance(output, SynapseBase):
                    out_args.append(output.plast)
                    conv_out[j] = True
                else:
                    out_args.append(output)
                    conv_out[j] = False
            kwargs['out'] = tuple(out_args)
        else:
            outputs = (None,) * ufunc.nout

        results = self._to_invert.__array_ufunc__(ufunc, method, *args,
                                                  **kwargs)

        if results is NotImplemented:
            return NotImplemented

        if ufunc.nout == 1:
            results = (results,)

        converted_results = []
        for result, output, cout in zip(results, outputs, conv_out):
            if output is None:
                if cout:
                    conv = self.copy()
                    conv.plast = result
                    converted_results.append(conv)
                else:
                    converted_results.append(result)
            else:
                converted_results.append(output)
        results = tuple(converted_results)

        return results[0] if len(results) == 1 else results

    @property
    def nplast(self) -> int:
        """number of plasticity types."""
        return self.frac.size

    @property
    def nstates(self) -> int:
        """Number of states."""
        return self.plast.shape[-1]

    def normalise(self):
        """Ensure that all attributes are valid.
        """
        bld.stochastify_c(self.plast)
        bld.stochastify_d(self.frac)

    def fix(self):
        """Complete frac vector
        """
        if len(self.frac) == len(self.plast) - 1:
            self.frac = np.concatenate((self.frac, [1. - self.frac.sum()]))

    def valid_shapes(self) -> bool:
        """Do attributes (plast, weight, frac) have correct shapes?"""
        vld = self.plast.shape[-2] == self.plast.shape[-1]
        vld &= len(self.plast) == len(self.frac)
        return vld

    def valid_values(self) -> bool:
        """Do attributes (plast, frac) have valid values?"""
        vld = bld.isstochastic_c(self.plast, self.StochThresh)
        vld &= (self.frac >= 0.).all()
        vld &= bld.isstochastic_d(self.frac, self.StochThresh)
        return vld

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
    def zero(cls, *args, **kwargs) -> 'SynapseBase':
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
    def empty(cls, *args, **kwargs) -> 'SynapseBase':
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
    def rand(cls, *args, **kwargs) -> 'SynapseBase':
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

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = " Wpm = {}\n fpm = {}\n   w = {}"
        rpr = rpr.format(str(self.plast).replace("\n", "\n" + " " * 7),
                         self.frac, self.weight)
        return rpr

    def __str__(self) -> str:
        """Short representation of object"""
        rpr = "{} with {} states, fpm = {}, w = {}"
        rpr = rpr.format(type(self).__name__, self.nstates, self.frac,
                         self.weight)
        return rpr
