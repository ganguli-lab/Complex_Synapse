# -*- coding: utf-8 -*-
# Created on Mon Sep 18 15:49:42 2017
# @author: Subhy
"""Base class for complex synapse models, andd utility functions
"""
from __future__ import annotations

from numbers import Number
from typing import (Any, Callable, ClassVar, Dict, Optional, Sequence, Tuple,
                    Type, TypeVar, Union)

import numpy as np

import numpy_linalg as la
import numpy_linalg.convert as _cvl
import sl_py_tools.arg_tricks as _ag
import sl_py_tools.numpy_tricks.markov as _ma

from . import builders as _bld

# =============================================================================


class SynapseBase(np.lib.mixins.NDArrayOperatorsMixin):
    """Base class for complex synapses.

    Contains methods that modify instance variables, those that do not
    perform calculations and overloads of arithmetic operators, str and repr.
    Subclasses should override: `__init__`, `dict_copy` if necessary,
    including calls to super().

    Parameters (and attributes)
    ---------------------------
    plast : array_like, (P,M,M), float[0:1]
        potentiation/depression transition rate matrix.
    frac : array_like, (P,), float[0:1]
        fraction of events that are potentiating/depressing.

    Properties
    ----------
    nstate : int
        number of states, M.
    nplast : int
        number of plasticity types, P.
    nmodel : Tuple[int]
        Number and shape of models being broadcast.
    """
    # Attributes

    # potentiation/depression transition rate matrices.
    plast: la.lnarray
    # fraction of events that are potentiating./depressing
    frac: la.lnarray
    # largest row sum for valid plast & frac
    StochThresh: ClassVar[float] = 1e-5
    # largest condition number for inverting zinv
    CondThresh: ClassVar[float] = 1e-5

    def __init__(self, plast: ArrayLike,
                 frac: ArrayLike = 0.5):
        """Class for complex synapse models.

        Parameters (and attributes)
        ---------------------------
        plast : array_like, (P,M,M), float[0:1]
            potentiation/depression transition rate matrix.
        frac : array_like, (P,), float[0:1]
            fraction of events that are potentiating/depressing.
        """
        # store inputs
        self.plast = la.asarray(plast)
        self.frac = append_frac(frac, self.nplast)
        self.frac = trim_frac(self.frac, self.nplast)

    def __array_ufunc__(self, ufunc: np.ufunc, method: str, *inputs, **kwargs):
        """Handling ufuncs with SynapseBases
        """
        args, _ = _cvl.conv_in_attr('plast', SynapseBase, inputs)

        conv = [True] + [False] * (ufunc.nout-1)
        outputs, conv = _cvl.conv_in_attr('plast', SynapseBase, kwargs, conv)

        results = self.plast.__array_ufunc__(ufunc, method, *args, **kwargs)

        return _cvl.conv_out_attr(self, 'plast', results, outputs, conv)

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
        return f"{type(self).__name__} with M={self.nstate}, f={self.frac}"

    # -------------------------------------------------------------------------
    # Utility methods
    # -------------------------------------------------------------------------

    def view(self, astype: Optional[Type[Syn]], **kwds) -> Syn:
        """Copy of object, with views of array attributes

        Parameters
        ----------
        astype : type|None, optional
            The returned view is cast to this type. Requires the names of
            attributes/parameters to be the same in both classes.
            By default `None` -> `type(self)`.

        Requires `__init__` parameter names to be the same as attribute names.
        """
        astype = _ag.default(astype, type(self))
        attrs = array_attrs(self)
        attrs.update(kwds)
        return astype(**attrs)

    def copy(self, astype: Optional[Type[Syn]], order: str = 'C', **kwargs
             ) -> Syn:
        """Copy of object, with copies of array attributes

        Parameters
        ----------
        astype : type|None, optional
            The returned view is cast to this type. Requires the names of
            attributes/parameters to be the same in both classes.
            By default `None` -> `type(self)`.
        order : str, optional
            Memory order of new arrays. By default `'C'`.

        Requires `__init__` parameter names to be the same as attribute names.
        """
        astype = _ag.default(astype, type(self))
        attrs = array_attrs(self)
        for k in attrs:
            if k not in kwargs.keys():
                kwargs[k] = attrs[k].copy(order)
        return astype(**kwargs)

    @property
    def nplast(self) -> int:
        """number of plasticity types, P."""
        return self.plast.shape[-3]

    @property
    def nstate(self) -> int:
        """Number of states, M."""
        return self.plast.shape[-1]

    @property
    def nmodel(self) -> Tuple[int, ...]:
        """Number and shape of models being broadcast."""
        return self.plast.shape[:-3]

    # -------------------------------------------------------------------------
    # Factory methods
    # -------------------------------------------------------------------------

    @classmethod
    def build(cls, builder: Callable[..., Dict[str, la.lnarray]],
              nst: int, frac: ArrayLike = 0.5,
              extra_args=(), **kwargs) -> SynapseBase:
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
            extra arguments passed to `cls.build`, `builders.build_zero`.
            or `sl_py_tools.numpy_tricks.markov...`.

        Returns
        -------
        synobj
            SynapseBase instance
        """
        return cls.build(_bld.build_zero, *args, **kwargs)

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
            extra arguments passed to `cls.build`, `builders.build_empty`
            or `sl_py_tools.numpy_tricks.markov...`.

        Returns
        -------
        synobj
            SynapseBase instance
        """
        return cls.build(_bld.build_empty, *args, **kwargs)

    @classmethod
    def rand(cls, nst, *args, **kwargs) -> SynapseBase:
        """Random model

        Synapse model with random transition matrices

        Parameters
        ----------
            All passed to `cls.build`, `builders.build_rand`
            or `sl_py_tools.numpy_tricks.markov_param`:
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
            or `sl_py_tools.numpy_tricks.markov...`.

        Returns
        -------
        synobj
            SynapseBase instance
        """
        return cls.build(_bld.build_rand, nst, *args, **kwargs)


# =============================================================================
# Utility functions
# =============================================================================


def append_frac(frac: ArrayLike, npl: int) -> la.lnarray:
    """Add an extra element to the end of `frac` if it is one short or
    subnormalised

    Parameters
    ----------
    frac : ArrayLike
        Fraction of events of each plasticity type
    npl : int
        Number of plasticity types. If you only want to fix subnormalised
        `frac`, set this to 0.

    Returns
    -------
    frac : la.lnarray
        Input with extra elements, if needed
    """
    frac = np.atleast_1d(la.asarray(frac))
    stoch_thresh = 1 - SynapseBase.StochThresh
    if frac.shape[-1] == npl - 1 or (frac.sum(-1) < stoch_thresh).any():
        total = frac.sum(axis=-1, keepdims=True)
        frac = np.concatenate((frac, 1. - total), axis=-1)
    return frac


def trim_frac(frac: ArrayLike, npl: int) -> la.lnarray:
    """Remove elements from the end of `frac` if it has too many

    Parameters
    ----------
    frac : ArrayLike
        Fraction of events of each plasticity type
    npl : int
        Number of plasticity types.

    Returns
    -------
    frac : la.lnarray
        Input with extra elements, if needed
    """
    frac = np.atleast_1d(la.asarray(frac))
    frac = frac[..., :npl]
    _ma.stochastify_d(frac)
    return frac


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


def instance_attrs(obj: Any) -> Dict[str, Any]:
    """Dict of instance attributes that are not class attributes

    Parameters
    ----------
    obj : Any
        The instance whose attributes we want

    Returns
    -------
    attrs : Dict[str, Any]
        Instance attributes
    """
    names = set(dir(obj)) - set(dir(type(obj)))
    return {name: getattr(obj, name) for name in names}


def array_attrs(obj: Any) -> Dict[str, np.ndarray]:
    """Dict of instance array attributes

    Parameters
    ----------
    obj : Any
        The instance whose attributes we want

    Returns
    -------
    attrs : Dict[str, Any]
        Instance attributes
    """
    attrs = instance_attrs(obj)
    return {k: v for k, v in attrs.items() if isinstance(v, np.ndarray)}


# =============================================================================
# Hinting variables and aliases
# =============================================================================
# types that can multiply with/add to a matrix
ArrayLike = Union[Number, Sequence[Number], np.ndarray]
Syn = TypeVar('Syn', bound=SynapseBase)
