# -*- coding: utf-8 -*-
# Created on Mon Sep 18 15:49:42 2017
# @author: Subhy
"""Base class for complex synapse models, andd utility functions
"""
from __future__ import annotations
from abc import abstractmethod, ABC

from numbers import Number
from typing import (Any, Callable, ClassVar, Dict, Optional, Sequence, Tuple,
                    Type, TypeVar, Union)

import numpy as np

import numpy_linalg as la
import numpy_linalg.convert as _cvl
import sl_py_tools.arg_tricks as _ag
import sl_py_tools.numpy_tricks.markov as _ma
import sl_py_tools.numpy_tricks.markov.params as _mp
from sl_py_tools.graph_tricks import MultiDiGraph, param_to_graph

import complex_synapse.builders as _bld

# =============================================================================


class SynapseBase(np.lib.mixins.NDArrayOperatorsMixin, ABC):
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

    def reorder(self, inds: ArrayLike) -> None:
        """Put the states into a new order, in-place.

        Parameters
        ----------
        inds : ArrayLike[int] (M,)
            `inds[i] = j` means state `j` moves to position `i`.
        """
        self.plast = self.plast[(...,) + np.ix_(inds, inds)]

    def view(self, astype: Optional[Type[Syn]] = None, **kwds) -> Syn:
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

    def copy(self, astype: Optional[Type[Syn]] = None, order: str = 'C',
             **kwargs) -> Syn:
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
        """Number of plasticity types, P."""
        return self.plast.shape[-3]

    @property
    def nstate(self) -> int:
        """Number of states, M."""
        return self.plast.shape[-1]

    @property
    def nmodel(self) -> Tuple[int, ...]:
        """Number and shape of models being broadcast."""
        return self.plast.shape[:-3]

    def to_graph(self, graph: Optional[MultiDiGraph] = None, **kwds
                 ) -> MultiDiGraph:
        """Create/modify a directed graph from a parameterised synapse model.

        Parameters
        ----------
        graph : MultiDiGraph|None
            Graph to be modified, by default `None` -> create a new graph.
        node_v : ndarray (M,)
            Values that determine node size (area), by default `self.peq()`.
        edge_v : ndarray (P,Q)
            Values that determine edge size (width), by default off-diagonals.
            `Q in {PM(M-1), P(M-1), P}`.
        node_k : ndarray (M,)
            Values that determine node colour, by default `np.zeros(M)`.
        edge_k : ndarray (P,)
            Values that determine edge colour, by default `topology.directions`
            or `[P:0:-1]`.
        topology : TopologyOptions
            Options describing topology of model class.

        Returns
        -------
        graph : MultiDiGraph
            Graph describing model.
        """
        node_v = _ag.default_eval(kwds.pop('node_v', None), self.peq)
        edge_v = kwds.pop('edge_v', None)
        if edge_v is None:
            edge_v = _ma.params.gen_mat_to_params(self.plast,
                                                  (0,) * self.nplast)
        if graph is None:
            node_k = kwds.pop('node_k', None)
            edge_k = kwds.pop('edge_k', None)
            topology = kwds.pop('topology', _ma.TopologyOptions())
            return param_to_graph(edge_v, node_v, node_k, edge_k, topology)
        graph.set_edge_attr('value', edge_v.ravel())
        graph.set_node_attr('value', node_v.ravel())
        return graph

    # -------------------------------------------------------------------------
    # Markov quantities
    # -------------------------------------------------------------------------

    def markov(self) -> la.lnarray:
        """Average transition rate matrix.

        .. math:: W^f = fp Wp + fm Wm
        """
        return (self.frac.s * self.plast).sum(-3)

    @abstractmethod
    def zinv(self, rate: Optional[ArrayLike] = None,
             rowv: Optional[ArrayLike] = None) -> la.lnarray:
        """Inverse of fundamental matrix"""

    def cond(self, rate: Optional[ArrayLike] = None, *,
             rowv: Optional[ArrayLike] = None,
             order: Order = None,
             rate_max: bool = False) -> la.lnarray:
        r"""Condition number of generalised fundamental matrix.

        Parameters
        ----------
        rate : float, array, optional
            Parameter of Laplace transform, ``s``. Default: 0.
        rowv : la.lnarray, optional
            Arbitrary row vector, ``xi``. If `None`, use vector of ones.
            If `np.nan`, use  `peq`. By default: `None`.
        order : {None, 1, -1, 2, -2, inf, -inf, 'fro'}, optional
            Order of the norm. By default `None`.
        rate_max : bool, optional
            Compute `max(cond(rate), cond(None))`. By default, `False`.

        Returns
        -------
        r : la.lnarray
             condition number of ``Z``.

        Notes
        -----
        ``c`` is the condition number, :math:`c(Z)`

        .. math:: c(Z) = \Vert Z \Vert \cdot \Vert Z^{-1} \Vert,

        where :math:`Z` is the inverse of the matrix returned by `self.zinv()`.

        See Also
        --------
        zinv : fundamental matrix
        """
        zinv = self.zinv(rate, rowv)
        cond = np.linalg.cond(zinv, order)
        if rate_max and rate is not None:
            zinv = self.zinv(None, rowv)
            cond = max(cond, np.linalg.cond(zinv, order))
        return cond

    def peq(self) -> la.lnarray:
        """Steady state distribution.

        Returns
        -------
        peq : la.lnarray
            Steady state distribution, :math:`\\pi`.
            Solution of: :math:`\\pi W^f = 0`.
        """
        rowv = la.ones(self.nstate)
        fundi = self.zinv(rowv=rowv)
        return rowv @ fundi.inv

    def time_rev(self) -> SynapseBase:
        """Swap plast with time rev"""
        plast = _ma.adjoint(self.plast, self.peq())
        _ma.stochastify_c(plast)
        return self.copy(plast=np.flip(plast, -3))

    # -------------------------------------------------------------------------
    # Factory methods
    # -------------------------------------------------------------------------

    @classmethod
    def build(cls: Type[Syn], builder: Callable[..., Dict[str, la.lnarray]],
              nst: int, frac: ArrayLike = 0.5,
              extra_args=(), **kwargs) -> Syn:
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
    def zero(cls: Type[Syn], *args, **kwargs) -> Syn:
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
    def empty(cls: Type[Syn], *args, **kwargs) -> Syn:
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
    def rand(cls: Type[Syn], nst, *args, **kwargs) -> Syn:
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


class SynapseContinuousTime(SynapseBase):
    """Class for Continuous time Markovian dynamics of complex synapse models.

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

    def zinv(self, rate: Optional[ArrayLike] = None,
             rowv: Optional[ArrayLike] = None) -> la.lnarray:
        r"""Inverse of generalised fundamental matrix.

        Parameters
        ----------
        rate : float, array, optional
            Parameter of Laplace transform, ``s``. Default: 0.
        rowv : la.lnarray, optional
            Arbitrary row vector, ``xi``. If `None`, use vector of ones.
            If `np.nan`, use  `peq`. By default: `None`.

        Returns
        -------
        Zi : la.lnarray
            inverse of generalised fundamental matrix, ``Z``.

        Notes
        -----
        ``zinv`` is the inverse of :math:`Z`, defined as

        .. math:: Z = (e \xi - W^f)^{-1},

        where :math:`e` is a vector of ones and :math:`\xi` is any row vector
        with :math:`\xi.e \neq 0`.
        When we include :math:`s`

        .. math:: Z(s) = (s I + e \xi - W^f)^{-1} \simeq \int e^{t(W^f-sI)} dt,

        i.e. :math:`Z^{-1}` with :math:`s` added to the diagonal.
        Effectively the matrix inverse of the genarator's Laplace transform.
        """
        onev = la.ones(self.nstate)
        if rowv is None:
            rowv = onev
        elif np.isnan(rowv).all():
            rowv = self.peq()
        else:
            rowv = la.asarray(rowv)
        if rate is None:
            s_arr = 0
        else:
            # convert to lnarray, add singletons to broadcast with matrix
            s_arr = la.asarray(rate).s * la.eye(self.nstate)
        return onev.c * rowv - self.markov() + s_arr


class SynapseDiscreteTime(SynapseBase):
    """Synapse model in discrete time

    Parameters (and attributes)
    ---------------------------
    plast : array_like, (P,M,M), float[0:1]
        potentiation/depression transition probability matrix.
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

    def zinv(self, rate: Optional[ArrayLike] = None,
             rowv: Optional[ArrayLike] = None) -> la.lnarray:
        r"""Inverse of generalised fundamental matrix.

        Parameters
        ----------
        rate : float, array, optional
            Parameter of discrete Laplace transform, ``s``. Default: 0.
        rowv : la.lnarray, optional
            Arbitrary row vector, ``xi``. If `None`, use vector of ones.
            If `np.nan`, use  `peq`. By default: `None`.

        Returns
        -------
        Zi : la.lnarray
            inverse of generalised fundamental matrix, ``Z``.

        Notes
        -----
        ``zinv`` is the inverse of :math:`Z`, defined as

        .. math:: Z = (e \xi + I - M^f)^{-1},

        where :math:`e` is a vector of ones, :math:`\xi` is any row vector
        with :math:`\xi.e \neq 0` and :math:`I` is the identity matrix.
        When we include :math:`s`

        .. math:: Z(s) = (e \xi + I - e^{-s} M^f)^{-1}
                        \simeq \sum_n e^{-ns} (M^f)^n,

        Effectively the genarator's discrete Laplace transform.
        """
        onev = la.ones(self.nstate)
        if rowv is None:
            rowv = onev
        elif np.isnan(rowv).all():
            rowv = self.peq()
        else:
            rowv = la.asarray(rowv)
        if rate is None:
            s_arr = 1
        else:
            # convert to lnarray, add singletons to broadcast with matrix
            s_arr = np.exp(-la.asarray(rate)).s
        return onev.c * rowv + la.eye(self.nstate) - self.markov() * s_arr


# =============================================================================
# Class for parameterised synapses
# =============================================================================


class SynapseParam(SynapseBase):
    """Class for complex synapses, for working with independent parameters

    Subclass of SynapseModel.
    Same constructor & attributes, only methods added.
    """
    topology: _ma.TopologyOptions

    def __init__(self, *args, **kwds) -> None:
        self.topology = kwds.pop('topology', _ma.TopologyOptions())
        self.topology.pop_my_args(kwds)
        super().__init__(*args, **kwds)

    @property
    def nparam(self) -> int:
        """number of plasticity types."""
        return self.nplast * _mp.num_param(self.nstate, **self.directed(0))

    def get_params(self, ravel: bool = True) -> la.lnarray:
        """Independent parameters of transition matrices.

        Returns
        -------
        params : la.lnarray
            Vector of off-diagonal elements, potentiation before depression,
            each in order:
            Wp_01, Wp_02, ... Wp_0n-1, Wp_10, Wp_12, ... Wp_n-2,n-1,
            Wm_01, Wm_02, ... Wm_0n-1, Wm_10, Wm_12, ... Wm_n-2,n-1.
            Unless a special topology is specified, in which case see
            `sl_py_tools.numpy_tricks.markov.params`.

        See Also
        --------
        markov.mat_to_params
        """
        params = _mp.mat_to_params(self.plast, **self.directed(daxis=-3))
        return params.ravelaxes(-2) if ravel else params

    def set_params(self, params: np.ndarray, *args, **kwds) -> None:
        """Update transition matrices from independent parameters.

        Does not update if parameters are unchanged Optional arguments are
        passed on to `numpy.allclose`.

        Parameters
        ----------
        params : np.ndarray
            Vector of off-diagonal elements, potentiation before depression,
            each in order:
            Wp_01, Wp_02, ... Wp_0n-1, Wp_10, Wp_12, ... Wp_n-2n-1,
            Wm_01, Wm_02, ... Wm_0n-1, Wm_10, Wm_12, ... Wm_n-2n-1.
            Unless a special topology is specified, in which case see
            `sl_py_tools.numpy_tricks.markov.params`.

        See Also
        --------
        markov.params_to_mat
        """
        if np.allclose(params, self.get_params(), *args, **kwds):
            return
        params = la.asarray(params).unravelaxis(-1, (self.nplast, -1))
        _mp.mat_update_params(self.plast, params, **self.directed(pdaxis=-2))

    def to_graph(self, graph: Optional[MultiDiGraph] = None, **kwds
                 ) -> MultiDiGraph:
        """Create/modify a directed graph from a parameterised synapse model.

        Parameters
        ----------
        graph : MultiDiGraph
            Graph to be modified.
        node_v : ndarray (M,)
            Values that determine node size (area), by default `self.peq()`.
        node_k : ndarray (M,)
            Values that determine node colour, by default `np.zeros(M)`.
        edge_k : ndarray (P,)
            Values that determine edge colour, by default `topology.directions`
            or `[P:0:-1]`.

        Returns
        -------
        graph : MultiDiGraph
            Graph describing model.
        """
        kwds.setdefault('edge_v', self.get_params(ravel=False))
        kwds.setdefault('topology', self.topology)
        return super().to_graph(graph, **kwds)

    def directed(self, which: Union[int, slice] = slice(None),
                 **kwds) -> Dict[str, Any]:
        """Markov parameter options

        Parameters
        ----------
        which : int, slice, optional
            Which element of `self.topology.directions` to use as the `'drn'`
            item, where `None` -> omit `'drn'` item. By default `slice(None)`.
        Extra arguments are default values or unknown keys in `opts`

        Returns
        -------
        opts : Dict[str, Any]
            All options for `sl_py_tools.numpy_tricks.markov.params`.
        """
        return self.topology.directed(which, **kwds)

    @classmethod
    def from_params(cls: Type[Syn], params: np.ndarray, *args, **kwds) -> Syn:
        """Builds SynapseOpt object from independent parameters

        Parameters
        ----------
        params : np.ndarray
            Vector of off-diagonal elements, potentiation before depression,
            each in order:
            mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2n-1.
        frac : float
            Fraction of events of each plasticity type, An extra element is
            added if it is subnormalised. By default 0.5.
        binary : bool
            is the weight vector binary? Otherwise it's linear. Default: False
        topology : TopologyOptions, optional keyword
            Topology specifying options. By default `TopologyOptions()`.
        npl : int, optional keyword
            Total number of plasticity types. By default `topology.npl`
        nst : int or None, optional keyword
            Total number of states. Calculate from `len(params)` if `None`.
        npar : int or None, optional keyword
            Total number of parameters. Ignored.
        ...
            extra arguments passed to `cls.zero`.

        Returns
        -------
        synobj
            SynapseParamModel instance
        """
        kwds.pop('rng', None)
        topology = kwds.pop('topology', _ma.TopologyOptions())
        topology.pop_my_args(kwds)
        npl = kwds.setdefault('npl', topology.npl)
        nst = kwds.pop('nst', None)
        kwds.pop('npar', None)
        if nst is None:
            nst = _mp.num_state(params.shape[-1]//npl, **topology.directed(0))
        obj = cls.zero(nst, *args, **kwds)
        obj.topology = topology
        obj.set_params(params)
        return obj

    @classmethod
    def rand(cls: Type[Syn], nst: Optional[int], *args, **kwds) -> Syn:
        """Random model

        Synapse model with random transition matrices

        Parameters
        ----------
            All passed to `cls.build`, `builders.build_rand`
            or `sl_py_tools.numpy_tricks.markov_param`:
        nst : int or None
            Total number of states. Calculate from `npar` if `None`.
        npar : int or None, optional keyword
            Total number of parameters. By default calculate from `nst`,
            (or if `None`).
        npl : int, optional keyword
            Total number of plasticity types. By default `topology.npl`.
        frac : float, optional keyword
            fraction of events that are potentiating, by default `0.5`.
        binary : bool, optional keyword
            is the weight vector binary? Otherwise it's linear. Default: False
        topology : TopologyOptions, optional keyword
            Topology specifying options. By default `TopologyOptions()`.
        rng : np.random.Generator, optional keyword
            Source of random numbers. By default, `builders.RNG`.
        ...
            extra arguments passed to `cls.from_params`.

        One of 'nst/npar' must be provided

        Returns
        -------
        synobj
            SynapseParamModel instance
        """
        rng = kwds.pop('rng', _bld.RNG)
        topology = kwds.pop('topology', _ma.TopologyOptions())
        topology.pop_my_args(kwds)
        npl = kwds.setdefault('npl', topology.npl)
        npar = kwds.pop('npar', None)
        if npar is None:
            if nst is None:
                raise TypeError("Must specify one of [nst, npar]")
            npar = npl * _mp.num_param(nst, **topology.directed(0))
        elif nst is None:
            nst = _mp.num_state(npar // npl, **topology.directed(0))
        kwds.update(nst=nst, topology=topology)
        params = rng.random(npar)
        if not topology.constrained:
            cnstr = (cls.constraint_coeff(npl, nst) @ params).max()
            if cnstr > 1:
                params /= cnstr
        return cls.from_params(params, *args, **kwds)

    @staticmethod
    def constraint_coeff(npl: int, nst: int) -> la.lnarray:
        """Coefficient matrix for upper bound on off-diagonal row sums.

        Parameters
        ----------
        npl : int
            Number of types of plasticity
        nst: int
            Number of states.

        Returns
        -------
        coeffs : la.lnarray (2*nst, 2*nst(nst-1))
            matrix of coefficients s.t ``coeffs @ params <= 1``.
        """
        rows, cols = la.arange(npl*nst), la.arange(nst-1)
        cols = cols + rows.c * (nst-1)
        coeffs = la.zeros((npl*nst, npl*nst*(nst-1)))
        coeffs[rows.c, cols] = 1
        return coeffs


# =============================================================================
# Utility functions
# =============================================================================


def append_frac(frac: ArrayLike, npl: int = 0) -> la.lnarray:
    """Add an element to the end of `frac` if it is one short or subnormalised

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
    total = frac.sum(axis=-1, keepdims=True)
    if frac.shape[-1] == npl - 1 or (total < stoch_thresh).any():
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


def constraint_coeff(npl: int, nst: int) -> la.lnarray:
    """Coefficient matrix for upper bound on off-diagonal row sums.

    Parameters
    ----------
    npl : int
        Number of types of plasticity
    nst: int
        Number of states.

    Returns
    -------
    coeffs : la.lnarray (2*nst, 2*nst(nst-1))
        matrix of coefficients s.t ``coeffs @ params <= 1``.
    """
    rows, cols = la.arange(npl*nst), la.arange(nst-1)
    cols = cols + rows.c * (nst-1)
    coeffs = la.zeros((npl*nst, npl*nst*(nst-1)))
    coeffs[rows.c, cols] = 1
    return coeffs


# =============================================================================
# Hinting variables and aliases
# =============================================================================
# types that can multiply with/add to a matrix
ArrayLike = Union[Number, Sequence[Number], np.ndarray]
Syn = TypeVar('Syn', bound=SynapseBase)
Order = Union[int, float, str, None]
Getter = Callable[[SynapseBase], np.ndarray]
