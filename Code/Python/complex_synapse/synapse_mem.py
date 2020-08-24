# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 18:22:05 2017

@author: Subhy
"""
from contextlib import contextmanager
from typing import ClassVar, Optional, Tuple, Union

import numpy as np

import numpy_linalg as la
import sl_py_tools.arg_tricks as _ag
import sl_py_tools.numpy_tricks.markov as _ma

from . import builders as _bld
from . import synapse_base as _sb

Order = Union[int, float, str, None]
# =============================================================================


class SynapseMemoryModel(_sb.SynapseBase):
    """Class for memory curves of complex synapse models.

    Parameters (and attributes)
    ---------------------------
    plast : array_like, (P,M,M), float[0:1]
        potentiation/depression transition rate matrix.
    frac : array_like, (P,), float[0:1]
        fraction of events that are potentiating/depressing.
    weight : array_like, (M,), float[-1:1]
        synaptic weight of each state.
    signal : array_like, (P,), float[-1:1]
        desired signal contribution from each plasticity type.

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

    # synaptic weight of each state.
    weight: la.lnarray
    # desired signal contribution from each plasticity type
    signal: la.lnarray

    # Common constatnts / parameters

    # degeneracy threshold, for evals or eta^+
    DegThresh: ClassVar[float] = 1e-3
    # # largest row sum for valid plast & frac
    # StochThresh: ClassVar[float] = 1e-5
    # # largest condition number for inverting zinv
    # CondThresh: ClassVar[float] = 1e-5
    # # smallest singular value for Split models
    # SingValThresh: ClassVar[float] = 1e-10
    # # threshold for lumpability test
    # LumpThresh: ClassVar[float] = 1e-3
    # # threshold for orphaned states
    # OrphanThresh: ClassVar[float] = 1e-3

    def __init__(self, plast: la.lnarray,
                 frac: _sb.ArrayLike = 0.5,
                 weight: Optional[_sb.ArrayLike] = None,
                 signal: Optional[_sb.ArrayLike] = None):
        """Class for complex synapse models.

        Parameters (and attributes)
        ---------------------------
        plast : la.lnarray
            potentiation/depression transition rate matrix.
        frac : array_like
            fraction of events that are potentiating.
        weight : array_like
            synaptic weight of each state.
        signal : array_like
            desired signal contribution from each plasticity type.
        """
        super().__init__(plast, frac)
        # store inputs
        self.weight = _ag.default_non_eval(weight, la.asarray,
                                           _bld.linear_weights(self.nstate))
        self.signal = _ag.default_non_eval(signal, la.asarray,
                                           la.linspace(1, -1, self.nplast))
    # -------------------------------------------------------------------------
    # Housekeeping
    # -------------------------------------------------------------------------

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = super().__repr__()
        insert = f"    weight = {self.weight!r},\n"
        insert += f"    signal = {self.signal!r},\n"
        rpr = (rpr[:-1] + insert + rpr[-1])
        return rpr

    def __str__(self) -> str:
        """Short representation of object"""
        rpr = super().__str__()
        rpr += f", w = {self.weight}"
        return rpr

    # -------------------------------------------------------------------------
    # Markov quantities
    # -------------------------------------------------------------------------

    def markov(self) -> la.lnarray:
        """Average transition rate matrix.

        .. math:: W^f = fp Wp + fm Wm
        """
        return (self.frac.s * self.plast).sum(-3)

    def enc(self) -> la.lnarray:
        """Average transition rate matrix, weighted by initial desired weight.

        .. math:: Q = fp Wp - fm Wm
        """
        return (self.signal.s * self.frac.s * self.plast).sum(-3)

    def zinv(self, rate: Optional[_sb.ArrayLike] = None,
             rowv: Optional[_sb.ArrayLike] = None) -> la.lnarray:
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
        onev = la.ones_like(self.weight)
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

    def cond(self, rate: Optional[_sb.ArrayLike] = None, *,
             rowv: Optional[_sb.ArrayLike] = None,
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

        where :math:`Z` is defined as

        .. math:: Z = (e \xi - W^f)^{-1},

        where :math:`e` is a vector of ones and :math:`\xi` is any row vector
        with :math:`\xi.e \neq 0`.
        When we include :math:`s`

        .. math:: Z(s) = (s I + e \xi - W^f)^{-1} \simeq \int e^{t(W^f-sI)} dt,

        i.e. :math:`Z^{-1}` with :math:`s` added to the diagonal.
        Effectively the matrix inverse of the genarator's Laplace transform.

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
        rowv = np.ones_like(self.weight)
        fundi = self.zinv(rowv=rowv)
        return rowv @ fundi.inv

    def mean_weight(self) -> float:
        """Mean synaptic weight in steady state distribution.

        Returns
        -------
        wbar : la.lnarray
            mean synaptic weight :math:`\\pi w`.
        """
        return self.peq() @ self.weight

    def deta(self, rate: Optional[_sb.ArrayLike] = None) -> la.lnarray:
        """Weighted Laplacian mean first passage time

        A measure of how far a state is from the potentiated states, related to
        Kemeney's constant:

        Parameters
        ----------
        rate : float, array, optional
            Parameter of Laplace transform, ``s``. Default: 0.

        Returns
        -------
        deta : la.lnarray
            Matrix of differences in ``eta(s)``

        Notes
        -----

        .. math:: \\eta_i = \\sum_j T_{ij}(s) \\pi_j w_j.
        .. math:: \\delta \\eta_{ij}(s) = \\eta_i(s) - \\eta_j(s).
        """
        eta = - self.zinv(rate).inv @ self.weight.c
        return eta - eta.t

    # -------------------------------------------------------------------------
    # Memory curves
    # -------------------------------------------------------------------------

    def spectrum(self) -> Tuple[la.lnarray, la.lnarray]:
        """Timescales and coefficients of eigenmodes.

        Returns
        -------
        taua : la.larray
            Decay timescales of eigenmodes.
        inita : la.lnarray
            Initial SNR from eigenmodes.
        """
        qas, evs = np.linalg.eig(-self.markov())
        mask = la.ones_like(qas, dtype=bool).put(qas.real.argmin(), False)
        taua = 1. / qas[mask]
        inita = self.peq() @ self.enc() @ evs[:, mask]
        inita *= (evs.inv @ self.weight)[mask]
        return taua, inita

    def snr(self, time: _sb.ArrayLike) -> la.lnarray:
        """Memory SNR as a function of recall time (memory curve).

        Parameters
        ----------
        time : float, array
            Recall time.
        """
        # convert to lnarray if needed, add singleton to broadcast with vector
        t_arr = la.array(time).c
        taua, inita = self.spectrum()
        return _sb.scalarise(np.sum(inita * np.exp(- t_arr / taua), axis=-1))

    def snr_laplace(self, rate: Optional[_sb.ArrayLike] = None) -> la.lnarray:
        """Laplace transform of SNR memory curve.

        Parameters
        ----------
        rate : float, array, optional
            Parameter of Laplace transform, ``s``. Default: 0.
        """
        a_drv = (self.enc() * self.deta(rate)).sum(-1)
        return _sb.scalarise((self.peq() * a_drv).sum(-1))

    def snr_exp_ave(self, tau: _sb.ArrayLike) -> la.lnarray:
        """Exponential running average of SNR memory curve.

        Mean SNR under exponential distribution of recall times.

        Parameters
        ----------
        tau : float, array
            Mean recall time.
        """
        return _sb.scalarise(self.snr_laplace(1. / tau) / tau)

    def snr_area(self) -> float:
        """Area under SNR memory curve.
        """
        return self.snr_laplace(None)

    def snr_init(self) -> float:
        """Initial value of SNR memory curve.
        """
        return _sb.scalarise(self.peq() @ self.enc() @ self.weight)

    @contextmanager
    def shifted(self, rate: _sb.ArrayLike, rowv: Optional[_sb.ArrayLike] = None
                ) -> None:
        """Move rate parametr of Laplace transform to plast
        """
        rowv = self.peq() if rowv is None else rowv
        # convert to lnarray, add singletons to broadcast with plast
        sss = _sb.insert_axes(la.asarray(rate), 3)
        # sss = la.asarray(rate).expand_dims((-3, -2, -1)) / (2 * self.frac.s)
        try:
            old_plast = self.plast.copy()
            self.plast = old_plast + sss * (rowv - la.eye(self.nstate))
            yield
        finally:
            self.plast = old_plast


# =============================================================================
# Helper functions
# =============================================================================


def sign_fix(model: SynapseMemoryModel):
    """Swap plasticity matrices if snr is negative"""
    if model.snr_init() < 0:
        model.signal *= -1


def normalise(model: SynapseMemoryModel):
    """Ensure that all attributes are valid.
    """
    _ma.stochastify_c(model.plast)
    scale = -np.diagonal(model.plast).min()
    if scale > 1:
        model /= scale
    _ma.stochastify_d(model.frac)


def valid_shapes(model: SynapseMemoryModel) -> bool:
    """Do attributes (plast, weight, frac) have correct shapes?"""
    vld = model.plast.shape[-2] == model.nstate
    vld &= model.frac.shape[-1] == model.nplast
    vld &= model.weight.shape[-1] == model.nstate
    vld &= model.signal.shape[-1] == model.nplast
    return vld


def valid_values(model: SynapseMemoryModel) -> bool:
    """Do attributes (plast, frac) have valid values?"""
    vld = _ma.isstochastic_c(model.plast, model.StochThresh)
    vld &= _ma.isstochastic_d(model.frac, model.StochThresh)
    return vld


def well_behaved(model: SynapseMemoryModel,
                 rate: Optional[float] = None, cond: bool = False) -> bool:
    """Do attributes plast have finite values, and is Zinv well conditioned?"""
    vld = np.isfinite(model.plast).all()
    if cond:
        vld &= model.cond(rate, rate_max=True) < model.CondThresh
    return vld
