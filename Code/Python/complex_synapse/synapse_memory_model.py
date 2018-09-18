# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 18:22:05 2017

@author: Subhy
"""


from typing import ClassVar, Optional
import numpy as np
from .synapse_base import SynapseBase, ArrayLike, la


class SynapseMemoryModel(SynapseBase):
    """Class for memory curves of complex synapse models.

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

    # synaptic weight of each state.
    weight: la.lnarray
    # desired signal contribution from each plasticity type
    signal: la.lnarray

    # Common constatnts / parameters

    # smallest reciprocal condition number for inverting zinv
    RCondThresh: ClassVar[float] = 1e-5
    # degeneracy threshold, for evals or eta^+
    DegThresh: ClassVar[float] = 1e-3
    # smallest singular value for Split models
    SingValThresh: ClassVar[float] = 1e-10
    # threshold for lumpability test
    LumpThresh: ClassVar[float] = 1e-3
    # threshold for orphaned states
    OrphanThresh: ClassVar[float] = 1e-3

    def __init__(self, plast: la.lnarray,
                 weight: la.lnarray,
                 frac: ArrayLike = 0.5,
                 signal: ArrayLike = np.nan,
                 **kwargs):
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
        # store inputs
        super().__init__(plast, frac)
        self.weight = la.asarray(weight)
        self.signal = la.asarray(signal).ravel()
        self.fix()

    def copy(self) -> 'SynapseBase':
        """Make copy of object
        """
        return type(self)(plast=self.plast.copy(),
                          weight=self.weight.copy(),
                          frac=self.frac.copy(),
                          signal=self.signal.copy())

    def fix(self):
        """Complete frac and signal vectors
        """
        super().fix()
        try:
            if np.isnan(self.signal).all():
                self.signal = np.linspace(1, -1, self.nplast)
        except AttributeError:
            pass  # probably called in super constructor

    def valid_shapes(self) -> bool:
        """Do attributes (Wp, Wm, w) have correct shapes?"""
        vld = super().valid_shapes()
        vld &= len(self.weight) == self.plast.shape[-1]
        vld &= len(self.frac) == len(self.signal)
        return vld

    def markov(self) -> la.lnarray:
        """Average transition rate matrix.

        .. math:: W^f = fp Wp + fm Wm
        """
        return (self.frac.s * self.plast).sum(0)

    def enc(self) -> la.lnarray:
        """Average transition rate matrix, weighted by initial desired weight.

        .. math:: Q = fp Wp - fm Wm
        """
        return (self.signal.s * self.frac.s * self.plast).sum(0)

    def zinv(self, rate: Optional[ArrayLike] = None,
             rowv: Optional[la.lnarray] = None) -> la.lnarray:
        r"""Inverse of generalised fundamental matrix.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, :math:`s`. Default: 0.
        rowv : la.lnarray, optional
            Arbitrary row vector, :math:`\xi`. Default: vector of ones.

        Returns
        -------
        Zi : la.lnarray
            inverse of generalised fundamental matrix, :math:`Z`.

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

        See Also
        --------
        zinv_s : specialised fundamental matrix
        """
        if rowv is None:
            rowv = la.ones_like(self.weight)
        fundi = la.ones((self.nstates, 1)) * rowv - self.markov()
        if rate is not None:
            # convert to ndarray, add singletons to broadcast with matrix
            s_arr = la.asarray(rate).s
            fundi += s_arr * np.eye(self.nstates)
        return fundi

    def peq(self) -> la.lnarray:
        r"""Steady state distribution.

        Returns
        -------
        peq : la.lnarray
            Steady state distribution, :math:`\pi`.
            Solution of: :math:`\pi W^f = 0`.
        """
        rowv = la.ones_like(self.weight)
        fundi = self.zinv(rowv=rowv)
        return rowv @ fundi.inv

    def zinv_s(self, rate: Optional[ArrayLike] = None) -> la.lnarray:
        r"""Inverse of special fundamental matrix.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, :math:`s`. Default: 0.

        Returns
        -------
        Zi : la.lnarray
            inverse of special fundamental matrix, :math:`Z`.

        Notes
        -----
        ``zinv_s`` is the inverse of :math:`Z`, defined as

        .. math:: Z = (e \pi - W^f)^{-1},

        i.e. ``zinv`` with special choice of :math:`\xi = \pi`,
        the steady state distribution.
        When we include :math:`s`

        .. math:: Z(s) = (s I + e \pi - W^f)^{-1} \simeq \int e^{t(W^f-sI)} dt,

        i.e. :math:`Z^{-1}` with :math:`s` added to the diagonal.
        Effectively the matrix inverse of the genarator's Laplace transform.

        See Also
        --------
        zinv : generalised fundamental matrix
        """
        return self.zinv(rate, self.peq())

    def mean_weight(self) -> float:
        """Mean synaptic weight in steady state distribution.

        Returns
        -------
        wbar : la.lnarray
            mean synaptic weight :math:`\\pi w`.
        """
        return self.peq() @ self.weight

    def deta(self, rate: Optional[ArrayLike] = None) -> la.lnarray:
        r"""Weighted Laplacian mean first passage time

        A measure of how far a state is from the potentiated states, related to
        Kemeney's constant:

        .. math:: \eta_i = \sum_j T_{ij} \pi_j w_j.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, :math:`s`. Default: 0.

        Returns
        -------
        deta : la.lnarray
            :math:`\delta \eta_{ij} = \eta_i - \eta_j`.
        """
        eta = - self.zinv(rate).inv @ self.weight
        return eta.c - eta

    def spectrum(self) -> (la.lnarray, la.lnarray):
        """Timescales and coefficients of eigenmodes.

        Returns
        -------
        taua : la.larray
            Decay timescales of eigenmodes.
        inita : la.lnarray
            Initial SNR from eigenmodes.
        """
        eig = la.wrappers.wrap_several(np.linalg.eig)
        qas, evs = eig(-self.markov())
        mask = np.ones_like(qas, dtype=bool).put(qas.real.argmin(), False)
        taua = 1. / qas[mask]
        inita = self.peq() @ self.enc() @ evs[:, mask]
        inita *= (evs.inv @ self.weight)[mask]
        return taua, inita

    def snr(self, time: ArrayLike) -> la.lnarray:
        """SNR memory curve

        Parameters
        ----------
        time : float, la.lnarray, optional
            Recall time, :math:`t`.
        """
        # convert to lnarray if needed, add singleton to broadcast with vector
        t_arr = la.array(time).c
        taua, inita = self.spectrum()
        return np.sum(inita * np.exp(- t_arr / taua), axis=-1)

    def snr_laplace(self, rate: Optional[ArrayLike] = None) -> la.lnarray:
        """Laplace transform of SNR memory curve.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, :math:`s`.
        """
        return np.sum(self.peq() @ (self.enc() * self.deta(rate)))

    def snr_exp_ave(self, tau: ArrayLike) -> la.lnarray:
        """Exponential running average of SNR memory curve.

        Mean SNR under exponential distribution of recall times.

        Parameters
        ----------
        tau : float, optional
            Mean recall time, :math:`\\tau`
        """
        return self.snr_laplace(1. / tau) / tau

    def snr_area(self) -> float:
        """Area under SNR memory curve.
        """
        return self.snr_laplace(None)

    def snr_init(self) -> float:
        """Initial value of SNR memory curve.
        """
        return self.peq() @ self.enc() @ self.weight
