# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:06:34 2017

@author: Subhy
"""

from numbers import Number
from typing import Tuple, Optional, Union
import numpy as np
from .builders import la
from .builders import stochastify_c as _stochastify_c
from .synapse_memory_model import SynapseMemoryModel as _SynapseMemoryModel
from .synapse_base import SynapseBase as _SynapseBase


class SynapseOptModel(_SynapseMemoryModel):
    """Class for complex synapses, suitable for optimisation

    Subclass of SynapseMemoryModel.
    Same constructor & attributes, only methods added.
    """

    def get_params(self) -> la.lnarray:
        """Independent parameters of transition matrices.

        Returns
        -------
        params : la.lnarray
            Vector of off-diagonal elements, potentiation before depression,
            each in order:
            Wp_01, Wp_02, ... Wp_0n-1, Wp_10, Wp_12, ... Wp_n-2,n-1,
            Wm_01, Wm_02, ... Wm_0n-1, Wm_10, Wm_12, ... Wm_n-2,n-1.

        See Also
        --------
        mat2params
        """
        return la.hstack((mat2params(self.plast[0]),
                          mat2params(self.plast[1])))

    def set_params(self, params: np.ndarray):
        """Transition matrices. from independent parameters

        Parameters
        ----------
        params : np.ndarray
            Vector of off-diagonal elements, potentiation before depression,
            each in order:
            Wp_01, Wp_02, ... Wp_0n-1, Wp_10, Wp_12, ... Wp_n-2n-1,
            Wm_01, Wm_02, ... Wm_0n-1, Wm_10, Wm_12, ... Wm_n-2n-1.

        See Also
        --------
        params2mat
        """
        num = params.shape[-1] // 2
        self.plast = la.stack((params2mat(params[:num]),
                               params2mat(params[num:])))

    def _derivs(self, rate: Optional[Number] = None,
                inv: bool = False) -> (Tuple[la.lnarray, la.lnarray],
                                       Tuple[la.lnarray, la.lnarray],
                                       Tuple[la.lnarray, ...]):
        """Gradients of Laplace transform of SNR memory curve.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, :math:`s`.
        inv : bool, default: False
            Should we compute matrix inverses?

        Returns
        -------
        rows : tuple(la.lnarray)
            (p,c),
        cols : tuple(la.lnarray)
            (eta,theta),
        mats : tuple(la.lnarray)
            if inv: (Z,Zs,ZQZs)
            if not inv: (Z^-1,Zs^-1)
        """
        fundi = self.zinv()
        fundis = self.zinv(rate)

        peq = self.peq()
        c_s = peq @ self.enc() @ fundis.inv

        eta = fundis.inv @ self.weight
        theta = fundi.inv @ self.enc() @ eta

        if inv:
            fund = fundi.inv()
            funds = fundis.inv()
            zqz = fundi.inv @ self.enc() @ fundis.inv
            return (peq, c_s), (eta, theta), (fund, funds, zqz)
        return (peq, c_s), (eta, theta), (fundi, fundis)

    def snr_grad(self, time: Number) -> (float, la.lnarray):
        """Gradient of SNR memory curve.

        Parameters
        ----------
        time : float
            Recall time.

        Returns
        -------
        func : float
            Value of ``SNR`` at ``time``.
        grad : la.lnarray (2n(n-1),)
            Gradient of ``SNR`` at ``time`` with respect to parameters.
        """
        peq = self.peq()

        eig = la.wrappers.wrap_several(np.linalg.eig)
        qas, evs = eig(self.markov)
        expqt = np.exp(qas * time)
        expwftw = (evs * expqt) @ (evs.inv @ self.weight)

        func = - peq @ self.enc() @ expwftw

        dsdq = _diagsub(peq.c * expwftw)
        dsdw = _diagsub(peq.c * (self.zinv().inv @ (self.enc() @ expwftw)))

        fab = qas.c - qas.r
        degenerate = np.fabs(fab) < self.DegThresh
        fab[degenerate] = 1.
        fab = (expqt.c - expqt.r) / fab
        fab[degenerate] = expqt[degenerate.nonzero()[0]] * time
        fab *= ((peq @ self.enc()) @ evs).c * (evs.inv @ self.weight)

        dsdw += _diagsub(evs.inv @ fab @ evs)

        grad = - np.hstack((mat2params(dsdq + self.frac[0] * dsdw),
                            mat2params(-dsdq + self.frac[1] * dsdw)))

        return func, grad

    def laplace_grad(self, rate: Optional[Number]) -> (float, la.lnarray):
        """Gradient of Laplace transform of SNR memory curve.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, ``s``.

        Returns
        -------
        func : float
            Value of ``snr_laplace`` at ``s``.
        grad : la.lnarray (2n(n-1),)
            Gradient of ``snr_laplace`` at ``s`` with respect to parameters.
        """
        # (p,c), (eta,theta)
        rows, cols = self._derivs(rate)[:2]
        func = - rows[1] @ self.weight

        dadq = _diagsub(rows[0].c * cols[0])
        dadw = _diagsub(rows[0].c * cols[1] + rows[1].c * cols[0])
        dadwp = dadq + self.frac[0] * dadw
        dadwm = -dadq + self.frac[1] * dadw
        grad = - np.hstack((mat2params(dadwp), mat2params(dadwm)))

        return func, grad

    def area_grad(self) -> (float, la.lnarray):
        """Gradient of Area under SNR memory curve.

        Returns
        -------
        func : float
            Value of ``snr_area``.
        grad : la.lnarray (2n(n-1),)
            Gradient of ``snr_area`` with respect to parameters.
        """
        return self.laplace_grad(None)

    def laplace_hess(self, rate: Number) -> la.lnarray:
        """Hessian of Laplace transform of SNR memory curve.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, ``s``.

        Returns
        -------
        hess : la.lnarray (2n(n-1),2n(n-1))
            Hessian of ``snr_laplace`` at ``s`` with respect to parameters.
        """
        # (p,c), (eta,theta), (Z,Zs,ZQZs)
        rows, cols, mats = self._derivs(rate, True)

        # (n,n,n,n)
        hessww = _dbl_diagsub(_outer3(rows[0], mats[0], cols[1])
                              + _outer3(rows[0], mats[2], cols[0])
                              + _outer3(rows[1], mats[1], cols[0])).sum(0)
        # (2,n,n,n,n)
        hesswq = _dbl_diagsub(_outer3(rows[0], mats[0], cols[0])
                              + _trnsp4(_outer3(rows[0], mats[1], cols[0])))

        # (n(n-1),n(n-1))
        hesspp = tens2mat(hessww + hesswq.sum(0)/self.frac[0])*self.frac[0]**2
        hesspm = tens2mat(hessww - hesswq[0]/self.frac[1]
                          + hesswq[1]/self.frac[0]) * self.frac[0]*self.frac[1]
        hessmm = tens2mat(hessww - hesswq.sum(0)/self.frac[1])*self.frac[1]**2
        # (2n(n-1),2n(n-1))
        return - np.block([[hesspp, hesspm],
                           [hesspm.T, hessmm]])

    def laplace_hessp(self, rate: Number,
                      other: _SynapseBase) -> la.lnarray:
        """Matrix product of snr_laplace's hessian and change in parameters.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, ``s``.
        other : synapse_base
            CHange in parameters stored as a synapse model.

        Returns
        -------
        hessp : la.lnarray (2n(n-1),)
            Change in parameters matrix-multiplied by hessian of
            ``snr_laplace`` at ``s`` with respect to parameters.
        """
        # (p,c), (eta,theta), (Zi,Zis)
        rows, cols, mats = self._derivs(rate)
        # ZQZs
        mats += (mats[0].inv @ self.enc() @ mats[1].inv,)

        def _hesspr(vecm, frc):
            """Hess vec products
            Parameters: vecm: lnarray (n,n)
            Returns: mat3: lnarray (n,n), mats4: lnarray (2,n,n)
            """
            # p Z theta V  -> (2,n,n) -> (n,n)
            # c Zs eta V   -> (2,n,n) -> (n,n)
            # p ZQZs eta V -> (n,n)
            h_ww = (_outerdiv3p(rows[0], mats[0], cols[1], vecm).sum(0)
                    + _outerdiv3p(rows[1], mats[1], cols[0], vecm).sum(0)
                    + _outer3ps(rows[0], mats[2], cols[0], vecm)) * frc
            # p Z eta V  -> (2,n,n)
            # p Zs eta V -> (2,n,n)
            h_wq = (_outerdiv3p(rows[0], mats[0], cols[0], vecm)
                    + np.flip(_outerdiv3p(rows[0], mats[1], cols[0], vecm), 0))
            # (n,n), (2,n,n)
            return h_ww, h_wq

        # (n,n), (2,n,n)
        hwwp, hwqp = _hesspr(other.plast[0], self.frac[0])
        hwwm, hwqm = _hesspr(other.plast[1], self.frac[1])

        frq = self.frac[0] / self.frac[1]
        # (n(n-1),)
        hessp = mat2params(hwwm - hwqm[0] + hwqm[1]/frq + hwwp + hwqp.sum(0))
        hessm = mat2params(hwwp + hwqp[0] - hwqp[1]*frq + hwwm - hwqm.sum(0))
        # (2n(n-1),)
        return - np.hstack((hessp * self.frac[0], hessm * self.frac[1]))

    @classmethod
    def from_params(cls, params: np.ndarray, *args, **kwdargs):
        """Buils SynapseOpt object from independent parameters

        Parameters
        ----------
        params : np.ndarray
            Vector of off-diagonal elements, potentiation before depression,
            each in order:
            mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2n-1.
        frac : float
            fraction of events that are potentiating, default=0.5.
        binary : bool
            is the weight vector binary? Otherwise it's linear. Default: False
        ...
            extra arguments passed to `cls.build` or `builders.build_empty`.

        Returns
        -------
        synobj
            SynapseOpt instance
        """
        nst = int(np.sqrt(2*params.size + 1) + 1) // 2
        self = cls.empty(nst, *args, npl=2, **kwdargs)
        self.set_params(params)
        return self


# =============================================================================
# %%* Independent parameter helper functions
# =============================================================================


def offdiag_inds(nst: int) -> la.lnarray:
    """Indices of independent parameters of transition matrix.

    Parameters
    ----------
    nst : int
        Number of states.

    Returns
    -------
    K : la.lnarray (n(n-1),)
        Vector of ravel indices of off-diagonal elements, in order:
        mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2,n-1.
    """
    # when put into groups of size nst+1:
    # M[0,0], M[0,1], M[0,2], ..., M[0,n-1], M[1,0],
    # M[1,1], M[1,2], ..., M[1,n-1], M[2,0], M[2,1],
    # ...
    # M[n-2,n-2], M[n-2,n-1], M[n-1,0], ..., M[n-1,n-2],
    # unwanted elements are 1st in each group
    k_1st = (nst+1) * np.arange(nst-1)  # ravel ind of 1st element in group
    k = np.arange(nst**2 - 1)  # exclude final unwanted element
    return np.delete(k, k_1st)


def mat2params(mat: la.lnarray) -> la.lnarray:
    """Independent parameters of transition matrix.

    Parameters
    ----------
    mat : la.lnarray (n,n)
        Continuous time stochastic matrix.

    Returns
    -------
    params : la.lnarray (n(n-1),)
        Vector of off-diagonal elements, in order:
        mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2,n-1.
    """
    nst = mat.shape[0]
    # m_00, m_01, m_02, ..., m_0n-1, m_10,
    # mat_12, ..., mat_n-2,n
    param = mat.flatten()
    return param[offdiag_inds(nst)]


@la.wrappers.wrap_one
def params2mat(params: la.lnarray) -> la.lnarray:
    """Transition matrix from independent parameters.

    Parameters
    ----------
    params : la.lnarray (n(n-1),)
        Vector of off-diagonal elements, in order:
        mat_01, mat_02, ..., mat_0n-1, mat10, mat_12, ..., mat_n-2,n-1.

    Returns
    -------
    mat : la.lnarray (n,n)
        Continuous time stochastic matrix.
    """
    nst = ((1. + np.sqrt(1. + 4. * params.size)) / 2.).astype(int)
    mat = np.empty(nst**2)
    mat[offdiag_inds(nst)] = params
    mat.resize((nst, nst))
    _stochastify_c(mat)
    return mat


def tens2mat(tens: la.lnarray) -> la.lnarray:
    """Convert 4th rank tensor to matrix for independent elements vector

    Parameters
    ----------
    tens : la.lnarray (n,n,n,n)
        4D tensor in state space.

    Returns
    -------
    mat : la.lnarray (n(n-1),n(n-1))
        matrix in independent parameter space (off-diagonal elements).
    """
    nst = tens.shape[0]
    mat = tens.reshape((nst**2, nst**2))
    k = offdiag_inds(nst)
    return mat[np.ix_(k, k)]


# =============================================================================
# %%* Private helper functions
# =============================================================================


def _diagsub(mat: la.lnarray) -> la.lnarray:
    """Subtract diagonal elements from each element of corresponding row
    Returns: mat2: lnarray (n,n)
    """
    return mat - mat.diagonal(0, -2, -1).c


def _outer3(vec1: la.lnarray, mat: la.lnarray, vec2: la.lnarray) -> la.lnarray:
    """Outer product of vector, matrix and vector.
    Returns: tens: lnarray (n,n,n,n)
    """
    return vec1[..., None, None, None] * mat[..., None] * vec2


def _outer3p(vec1: la.lnarray, mat1: Union[la.lnarray, la.invarray],
             vec2: la.lnarray, mat2: la.lnarray) -> la.lnarray:
    """Outer product of vector, matrix and vector, multiplying matrix
    Returns: mat1: lnarray (n,n), mat2: lnarray (n,n)
    """
    return (_diagsub(vec1.c * (mat1 @ (mat2 @ vec2))),
            _diagsub(((vec1 @ mat2) @ mat1).c * vec2))


def _outer3ps(vec1: la.lnarray, mat1: la.lnarray, vec2: la.lnarray,
              mat2: la.lnarray) -> la.lnarray:
    """Outer product of vector, matrix and vector, multiplying matrix
    Returns: mat3: lnarray (n,n)
    """
    mats = _outer3p(vec1, mat1, vec2, mat2)
    return mats[0] + mats[1]


def _outerdiv3p(vec1: la.lnarray, mat1: la.lnarray, vec2: la.lnarray,
                mat2: la.lnarray) -> la.lnarray:
    """Outer product of vector, inverse matrix and vector, multiplying matrix
    Returns: mats: lnarray (2,n,n)
    """
    return la.stack(_outer3p(vec1, mat1.inv, vec2, mat2))


def _trnsp4(tens: la.lnarray) -> la.lnarray:
    """Swap 1st two and last two indices of 4th rank tensor
    Returns: tens: lnarray (n,n,n,n)
    """
    return tens.transpose(2, 3, 0, 1)


def _dbl_diagsub(tens: la.lnarray) -> la.lnarray:
    """Subtract diagonal elements from each element of corresponding row
    for 1st two and last two indices.
    Returns: tens: lnarray (2,n,n,n,n)
    """
    tens_trn = _diagsub(_trnsp4(_diagsub(tens)))
    return la.stack(_trnsp4(tens_trn), tens_trn)
