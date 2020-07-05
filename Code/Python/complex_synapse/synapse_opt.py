# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:06:34 2017

@author: Subhy
"""
from __future__ import annotations

from numbers import Number
from typing import ClassVar, Dict, Optional, Tuple, Union

import numpy as np

import numpy_linalg as la
from sl_py_tools.numpy_tricks import markov_param as mp

from . import builders as bld
from .synapse_base import SynapseBase as _SynapseBase
from .synapse_memory_model import SynapseMemoryModel as _SynapseMemoryModel
from .synapse_memory_model import well_behaved

wrap = la.wrappers.Wrappers(la.lnarray)
eig = wrap.several(np.linalg.eig)


class SynapseOptModel(_SynapseMemoryModel):
    """Class for complex synapses, suitable for optimisation

    Subclass of SynapseMemoryModel.
    Same constructor & attributes, only methods added.
    """
    _saved: Tuple[Tuple[la.lnarray, ...], ...]
    _unchanged: bool
    type: ClassVar[Dict[str, bool]] = {'serial': False,
                                       'ring': False,
                                       'uniform': False}
    directions: ClassVar[Tuple[int, ...]] = (0, 0)

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self._saved = ((), (), (), ())
        self._unchanged = False

    @property
    def nparam(self) -> int:
        """number of plasticity types."""
        npr = mp.num_param(self.nstates, drn=self.directions[0], **self.type)
        return self.nplast * npr

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
        markov.mat_to_params
        """
        params = [mp.mat_to_params(self.plast[i], drn=d, **self.type)
                  for i, d in enumerate(self.directions)]
        return la.hstack(params)

    def set_params(self, params: np.ndarray, *args, **kwds):
        """Transition matrices. from independent parameters

        Does not update if parameters are unchanged Optional arguments are
        passed on to `numpy.allclose`.

        Parameters
        ----------
        params : np.ndarray
            Vector of off-diagonal elements, potentiation before depression,
            each in order:
            Wp_01, Wp_02, ... Wp_0n-1, Wp_10, Wp_12, ... Wp_n-2n-1,
            Wm_01, Wm_02, ... Wm_0n-1, Wm_10, Wm_12, ... Wm_n-2n-1.

        See Also
        --------
        markov.params_to_mat
        """
        if np.allclose(params, self.get_params(), *args, **kwds):
            return
        plast = np.split(params, len(self.directions))
        for i, drn in enumerate(self.directions):
            mp.mat_update_params(self.plast[i], plast[i], drn=drn, **self.type)
        self._unchanged = False

    def _derivs(self, rate: Optional[Number] = None,
                inv: bool = False) -> Tuple[Tuple[la.lnarray, ...], ...]:
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
            if not inv: (Z^{-1},Zs^{-1})
        """
        if self._unchanged and (rate, inv) == self._saved[0]:
            return self._saved[1:]
        self._unchanged = True

        fundi = self.zinv()
        peq = self.peq()

        if rate is None:
            c_s = peq @ self.enc() @ fundi.inv
            eta = fundi.inv @ self.weight
            theta = fundi.inv @ self.enc() @ eta
            if inv:
                mats = (fundi.inv(),)
            else:
                mats = (fundi,)
            self._saved = (rate, inv), (peq, c_s), (eta, theta), mats
            return self._saved[1:]

        fundis = self.zinv(rate)
        c_s = peq @ self.enc() @ fundis.inv

        eta = fundis.inv @ self.weight
        theta = fundi.inv @ self.enc() @ eta

        if inv:
            zqz = fundi.inv @ self.enc() @ fundis.inv
            mats = (fundi.inv(), fundis.inv(), zqz)
        else:
            mats = (fundi, fundis)

        self._saved = (rate, inv), (peq, c_s), (eta, theta), mats
        return self._saved[1:]

    def _calc_fab(self, time, qas, evs, peq_enc):
        """Calculate eigenvalue differences
        """
        expqt = np.exp(qas * time)
        fab = qas.c - qas.r
        degenerate = np.fabs(fab) < self.DegThresh
        fab[degenerate] = 1.
        fab = (expqt.c - expqt.r) / fab
        fab[degenerate] = expqt[degenerate.nonzero()[0]] * time
        fab *= evs.inv @ self.weight
        fab *= (peq_enc @ evs).c
        expwftw = (evs * expqt) @ (evs.inv @ self.weight)
        return (evs.inv @ fab @ evs), expwftw

    def snr_grad(self, time: Number) -> (float, la.lnarray):
        """Gradient of instantaneous SNR memory curve.

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
        peq_enc = peq @ self.enc()

        evd = eig(self.markov)
        fab, expqtw = self._calc_fab(time, *evd, peq_enc)

        func = - peq_enc @ expqtw

        dsdq = - _diagsub(peq.c * expqtw)
        dsdw = - _diagsub(peq.c * (self.zinv().inv @ (self.enc() @ expqtw)))
        dsdw -= _diagsub(fab)

        dsdwp = self.frac[0] * (dsdw + dsdq)
        dsdwm = self.frac[1] * (dsdw - dsdq)

        grad = [mp.mat_to_params(m, drn=d, **self.type)
                for m, d in zip([dsdwp, dsdwm], self.directions)]
        grad = np.hstack(grad)

        return bld.scalarise(func), grad

    def laplace_fun(self, rate: Number, inv: bool = False) -> float:
        """Gradient of Laplace transform of SNR memory curve.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, ``s``.
        inv : bool, default: False
            Should we compute matrix inverses?

        Returns
        -------
        func : float
            Value of ``snr_laplace`` at ``s``.
        """
        if not well_behaved(self, rate):
            return np.array(1.e10), bld.RNG.random(self.nparam)
        # (p,c), (eta,theta)
        rows, cols, mats = self._derivs(rate, inv)

        func = - rows[1] @ self.weight
        # afunc = -rows[0] @ self.enc() @ cols[0]

        return bld.scalarise(func)

    def laplace_grad(self, rate: Number,
                     inv: bool = False) -> (float, la.lnarray):
        """Gradient of Laplace transform of SNR memory curve.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, ``s``.
        inv : bool, default: False
            Should we compute matrix inverses?

        Returns
        -------
        func : float
            Value of ``snr_laplace`` at ``s``.
        grad : la.lnarray (2n(n-1),)
            Gradient of ``snr_laplace`` at ``s`` with respect to parameters.
        """
        if not well_behaved(self, rate):
            return np.array(1.e10), bld.RNG.random(self.nparam)
        # (p,c), (eta,theta)
        rows, cols, mats = self._derivs(rate, inv)

        func = - rows[1] @ self.weight
        # afunc = -rows[0] @ self.enc() @ cols[0]

        dadq = - _diagsub(rows[0].c * cols[0])
        dadw = - _diagsub(rows[0].c * cols[1] + rows[1].c * cols[0])
        dadwp = self.frac[0] * (dadw + dadq)
        dadwm = self.frac[1] * (dadw - dadq)
        grad = [mp.mat_to_params(m, drn=d, **self.type)
                for m, d in zip([dadwp, dadwm], self.directions)]
        grad = np.hstack(grad)

        return bld.scalarise(func), grad

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
        hesspp = mp.tens_to_mat(hessww + hesswq.sum(0),
                                drn=(self.directions[0],) * 2,
                                **self.type) * self.frac[0]**2
        hesspm = mp.tens_to_mat(hessww - hesswq[0] + hesswq[1],
                                drn=self.directions,
                                **self.type) * self.frac[0]*self.frac[1]
        hessmm = mp.tens_to_mat(hessww - hesswq.sum(0),
                                drn=(self.directions[1],) * 2,
                                **self.type) * self.frac[1]**2
        # (2n(n-1),2n(n-1))
        return - np.block([[hesspp, hesspm], [hesspm.T, hessmm]])

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
        rows, cols, mats = self._derivs(rate, False)
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
                    + _outer3ps(rows[0], mats[2], cols[0], vecm))
            # p Z eta V  -> (2,n,n)
            # p Zs eta V -> (2,n,n)
            h_wq = (_outerdiv3p(rows[0], mats[0], cols[0], vecm)
                    + np.flip(_outerdiv3p(rows[0], mats[1], cols[0], vecm), 0))
            # (n,n), (2,n,n)
            return -frc * h_ww, -frc * h_wq

        # (n,n), (2,n,n)
        hwwp, hwqp = _hesspr(other.plast[0], self.frac[0])
        hwwm, hwqm = _hesspr(other.plast[1], self.frac[1])

        # (n,n)
        hessp = hwwm - hwqm[0] + hwqm[1] + hwwp + hwqp.sum(0)
        hessm = hwwp + hwqp[0] - hwqp[1] + hwwm - hwqm.sum(0)
        # (2n(n-1),)
        hsv = [mp.mat_to_params(f * m, drn=d, **self.type)
               for f, m, d in zip(self.frac, [hessp, hessm], self.directions)]
        return np.hstack(hsv)

    def peq_min_fun(self, rate: Number) -> la.lnarray:
        """Gradient of lower bound of shifted Laplace problem.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, ``s``.
        inv : bool, default: False
            Should we compute matrix inverses?

        Returns
        -------
        func : la.lnarray (2n**2,)
            Value of ``W - s * e @ peq``.
        """
        if not np.isfinite(self.plast).all():
            return -np.ones(self.nparam)
        rate = rate / (2 * self.frac.s)
        # (p,c), (eta,theta), (Z,Zs,ZQZs)
        rows, cols, mats = self._derivs(None, True)
        func = self.plast - rate * rows[0] + (1 + rate) * np.eye(self.nstates)
        return func.flattish(-3)

    def peq_min_grad(self, rate: Number) -> la.lnarray:
        """Gradient of lower bound of shifted Laplace problem.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, ``s``.
        inv : bool, default: False
            Should we compute matrix inverses?

        Returns
        -------
        grad : la.lnarray (2n**2,2n(n-1),)
            Gradient of ``func`` with respect to parameters.
        """
        if not np.isfinite(self.plast).all():
            return -bld.RNG.random((self.nplast*self.nstates**2, self.nparam,))
        # (p,c), (eta,theta), (Z,Zs,ZQZs)
        rows, cols, mats = self._derivs(None, True)
        peq, fund = rows[0], mats[0]
        # (n,1,n,n)
        grad = np.moveaxis(_diagsub(- rate/2 * peq[None, :].s * fund), -1, -4)
        # (2,1,n,2,n,n)
        grad = grad * (self.frac / self.frac.c).expand_dims((1, 2, 4, 5))
        # (2,n,n,2,n,n)
        shape = ((self.nplast,) + (self.nstates,) * 2) * 2
        gradd = _diagsub(la.eye(self.nplast * self.nstates**2).reshape(shape))
        # (2n**2,2,n,n)
        grad = (gradd + grad).flattish(0, -3)
        # (2n**2,n(n-1)), (2n**2,n(n-1))
        opts = {'grad': True, 'axes': (-2, -1)}
        opts.update(self.type)
        gradp = mp.mat_to_params(grad[:, 0], drn=self.directions[0], **opts)
        gradm = mp.mat_to_params(grad[:, 1], drn=self.directions[1], **opts)
        # (2n**2,2n(n-1))
        return np.c_[gradp, gradm]

    def peq_min_hess(self, rate: Number, lag: np.ndarray) -> la.lnarray:
        """Hessian of lower bound of shifted Laplace problem.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, ``s``.
        lag : ndarray
            Lagrange multipliers

        Returns
        -------
        hess : la.lnarray (2n(n-1),2n(n-1))
            Hessian of ``func @ lag`` with respect to parameters.
        """
        if np.isnan(self.plast).any():
            return -bld.RNG.random((self.nparam,) * 2)
        # (p,c), (eta,theta), (Z,Zs,ZQZs)
        rows, cols, mats = self._derivs(None, True)
        peq, fund = rows[0], mats[0]
        # (n,n,n,n,1,1,n), constraint axes last three
        hess = (((-rate * peq).s * fund).s * fund).expand_dims((-3, -2))
        # (n,n,n,n,2,1,n)
        hess = hess / (2 * self.frac.s)
        # (n,n,n,n,2,n,n)
        hess = np.broadcast_to(hess, hess.shape[:-2] + (self.nstates,) * 2)
        # (n,n,n,n,2n**2) @ (2n**2,) -> (n,n,n,n)
        hess = hess.flattish(-3) @ lag
        hess = _dbl_diagsub(hess).sum(0)
        # (n(n-1),n(n-1))
        hesspp = mp.tens_to_mat(hess, drn=(self.directions[0],) * 2,
                                **self.type) * self.frac[0]**2
        hesspm = mp.tens_to_mat(hess, drn=self.directions,
                                **self.type) * self.frac[0]*self.frac[1]
        hessmm = mp.tens_to_mat(hess, drn=(self.directions[1],) * 2,
                                **self.type) * self.frac[1]**2
        # (2n(n-1),2n(n-1))
        return np.block([[hesspp, hesspm], [hesspm.T, hessmm]])

    @classmethod
    def from_params(cls, params: np.ndarray, *args, **kwds) -> SynapseOptModel:
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
            extra arguments passed to `cls.zero`.

        Returns
        -------
        synobj
            SynapseOpt instance
        """
        npl = kwds.pop('npl', len(cls.directions))
        nst = kwds.pop('nst', None)
        if nst is None:
            mat_type = cls.type.copy()
            del mat_type['uniform']
            mat_type['drn'] = cls.directions[0]
            nst = mp.num_state(params.size // npl, **mat_type)
        self = cls.zero(nst, *args, npl=npl, **kwds)
        self.set_params(params)
        return self

    @classmethod
    def rand(cls, nst: Optional[int], *args, **kwds) -> SynapseOptModel:
        """Random model

        Synapse model with random transition matrices

        Parameters
        ----------
            All passed to `cls.build`, `builders.build_rand`
            or `builders.rand_trans`:
        nst : int or None
            total number of states. Use `npar` if `None`.
        npar : int or None
            total number of parameters. Use `nst` if `None` (default).
        npl : int
            total number of plasticity types
        frac : float
            fraction of events that are potentiating, default=0.5.
        binary : bool
            is the weight vector binary? Otherwise it's linear. Default: False
        sp : float
            sparsity, default: 1.
        ...
            extra arguments passed to `cls.from_params`.

        One of 'nst/npar' must be provided

        Returns
        -------
        synobj
            SynapseOpt instance
        """
        npar = kwds.pop('npar', None)
        if npar is None:
            npl = kwds.get('npl', len(cls.directions))
            npar = npl * mp.num_param(nst, drn=cls.directions[0], **cls.type)
        return cls.from_params(bld.RNG.random(npar), *args, nst=nst, **kwds)


# =============================================================================
# Private helper functions
# =============================================================================


def _diagsub(mat: la.lnarray) -> la.lnarray:
    """Subtract diagonal elements from each element of corresponding row
    Returns: mat2: lnarray (...,n,n)
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
    return la.stack((_trnsp4(tens_trn), tens_trn))


if __name__ == "__main__":
    pass
