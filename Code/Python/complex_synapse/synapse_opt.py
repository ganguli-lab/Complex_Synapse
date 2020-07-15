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
from .builders import insert_axes
from .synapse_base import SynapseBase as _SynapseBase
from .synapse_memory_model import SynapseMemoryModel as _SynapseMemoryModel
from .synapse_memory_model import well_behaved
# =============================================================================


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
        return self.nplast * mp.num_param(self.nstates, **self.directed(0))

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
        params = mp.mat_to_params(self.plast, **self.directed(daxis=-3))
        return params.flattish(-2)

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
        params = la.asarray(params).foldaxis(-1, (self.nplast, -1))
        mp.mat_update_params(self.plast, params, **self.directed(pdaxis=-2))
        self._unchanged = False

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

        evd = np.linalg.eig(self.markov())
        fab, expqtw = self._calc_fab(time, *evd, peq_enc)

        func = - peq_enc @ expqtw

        # (n,n)
        dsdq = - _diagsub(peq.c * expqtw)
        dsdw = - _diagsub(peq.c * (self.zinv().inv @ (self.enc() @ expqtw)))
        dsdw -= _diagsub(fab)

        # (2,n,n)
        grad = self.frac.s * (dsdw + self.signal.s * dsdq)
        # (2,n(n-1))
        grad = mp.mat_to_params(grad, **self.directed(daxis=-3))
        # (2n(n-1),)
        grad = grad.flattish(-2)

        return bld.scalarise(func), grad

    def _derivs(self, rate: Optional[Number] = None, inv: bool = False,
                svd: bool = False) -> Tuple[Tuple[la.lnarray, ...], ...]:
        """Gradients of Laplace transform of SNR memory curve.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, :math:`s`.
        inv : bool, default: False
            Should we compute matrix inverses?
        svd : bool, default: False
            Should we compute SVD (1st & last only)

        Returns
        -------
        rows : tuple(la.lnarray)
            (p,c) + (U.h if svd)
        cols : tuple(la.lnarray)
            (eta,theta) + (V if svd)
        mats : tuple(la.lnarray)
            if not inv: (Z^{-1},Zs^{-1}) + (S if svd)
            if inv: (Z,Zs,ZQZs) + (if inv)
        """
        if self._unchanged and (rate, inv) == self._saved[0]:
            return self._saved[1:]
        self._unchanged = True

        args = (rate, inv, svd)
        fundi = self.zinv()
        peq = self.peq()

        if svd:
            extract = np.s_[::self.nstates-1]
            swap = la.array([[0, 1], [-1, 0]])
            fund_u, fund_s, fund_v = np.linalg.svd(fundi)
            # svs = diag(smin, -smax) / smin**2
            svs = (np.diagflat(swap @ fund_s[extract] / fund_s[-1]**2),)
            srow = (fund_u.h[extract],)
            scol = (fund_v.h[:, extract],)
        else:
            svs, scol, srow = (), (), ()

        if rate is None:
            c_s = peq @ self.enc() @ fundi.inv
            eta = fundi.inv @ self.weight
            theta = fundi.inv @ self.enc() @ eta
            mats = (fundi, fundi) + svs
            if inv:
                zqz = fundi.inv @ self.enc() @ fundi.inv
                mats = (fundi.inv(),) * 2 + (zqz,) + mats
            self._saved = args, (peq, c_s) + srow, (eta, theta) + scol, mats
            return self._saved[1:]

        fundis = self.zinv(rate)
        c_s = peq @ self.enc() @ fundis.inv

        if svd:
            fund_u, fund_s, fund_v = np.linalg.svd(fundis)
            # svs = diag(smin, -smax) / smin**2
            svs += (np.diagflat(swap @ fund_s[extract] / fund_s[-1]**2),)
            srow += (fund_u.h[extract],)
            scol += (fund_v.h[:, extract],)
            svs, srow, scol = [(np.stack(arr),) for arr in (svs, srow, scol)]

        eta = fundis.inv @ self.weight
        theta = fundi.inv @ self.enc() @ eta

        mats = (fundi, fundis) + svs
        if inv:
            zqz = fundi.inv @ self.enc() @ fundis.inv
            mats = (fundi.inv(), fundis.inv(), zqz) + mats

        self._saved = args, (peq, c_s) + srow, (eta, theta) + scol, mats
        return self._saved[1:]

    def laplace_fun(self, rate: Number, **kwds) -> float:
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
            return np.array(1.e10)
        # (p,c), (eta,theta)
        rows = self._derivs(rate, **kwds)[0]

        func = - rows[1] @ self.weight
        # afunc = -rows[0] @ self.enc() @ cols[0]

        return bld.scalarise(func)

    def laplace_grad(self, rate: Number, **kwds) -> la.lnarray:
        """Gradient of Laplace transform of SNR memory curve.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, ``s``.
        inv : bool, default: False
            Should we compute matrix inverses?

        Returns
        -------
        grad : la.lnarray (2n(n-1),)
            Gradient of ``snr_laplace`` at ``s`` with respect to parameters.
        """
        if not well_behaved(self, rate):
            return bld.RNG.random(self.nparam)
        # (p,c), (eta,theta)
        rows, cols, _ = self._derivs(rate, **kwds)

        # (n,n)
        dadq = - _diagsub(rows[0].c * cols[0])
        dadw = - _diagsub(rows[0].c * cols[1] + rows[1].c * cols[0])
        # (2,n,n)
        grad = self.frac.s * (dadw + self.signal.s * dadq)
        # (2,n(n-1))
        grad = mp.mat_to_params(grad, **self.directed(daxis=-3))
        # (2n(n-1),)
        return grad.flattish(-2)

    def laplace_hess(self, rate: Number, **kwds) -> la.lnarray:
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
        kwds['inv'] = True
        rows, cols, mats = self._derivs(rate, **kwds)

        # (...,n,n,n,n)
        hessww = (_outer3(rows[0], mats[0], cols[1])
                  + _outer3(rows[0], mats[2], cols[0])
                  + _outer3(rows[1], mats[1], cols[0])).sum(0)
        # (2,...,n,n,n,n) - axis 0: h_wq, h_qw
        hesswq = (_outer3(rows[0], mats[0], cols[0])
                  + _trnsp4(_outer3(rows[0], mats[1], cols[0])))
        # (...,1,n,n,1,n,n), (2,...,1,n,n,1,n,n)
        nax = (-6, -3)
        hessww, hesswq = hessww.expand_dims(nax), hesswq.expand_dims(nax)
        # (2,1,1,2,1,1)
        fmfn = insert_axes(self.frac, 5) * self.frac.s
        snu = self.signal.s
        smu = insert_axes(snu, 3)
        # (...,2,n,n,2,n,n)
        hess = fmfn * (hessww + snu * hesswq[0] + smu * hesswq[1])
        # (...,2,n(n-1),2,n(n-1))
        hess = mp.mat_to_params(hess, **self.directed(), daxis=(-6, -3),
                                axes=((-5, -4), (-2, -1)))
        return -hess.flattish(-4, -2).flattish(-2)

    def laplace_hessp(self, rate: Number, other: _SynapseBase,
                      **kwds) -> la.lnarray:
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
        kwds.setdefault('inv', False)
        rows, cols, mats = self._derivs(rate, **kwds)
        if kwds['inv']:
            mats = mats[3:5] + mats[2:3]
        else:
            # ZQZs
            mats = mats[:2] + (mats[0].inv @ self.enc() @ mats[1].inv,)

        # (2,n,n)
        vecm = self.frac.s * other.plast
        # p Z theta V  -> (2,2,n,n) -> (n,n)
        # c Zs eta V   -> (2,2,n,n) -> (n,n)
        # p ZQZs eta V -> (2,2,n,n) -> (n,n)
        h_ww = (_outer3p(rows[0], mats[0].inv, cols[1], vecm)
                + _outer3p(rows[1], mats[1].inv, cols[0], vecm)
                + _outer3p(rows[0], mats[2], cols[0], vecm)).sum((0, -3))
        # p Z eta V  -> (2,2,n,n) - h_wq, h_qw
        # p Zs eta V -> (2,2,n,n) - h_wq, h_qw
        h_wqqw = (_outer3p(rows[0], mats[0].inv, cols[0], vecm)
                  + _outer3p(rows[0], mats[1].inv, cols[0], vecm)[::-1])
        # ^axis -3: hess @ other.wp, hess @ other.wm
        # ^axis 0: hess_wq @ other.w_, hess_qw @ other.w_
        # (n,n)
        h_wq = (self.signal.s * h_wqqw[0]).sum(-3)
        h_qw = h_wqqw[1].sum(-3)

        # (2,n,n)
        hsv = self.frac.s * (h_ww + h_wq + self.signal.s * h_qw)
        # (2,n(n-1))
        hsv = mp.mat_to_params(hsv, **self.directed(daxis=-3))
        # (2n(n-1),)
        return hsv.flattish(-2)

    def area_fun(self, **kwds) -> float:
        """Area under SNR memory curve.

        Returns
        -------
        func : float
            Value of ``snr_area``.
        """
        return self.laplace_fun(None, **kwds)

    def area_grad(self, **kwds) -> la.lnarray:
        """Gradient of Area under SNR memory curve.

        Returns
        -------
        grad : la.lnarray (2n(n-1),)
            Gradient of ``snr_area`` with respect to parameters.
        """
        return self.laplace_grad(None, **kwds)

    def area_hess(self, **kwds) -> la.lnarray:
        """Hessian of Area under SNR memory curve.

        Returns
        -------
        hess : la.lnarray (2n(n-1),2n(n-1))
            Hessian of ``snr_area`` with respect to parameters.
        """
        return self.laplace_hess(None, **kwds)

    def peq_min_fun(self, rate: Number, **kwds) -> la.lnarray:
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
        if not well_behaved(self, rate):
            return -np.ones(self.nplast*self.nstates**2)
        rate = rate / (2 * self.frac.s)
        # (p,c), (eta,theta), (Z,Zs,ZQZs)
        kwds['inv'] = True
        rows, cols, mats = self._derivs(None, **kwds)
        func = self.plast - rate * rows[0] + (1 + rate) * np.eye(self.nstates)
        return func.flattish(-3)

    def peq_min_grad(self, rate: Number, **kwds) -> la.lnarray:
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
        if not well_behaved(self, rate):
            return -bld.RNG.random((self.nplast*self.nstates**2, self.nparam,))
        # (p,c), (eta,theta), (Z,Zs,ZQZs)
        kwds['inv'] = True
        rows, cols, mats = self._derivs(None, **kwds)
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
        # (2n**2,2,n(n-1))
        grad = mp.mat_to_params(grad, **self.directed(daxis=-3))
        # (2n**2,2n(n-1))
        return grad.flattish(-2)
        # # (2n**2,n(n-1)), (2n**2,n(n-1))

    def peq_min_hess(self, rate: Number, lag: np.ndarray,
                     **kwds) -> la.lnarray:
        """Hessian of lower bound of shifted Laplace problem.

        Parameters
        ----------
        rate : float, optional
            Parameter of Laplace transform, ``s``.
        lag : ndarray (2n**2,)
            Lagrange multipliers

        Returns
        -------
        hess : la.lnarray (2n(n-1),2n(n-1))
            Hessian of ``peq_min_fun(rate) @ lag`` with respect to parameters.
        """
        if not well_behaved(self, rate):
            return bld.RNG.random((self.nparam,) * 2)
        # (p,c), (eta,theta), (Z,Zs,ZQZs)
        kwds['inv'] = True
        rows, cols, mats = self._derivs(None, **kwds)
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
        # (2,n,n,2,n,n)
        hess = hess.expand_dims(-3) * insert_axes(self.frac, 5) * self.frac.s
        # (2,n(n-1),2,n(n-1))
        hess = mp.mat_to_params(hess, **self.directed(
            axes=((-5, -4), (-2, -1)), daxis=(-6, -3)))
        # (2n(n-1),2n(n-1))
        return hess.flattish(-4, -2).flattish(-2)

    def cond_fun(self, rate: Number, **kwds) -> float:
        """Gap between condition number and threshold"""
        if not well_behaved(self, rate):
            nout = 1 if rate is None else 2
            return -np.ones(nout)
        if kwds.pop('svd', False):
            # (p,c,U.h), (eta,theta,V), (Z,Zs,ZQZs,S)
            svs = self._derivs(rate, svd=True, **kwds)[-1][-1]
            # svs = diag(smin, -smax) / smin**2
            return 1 + svs[..., 1, 1] / svs[..., 0, 0] / self.CondThresh
        if rate is None:
            return 1 - self.cond(None) / self.CondThresh
        conds = la.array([self.cond(val) for val in (None, rate)])
        return 1 - conds / self.CondThresh

    def cond_grad(self, rate: Number, **kwds) -> float:
        """Gap between condition number and threshold"""
        if not well_behaved(self, rate):
            nout = 1 if rate is None else 2
            return bld.RNG.random((nout, self.nparam))
        kwds['svd'] = True
        # (p,c,U.h), (eta,theta,V), (Z,Zs,ZQZs,S)
        row, col, svs = [drv[-1] for drv in self._derivs(rate, **kwds)]
        # svs = diag(smin, -smax) / smin**2
        grad = row.t @ svs @ col.t
        # (...,2,n,n)
        grad = _diagsub(grad).expand_dims(-3) * self.frac.s
        # (2,n(n-1))
        grad = mp.mat_to_params(grad, **self.directed(daxis=-3))
        # (2n(n-1),)
        return grad.flattish(-2) / self.CondThresh
        # grad = [mp.mat_to_params(grad[..., i, :, :], drn=d, **self.type)
        #         for i, d in enumerate(self.directions)]
        # return np.hstack(grad) / self.CondThresh

    @classmethod
    def directed(cls, which: Union[int, slice] = np.s_[:], **kwds) -> dict:
        """Markov parameter options

        Parameters
        ----------
        which : int, slice, optional
            Which element of `cls.directions` to use, by default
        extra arguments are default values or unknown keys in `opts`

        Returns
        -------
        opts
            Dictionary with `cls.type`, `drn = cls.directions[which]`.
        """
        kwds['drn'] = cls.directions[which]
        kwds.update(cls.type)
        return kwds

    @classmethod
    def from_params(cls, params: np.ndarray, *args, **kwds) -> SynapseOptModel:
        """Builds SynapseOpt object from independent parameters

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
            nst = mp.num_state(params.size // npl, **cls.directed(0))
        obj = cls.zero(nst, *args, npl=npl, **kwds)
        obj.set_params(params)
        return obj

    @classmethod
    def rand(cls, nst: Optional[int], *args, **kwds) -> SynapseOptModel:
        """Random model

        Synapse model with random transition matrices

        Parameters
        ----------
            All passed to `cls.build`, `builders.build_rand`
            or `sl_py_tools.numpy_tricks.markov_param`:
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
            npar = npl * mp.num_param(nst, **cls.directed(0))
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
    Returns: tens: lnarray (2,...,n,n,n,n) - 2nd array is transpose of 1st
    """
    return _dbl_diagsub(insert_axes(vec1, 3) * mat.expand_dims((-4, -1))
                        * insert_axes(vec2, 3, -2))


def _outer3p(vec1: la.lnarray, mat1: Union[la.lnarray, la.invarray],
             vec2: la.lnarray, mat2: la.lnarray) -> la.lnarray:
    """Outer product of vector, matrix and vector, multiplying matrix
    Returns: mats: lnarray (2,...,n,n) - tensor @ mat2, tensor.t* @ mat2
    where tensor = _outer3(vec1,mat1,vec2) and tensor.t* = _trnsp4(tensor)
    """
    return _diagsub(np.stack((vec1.c * (mat1 @ (mat2 @ vec2.c)).uc.r,
                              ((vec1.r @ mat2) @ mat1).ur.c * vec2.r)))


def _trnsp4(tens: la.lnarray) -> la.lnarray:
    """Swap 1st two and last two indices of 4th rank tensor
    Returns: tens: lnarray (...,n,n,n,n)
    """
    return tens.moveaxis((-4, -3), (-2, -1))


def _dbl_diagsub(tens: la.lnarray) -> la.lnarray:
    """Subtract diagonal elements from each element of corresponding row
    for 1st two and last two indices.
    Returns: tens: lnarray (2,...,n,n,n,n) - 2nd array is transpose of 1st
    """
    tens_trn = _diagsub(_trnsp4(_diagsub(tens)))
    return np.stack((_trnsp4(tens_trn), tens_trn))


if __name__ == "__main__":
    pass
