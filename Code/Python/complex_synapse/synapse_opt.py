# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 16:06:34 2017

@author: Subhy
"""
from __future__ import annotations
from numbers import Number
from typing import Tuple, Optional, Union, Dict, ClassVar
import numpy as np
import numpy_linalg as la
from sl_py_tools.numpy_tricks import markov_param as mp
from .builders import scalarise
from .synapse_memory_model import SynapseMemoryModel as _SynapseMemoryModel
from .synapse_base import SynapseBase as _SynapseBase

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
        self._saved = (0, (), (), ())
        self._unchanged = False

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
        # mats = [mp.params_to_mat(plast[i], drn=d, **self.type)
        #           for i, d in enumerate(self.directions)]
        # self.plast = la.stack(mats)
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
        if all((self._unchanged,
                len(self._saved[3]) == 2 + inv,
                np.isclose(rate, self._saved[0]))):
            return self._saved[1:]
        self._unchanged = True

        fundi = self.zinv()
        fundis = self.zinv(rate)

        peq = self.peq()
        c_s = peq @ self.enc() @ fundis.inv

        eta = fundis.inv @ self.weight
        theta = fundi.inv @ self.enc() @ eta

        if inv:
            zqz = fundi.inv @ self.enc() @ fundis.inv
            mats = (fundi.inv(), fundis.inv(), zqz)
        else:
            mats = (fundi, fundis)

        self._saved = rate, (peq, c_s), (eta, theta), mats
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

        return scalarise(func), grad

    def laplace_grad(self, rate: Number,
                     inv: bool = False) -> (float, la.lnarray):
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
        rows, cols, mats = self._derivs(rate, inv)
        rcond = 1. / np.linalg.cond(mats[:2]).max()
        if rcond < self.RCondThresh:
            func = np.array(np.inf)
        else:
            func = - rows[1] @ self.weight
        # afunc = -rows[0] @ self.enc() @ cols[0]

        dadq = - _diagsub(rows[0].c * cols[0])
        dadw = - _diagsub(rows[0].c * cols[1] + rows[1].c * cols[0])
        dadwp = self.frac[0] * (dadw + dadq)
        dadwm = self.frac[1] * (dadw - dadq)
        grad = [mp.mat_to_params(m, drn=d, **self.type)
                for m, d in zip([dadwp, dadwm], self.directions)]
        grad = np.hstack(grad)

        return scalarise(func), grad

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
            extra arguments passed to `cls.build` or `builders.build_empty`.

        Returns
        -------
        synobj
            SynapseOpt instance
        """
        mat_type = cls.type.copy()
        del mat_type['uniform']
        mat_type['drn'] = cls.directions[0]
        npl = kwds.pop('npl', len(cls.directions))
        nst = mp.num_state(params.size // npl, **mat_type)
        self = cls.zero(nst, *args, npl=npl, **kwds)
        self.set_params(params)
        return self

    @classmethod
    def rand(cls, nst, *args, **kwds) -> SynapseOptModel:
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
            SynapseOpt instance
        """
        mat_type = cls.type.copy()
        del mat_type['uniform']
        npl = kwds.pop('npl', len(cls.directions))
        npr = npl * mp.num_param(nst, drn=cls.directions[0], **cls.type)
        self = cls.zero(nst, *args, npl=npl, **kwds)
        self.set_params(la.random.rand(npr))
        return self


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
    return la.stack((_trnsp4(tens_trn), tens_trn))


if __name__ == "__main__":
    pass
