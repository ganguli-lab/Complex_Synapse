"""Baum-Welch update
"""
from __future__ import annotations

from typing import ClassVar, Dict, Optional, Tuple

import numpy as np

import numpy_linalg as la
import sl_py_tools.arg_tricks as ag
import sl_py_tools.iter_tricks as it

from .. import synapse_base as _sb
from . import plast_seq as ps
from . import synapse_id as si
from . import fit_synapse as _fs

# =============================================================================


def likelihood(model: si.SynapseIdModel, data: ps.PlasticitySequence) -> float:
    """Negative log-likelihood of observed readouts given model

    Parameters
    ----------
    model : SynapseIdModel
        The current model estimate
    data : PlasticitySequence
        The data to use for the update

    Returns
    -------
    nlike : float
        `-log(P(readouts|model))`
    """
    # (T,E,M,M),(E,M) - makes a copy :(
    update, initial = _get_updaters(model, data)[:2]
    # _,_,(T,E)->()
    return _sb.scalarise(np.log(_calc_alpha_beta(update, initial)[2]).sum())


def state_est(model: si.SynapseIdModel, data: ps.PlasticitySequence
              ) -> la.lnarray:
    """Marginal probability of state occupation at each time

    Parameters
    ----------
    model : SynapseIdModel
        The current model estimate
    data : PlasticitySequence
        The data to use for the update

    Returns
    -------
    state_probs : la.lnarray, (T,E,M)
        Current estimate of marginal occupation
    """
    # (T,E,M,M),(E,M) - makes a copy :(
    update, initial = _get_updaters(model, data)[:2]
    # (T,E,M),(T,E,M),_,
    alpha, beta = _calc_alpha_beta(update, initial)[:2]
    # (T,E,M)
    return alpha * beta


def model_est(model: si.SynapseIdModel, data: ps.PlasticitySequence,
              steady: bool = True, combine: bool = True) -> la.lnarray:
    """New, unnormalised estimate of the model

    Parameters
    ----------
    model : SynapseIdModel
        The current model estimate
    data : PlasticitySequence
        The data to use for the update
    steady : bool, optional
        Can we assume that `initial` is the steady-state?
        If `True`, we can use all times to estimate `initial`.
        If `False`, we can only use `t=0`. By default `True`.
    combine: bool
        Should we sum over experiments? By default `True`.

    Returns
    -------
    plast : array_like, (P,M,M), float[0:1]
        new estimate of transition probability matrix (unnormalised).
    initial : array_like, (M,) float[0:1]
        new estimate of distribution of initial state (unnormalised).
    """
    # (T-1,E,M,M),(E,M),(T-1,E) - makes a copy :)
    update, initial, plast_type = _get_updaters(model, data)
    # (T,E,M),(T,E,M),(T,E),
    alpha, beta, eta = _calc_alpha_beta(update, initial)

    # (P,M,M),(M,)
    return _calc_model(update, plast_type, alpha, beta, eta,
                       steady=steady, combine=combine, nplast=model.nplast)


# =============================================================================
# Helpers
# =============================================================================


def _get_updaters(model: si.SynapseIdModel, data: ps.PlasticitySequence,
                  ) -> Tuple[la.lnarray, la.lnarray, la.lnarray]:
    """Get updater matrices for each time-step

    Parameters
    ----------
    model : SynapseIdModel
        The model to update, modified in-place
    data : PlasticitySequence
        The data to use for the update

    Returns
    -------
    updaters : la.lnarray, (T-1,E,M,M)
        Plasticity matrices multiplied by readout indicators of 'to' state
    initial : la.lnarray, (E,M)
        Initial state distribution multiplied by readout indicators of state
    plast_type : la.lnarray, (T-1,E), int[0:P]
        id of plasticity type after each time-step
    """
    data = data.moveaxis(data.t_axis, 0)
    # (R,P,M,M),(R,M)
    updaters, initial = model.updaters()
    # (T-1,E,M,M) - makes a copy :(
    updaters = updaters[data.readouts[1:], data.plast_type]
    # (E,M)
    initial = initial[data.readouts[0]]
    # (T-1,E,M,M),(E,M),(T-1,E),
    return updaters, initial, data.plast_type


def _calc_alpha_beta(updaters: la.lnarray, initial: la.lnarray
                     ) -> Tuple[la.lnarray, la.lnarray, la.lnarray]:
    """Calculate BW forward/backward variables, (T,E,M),(T,E,M),(T,E),

    Parameters
    ----------
    updaters : la.lnarray, (T-1,E,M,M)
        Plasticity matrices multiplied by readout indicators of 'to' state
    initial : la.lnarray, (E,M)
        Initial state distribution multiplied by readout indicators of state

    Returns
    -------
    alpha : la.lnarray, (T,E,M)
        Normalised Baum-Welch forward variable
    beta : la.lnarray, (T,E,M)
        Scaled Baum-Welch backward variable
    eta : la.lnarray, (T,E)
        Norm of Baum-Welch forward variable
    """
    siz = (updaters.shape[0] + 1,) + updaters.shape[1:-1]
    # (T,E,1,M),(T,E,M,1),(T,E,1,1),
    alpha, beta, eta = la.empty(siz).r, la.empty(siz).c, la.empty(siz[:-1]).s

    def norm(ind: int):
        eta[ind] = 1. / alpha[ind].sum(-1, keepdims=True)
        alpha[ind] *= eta[ind]

    # (E,1,M),(E,1,1)
    alpha[0] = initial.r
    norm(0)
    for i, updater in it.zenumerate(updaters):
        alpha[i+1] = alpha[i] @ updater
        norm(i+1)
    # (E,M,1)
    beta[-1] = 1.
    for i, updater in it.rzenumerate(updaters):
        beta[i] = updater @ beta[i+1] * eta[i+1]
    return alpha.ur, beta.uc, eta.us


def _calc_model(updaters: la.lnarray, plast_type: la.lnarray,
                alpha: la.lnarray, beta: la.lnarray, eta: la.lnarray, *,
                steady: bool = True, combine: bool = True, normed: bool = True,
                nplast: Optional[int] = None) -> Tuple[la.lnarray, la.lnarray]:
    """One Baum-Welch/Rabiner-Juang update of the model

    Parameters
    ----------
    updaters : la.lnarray, (T-1,E,M,M), Modified
        Plasticity matrices multiplied by readout indicators of 'to' state.
    plast_type : la.lnarray, (T-1,E), int[0:P]
        id of plasticity type after each time-step
    alpha : la.lnarray, (T,E,M)
        Normalised Baum-Welch forward variable
    beta : la.lnarray, (T,E,M)
        Scaled Baum-Welch backward variable
    eta : la.lnarray, (T,E)
        Norm of Baum-Welch forward variable

    Keyword only:

    steady : bool, optional
        Can we assume that `initial` is the steady-state?
        If `True`, we can use all times to estimate `initial`.
        If `False`, we can only use `t=0`. By default `True`.
    combine: bool
        Should we sum over experiments? By default `True`.
    normed: ClassVar[bool] = True
        Should we normed the result? By default `True`.
    nplast : imt, optional
        number of plasticity types, P. If `None` calculate from `plast_type`.
        By default `True`.

    Returns
    -------
    plast : array_like, (P,M,M), float[0:1]
        new estimate of transition probability matrix.
    initial : array_like, (M,) float[0:1]
        new estimate of distribution of initial state.
    """
    nplast = ag.default_eval(nplast, plast_type.max)
    nexpt = plast_type.shape[1:]
    if not combine:
        # (E,P,M,M)
        plast = np.empty(nexpt + updaters.shape[-3:])
        # (E,M)
        initial = np.empty(nexpt + updaters.shape[-1:])
        for i in np.ndindex(*nexpt):
            j = np.s_[:,] + i
            plast[i], initial[i] = _calc_model(
                updaters[j], plast_type[j], alpha[j], beta[j], eta[j],
                steady=steady, combine=False, normed=normed, nplast=nplast)

    # (T-1,E,M,M)
    updaters *= alpha.c[:-1] * beta.r[1:] * eta.s[1:]
    # (P,M,M)
    plast = la.array([updaters[plast_type == i].sum(0) for i in range(nplast)])
    # (M,)
    if steady:
        initial = (alpha * beta).sum(tuple(range(len(nexpt) + 1)))
    else:
        initial = (alpha[0] * beta[0]).sum(tuple(range(len(nexpt))))

    if normed:
        plast /= plast.sum(axis=-1, keepdims=True)
        initial /= initial.sum(axis=-1, keepdims=True)

    return plast, initial


# =============================================================================
# Model fitter classes
# =============================================================================


class BaumWelchFitter(_fs.SynapseFitter):
    """Class that performs Baum-Welch/Rabiner-Juang EM updates.

    Parameters
    ----------
    data : SimPlasticitySequence
        The simulated data used to fit the synapse model.
    est : SynapseIdModel
        The initial guess/current estimate of the synapse model.
    Other keywords added to `self.opt` (see `SynapseFitter`).

    Attributes
    ----------
    alpha, beta, eta : la.lnarray
        Normalised Baum-Welch forward/backward variables and the normalisers.
    bw_opt : ClassVar[Dict[str, bool]]
        Optins for BW update
    See `SynapseFitter` for other attributes.

    Options
    -------
    steady : ClassVar[bool] = True
        Can we assume that `initial` is the steady-state?
        If `True`, we can use all times to estimate `initial`.
        If `False`, we can only use `t=0`. By default `True`.
    combine : ClassVar[bool] = True
        Should we sum over experiments? By default `True`.
    normed : ClassVar[bool] = True
        Should we normalise the result? By default `True`.
    All of the above are stored in `BaumWelchFitter.bw_opt`.
    See `SynapseFitter` for other options and attributes.

    See Also
    --------
    SynapseFitter
    """
    # Baum-Welch forward/backward variables
    # (T,E,M)
    alpha: la.lnarray
    # (T,E,M)
    beta: la.lnarray
    # (T,E)
    eta: la.lnarray
    # BW update options
    bw_opt: ClassVar[Dict[str, bool]] = {'steady': True,
                                         'combine': True,
                                         'normed': True}


    def __init__(self, data: ps.PlasticitySequence, est: si.SynapseIdModel,
                 **kwds) -> None:
        super().__init__(data, est, **kwds)
        # (T-1,E,M,M),(E,M),(T-1,E) - makes a copy :)
        update, initial, _ = _get_updaters(self.est, self.data)
        # (T,E,M),(T,E,M),(T,E),
        self.alpha, self.beta, self.eta = _calc_alpha_beta(update, initial)
        self.info['nlike'] = _sb.scalarise(np.log(self.eta).sum())

    def update_info(self) -> None:
        """Calculate stats for termination and display.
        """
        super().update_info()
        nlike = _sb.scalarise(np.log(self.eta).sum())
        self.info['dlike'] = self.info['nlike'] - nlike
        self.info['nlike'] = nlike

    def update_fit(self) -> None:
        """Perform a single update of the model"""
        # (T-1,E,M,M),(E,M),(T-1,E) - makes a copy :)
        update, initial, plast_type = _get_updaters(self.est, self.data)
        # (T,E,M),(T,E,M),(T,E),
        self.alpha, self.beta, self.eta = _calc_alpha_beta(update, initial)

        # (P,M,M),(M,)
        self.est.plast, self.est.initial = _calc_model(
            update, plast_type, self.alpha, self.beta, self.eta,
            nplast=self.est.nplast, **self.bw_opt)

        if self.bw_opt['normed']:
            self.est.sort(group=True)

    def plot_occ(self, handle: ps.Handle, ind: ps.Inds, **kwds) -> ps.Plot:
        """Plot current estimate of state occupation

        Parameters
        ----------
        handle : Union[Axes, Image, Line]
            Axes to plot on, or Image/Lines to update with new data
        ind : Tuple[Union[int, slice], ...]
            Time, experiment indices/slices to plot

        Returns
        -------
        imh : Union[Image, Line]
            Image/Line objects for the plots
        """
        trn = kwds.pop('trn', False)
        # (T.M)
        state_prob = self.alpha[ind] * self.beta[ind]
        state_prob = state_prob if trn else state_prob.T
        return ps.set_plot(handle, state_prob, **kwds)


# =============================================================================
# Fitter with ground truth
# =============================================================================


class GroundedBWFitter(_fs.GroundedFitter, BaumWelchFitter):
    """Class that performs BW/RJ EM updates with known ground-truth.

    Parameters
    ----------
    data : SimPlasticitySequence
        The simulated data used to fit the synapse model.
    est : SynapseIdModel
        The initial guess/current estimate of the synapse model.
    truth : SynapseIdModel
        The model used to generate `data`.
    Other keywords added to `self.opt` (see `SynapseFitter`).

    Attributes
    ----------
    alpha, beta, eta : la.lnarray
        Normalised Baum-Welch forward/backward variables and the normalisers.
    bw_opt : ClassVar[Dict[str, bool]]
        Options for BW update
    See `SynapseFitter` for other attributes.

    Statistics
    ----------
    true_like : float
        Negative log-likelihood of `data` given `truth`.
    true_dmodel : float
        Distance between `truth` and `model`.
    true_dlike : float
        `nlike - true_like`.
    All of the above are stored in `self.info`.
    See `SynapseFitter` for other statistics.

    Options
    -------
    steady : ClassVar[bool] = True
        Can we assume that `initial` is the steady-state?
        If `True`, we can use all times to estimate `initial`.
        If `False`, we can only use `t=0`. By default `True`.
    combine : ClassVar[bool] = True
        Should we sum over experiments? By default `True`.
    normed : ClassVar[bool] = True
        Should we normalise the result? By default `True`.
    All of the above are stored in `BaumWelchFitter.bw_opt`.
    See `SynapseFitter` for other options.

    See Also
    --------
    SynapseFitter, GroundedFitter, BaumWelchFitter.
    """

    def __init__(self, data: ps.SimPlasticitySequence, est: si.SynapseIdModel,
                 truth: si.SynapseIdModel, **kwds) -> None:
        super().__init__(data, est, truth=truth, **kwds)
        if 'truth' in self.info:
            self.truth = self.info.pop('truth')
        self.info['true_like'] = likelihood(self.truth, self.data)
        self.info['y_scale'] = self.info['true_like']
        self.info['true_dlike'] = self.info['nlike'] - self.info['true_like']

    def update_info(self) -> None:
        """Calculate stats for termination and display.
        """
        super().update_info()
        self.info['true_dlike'] = self.info['nlike'] - self.info['true_like']
