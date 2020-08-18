"""Baum-Welch update
"""
from __future__ import annotations

from typing import ClassVar, Dict, Optional, Tuple, Union

import numpy as np

import numpy_linalg as la
from sl_py_tools.arg_tricks import default_eval
from sl_py_tools.iter_tricks import rzenumerate, zenumerate

from ..builders import scalarise
from .fit_synapse import GroundedFitter, SynapseFitter
from .plast_seq import (Axes, Image, Line, PlasticitySequence,
                        SimPlasticitySequence, set_plot)
from .synapse_id import SynapseIdModel
# =============================================================================


def likelihood(model: SynapseIdModel, expt: PlasticitySequence) -> float:
    """Likelihood of observed readouts given model

    Parameters
    ----------
    model : SynapseIdModel
        The current model estimate
    expt : PlasticitySequence
        The data to use for the update

    Returns
    -------
    nlog_like : float
        `-log(P(readouts|model))`
    """
    # (T,E,M,M),(E,M) - makes a copy :(
    update, initial = get_updaters(model, expt)[:2]
    # _,_,(T,E)->()
    return scalarise(np.log(calc_alpha_beta(update, initial)[2]).sum())


def state_est(model: SynapseIdModel, expt: PlasticitySequence) -> la.lnarray:
    """Marginal probability of state occupation at each time

    Parameters
    ----------
    model : SynapseIdModel
        The current model estimate
    expt : PlasticitySequence
        The data to use for the update

    Returns
    -------
    state_probs : la.lnarray, (T,E,M)
        Current estimate of marginal occupation
    """
    # (T,E,M,M),(E,M) - makes a copy :(
    update, initial = get_updaters(model, expt)[:2]
    # (T,E,M),(T,E,M),_,
    alpha, beta = calc_alpha_beta(update, initial)[:2]
    # (T,E,M)
    return alpha * beta


def model_est(model: SynapseIdModel, expt: PlasticitySequence,
              steady: bool = True) -> la.lnarray:
    """New, unnormalised estimate of the model

    Parameters
    ----------
    model : SynapseIdModel
        The current model estimate
    expt : PlasticitySequence
        The data to use for the update

    Returns
    -------
    plast : array_like, (P,M,M), float[0:1]
        new estimate of transition probability matrix (unnormalised).
    initial : array_like, (M,) float[0:1]
        new estimate of distribution of initial state (unnormalised).
    """
    # (T-1,E,M,M),(E,M),(T-1,E) - makes a copy :)
    update, initial, plast_type = get_updaters(model, expt)
    # (T,E,M),(T,E,M),(T,E),
    alpha, beta, eta = calc_alpha_beta(update, initial)

    # (P,M,M),(M,)
    return calc_model(update, plast_type, alpha, beta, eta,
                      steady=steady, nplast=model.nplast)


# =============================================================================
# Helpers
# =============================================================================


def get_updaters(model: SynapseIdModel, expt: PlasticitySequence,
                 ) -> Tuple[la.lnarray, la.lnarray, la.lnarray]:
    """Get updater matrices for each time-step

    Parameters
    ----------
    model : SynapseIdModel
        The model to update, modified in-place
    expt : PlasticitySequence
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
    expt = expt.moveaxis(expt.t_axis, 0)
    # (R,P,M,M),(R,M)
    updaters, initial = model.updaters()
    # (T-1,E,M,M) - makes a copy :(
    updaters = updaters[expt.readouts[1:], expt.plast_type]
    # (E,M)
    initial = initial[expt.readouts[0]]
    # (T-1,E,M,M),(E,M),(T-1,E),
    return updaters, initial, expt.plast_type


def calc_alpha_beta(updaters: la.lnarray, initial: la.lnarray
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
    for i, updater in zenumerate(updaters):
        alpha[i+1] = alpha[i] @ updater
        norm(i+1)
    # (E,M,1)
    beta[-1] = 1.
    for i, updater in rzenumerate(updaters):
        beta[i] = updater @ beta[i+1] * eta[i+1]
    return alpha.ur, beta.uc, eta.us


def calc_model(updaters: la.lnarray, plast_type: la.lnarray,
               alpha: la.lnarray, beta: la.lnarray, eta: la.lnarray, *,
               steady: bool = True, combine: bool = True, normed: bool = True,
               nplast: Optional[int] = None
               ) -> Tuple[la.lnarray, la.lnarray]:
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
    nplast = default_eval(nplast, plast_type.max)
    nexpt = plast_type.shape[1:]
    if not combine:
        # (E,P,M,M)
        plast = np.empty(nexpt + updaters.shape[-3:])
        # (E,M)
        initial = np.empty(nexpt + updaters.shape[-1:])
        for i in np.ndindex(*nexpt):
            j = np.s_[:,] + i
            plast[i], initial[i] = calc_model(
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


class BaumWelchFitter(SynapseFitter):
    """Class that performs Baum-Welch/Rabiner-Juang EM updates.

    Parameters
    ----------
    data : SimPlasticitySequence
        The simulated data used to fit the synapse model.
    model : SynapseIdModel
        The initial guess/current estimate of the synapse model.
    Other keywords added to `self.stats` (see `SynapseFitter`).

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


    def __init__(self, data: PlasticitySequence, model: SynapseIdModel, **kwds
                 ) -> None:
        super().__init__(data, model, **kwds)
        # (T-1,E,M,M),(E,M),(T-1,E) - makes a copy :)
        update, initial, _ = get_updaters(self.model, self.data)
        # (T,E,M),(T,E,M),(T,E),
        self.alpha, self.beta, self.eta = calc_alpha_beta(update, initial)
        self.stats['nlog_like'] = scalarise(np.log(self.eta).sum())

    def update_stats(self) -> None:
        """Calculate stats for termination and display.
        """
        super().update_stats()
        nlog_like = scalarise(np.log(self.eta).sum())
        self.stats['dnlog_like'] = self.stats['nlog_like'] - nlog_like
        self.stats['nlog_like'] = nlog_like

    def update_model(self) -> None:
        """Perform a single update of the model"""
        # (T-1,E,M,M),(E,M),(T-1,E) - makes a copy :)
        update, initial, plast_type = get_updaters(self.model, self.data)
        # (T,E,M),(T,E,M),(T,E),
        self.alpha, self.beta, self.eta = calc_alpha_beta(update, initial)

        # (P,M,M),(M,)
        self.model.plast, self.model.initial = calc_model(
            update, plast_type, self.alpha, self.beta, self.eta,
            nplast=self.model.nplast, **self.bw_opt)

        if self.bw_opt['normed']:
            self.model.sort(group=True)

    def plot_occ(self, handle: Union[Axes, Image, Line],
                 ind: Tuple[Union[int, slice], ...],
                 **kwds) -> Union[Image, Line]:
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
        # (T.M)
        state_prob = self.alpha[ind] * self.beta[ind]
        return set_plot(handle, state_prob.T, **kwds)


# =============================================================================
# Fitter with ground truth
# =============================================================================


class GroundedBWFitter(GroundedFitter, BaumWelchFitter):
    """Class that performs BW/RJ EM updates with known ground-truth.

    Parameters
    ----------
    data : SimPlasticitySequence
        The simulated data used to fit the synapse model.
    model : SynapseIdModel
        The initial guess/current estimate of the synapse model.
    truth : SynapseIdModel
        The model used to generate `data`.
    Other keywords added to `self.stats` (see `SynapseFitter`).

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
        `nlog_like - true_like`.
    All of the above are stored in `self.stats`.
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

    def __init__(self, data: SimPlasticitySequence, model: SynapseIdModel,
                 truth: SynapseIdModel, **kwds) -> None:
        super().__init__(data, model, truth=truth, **kwds)
        if 'truth' in self.stats:
            self.truth = self.stats.pop('truth')
        self.stats['true_like'] = likelihood(self.truth, self.data)
        self.stats['y_scale'] = self.stats['true_like']
        self.stats['true_dlike'] = (self.stats['nlog_like']
                                    - self.stats['true_like'])

    def update_stats(self) -> None:
        """Calculate stats for termination and display.
        """
        super().update_stats()
        self.stats['true_dlike'] = (self.stats['nlog_like']
                                    - self.stats['true_like'])

    # def update_model(self) -> None:
    #     """Perform a single update of the model"""
    #     super().update_model()

    # def plot_occ(self, handle: Union[Axes, Image, Line],
    #              ind: Tuple[Union[int, slice], ...],
    #              **kwds) -> Union[Image, Line]:
    #     """Plot current estimate of state occupation

    #     Parameters
    #     ----------
    #     handle : Union[Axes, Image, Line]
    #         Axes to plot on, or Image/Lines to update with new data
    #     ind : Tuple[Union[int, slice], ...]
    #         Time, experiment indices/slices to plot

    #     Returns
    #     -------
    #     imh : Union[Image, Line]
    #         Image/Line objects for the plots
    #     """
    #     return super().plot_occ(handle, ind, **kwds)
