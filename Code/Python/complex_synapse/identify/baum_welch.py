"""Baum-Welch update
"""
from __future__ import annotations

from typing import List, Tuple, Union, Sequence

import numpy as np
import matplotlib as mpl

import numpy_linalg as la

from sl_py_tools.iter_tricks import zenumerate, rzenumerate

from .synapse_id import SynapseIdModel
from .plast_seq import PlasticitySequence
# =============================================================================


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
    return np.log(calc_alpha_beta(update, initial)[2]).sum()


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

    # (T,E,M)
    state_prob = alpha * beta
    # (M,)
    initial = state_prob.sum((0, 1)) if steady else state_prob[0].sum(0)
    # (T-1,E,M,M)
    update *= alpha.c[:-1] * beta.r[1:] * eta.s[1:]
    # (P,M,M)
    plast = la.array([update[plast_type == i].sum(0)
                      for i in range(model.nplast)])
    return plast, initial


def update_model(model: SynapseIdModel, expt: PlasticitySequence,
                 normalise: bool = True, steady: bool = True
                 ) -> Tuple[la.lnarray, la.lnarray]:
    """One Baum-Welch/Rabiner-Juang update of the model

    Parameters
    ----------
    model : SynapseIdModel
        The model to update, modified in-place
    expt : PlasticitySequence
        The data to use for the update
    normalise : bool, optional
        Should we normalise the result? By default `True`.
    steady : bool, optional
        Can we assume that `initial` is the steady-state?
        If `True`, we can use all times to estimate `initial`.
        If `False`, we can only use `t=0`. By default `True`.

    Returns
    -------
    state_probs : la.lnarray, (T,E,M)
        Current estimate of marginal occupation
    nlog_like : la.lnarray, (E,) or ()
        Negative log likelihood of data given old model.
        If `normalise`: sum over experiments.
    """
    # (T-1,E,M,M),(E,M),(T-1,E) - makes a copy :)
    update, initial, plast_type = get_updaters(model, expt)
    # (T,E,M),(T,E,M),(T,E),
    alpha, beta, eta = calc_alpha_beta(update, initial)

    # (T,E,M)
    state_prob = alpha * beta
    # (E,)
    nlog_like = np.log(eta).sum(0)
    # (M,)
    model.initial = state_prob.sum((0, 1)) if steady else state_prob[0].sum(0)
    # (T-1,E,M,M)
    update *= alpha.c[:-1] * beta.r[1:] * eta.s[1:]
    # (P,M,M)
    model.plast = la.array([update[plast_type == i].sum(0)
                            for i in range(model.nplast)])
    if normalise:
        model.normalise()
        model.sort(group=True)
        nlog_like = nlog_like.sum()

    return state_prob, nlog_like
