"""Baum-Welch update
"""
from __future__ import annotations

from typing import List, Tuple, Union, Sequence

import numpy as np
import matplotlib as mpl

import numpy_linalg as la

from sl_py_tools.iter_tricks import zenumerate, rzenumerate

from ..synapse_base import ArrayLike
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

    # (E,1,M),(E,1,1)
    alpha[0] = initial.r
    eta[0] = 1. / alpha[0].sum(-1, keepdims=True)
    alpha[0] *= eta[0]
    for i, updater in zenumerate(updaters):
        alpha[i+1] = alpha[i] @ updater
        eta[i+1] = 1. / alpha[i+1].sum(-1, keepdims=True)
        alpha[i+1] *= eta[i+1]
    # (E,M,1)
    beta[-1] = 1.
    for i, updater in rzenumerate(updaters):
        beta[i] = updater @ beta[i+1] * eta[i+1]
    return alpha.ur, beta.uc, eta.us


def update_model(model: SynapseIdModel, plast_seq: PlasticitySequence,
                 normalise: bool = True, steady: bool = True) -> Tuple[la.lnarray, la.lnarray]:
    """One Baum-Welch update of the model

    Parameters
    ----------
    model : SynapseIdModel
        The model to update, modified in-place
    plast_seq : PlasticitySequence
        The data to use for the update

    Returns
    -------
    state_probs : la.lnarray, (T,E,M)
        Current estimate of marginal occupation
    nlog_like : la.lnarray, (E,)
        Negative log likelihood of data given old model
    """
    plast_seq = plast_seq.moveaxis(plast_seq.t_axis, 0)
    # (R,P,M,M),(R,M)
    update, initial = model.updaters()
    # (T,E,M,M) - makes a copy :(
    update = update[plast_seq.readouts[1:], plast_seq.plast_type[:-1]]
    # (E,M)
    initial = initial[plast_seq.readouts[0]]
    # (T,E,M),(T,E,M),(T,E),
    alpha, beta, eta = calc_alpha_beta(update, initial)
    # (T,E,M)
    state_prob = alpha * beta
    # (E,)
    nlog_like = np.log(eta).sum(0)

    # (T,E,M,M)
    update *= alpha.c[:-1] * beta.r[1:] * eta.s[1:]
    # (P,M,M)
    model.plast = la.array([update[plast_seq.plast_type[:-1] == i].sum(0)
                            for i in range(model.nplast)])
    # (M,)
    model.initial = state_prob.sum((0, 1)) if steady else state_prob[0].sum(0)
    if normalise:
        model.normalise()
        model.sort(group=True)
        nlog_like = nlog_like.sum()

    return state_prob, nlog_like


def likelihood(model: SynapseIdModel, plast_seq: PlasticitySequence) -> float:
    """Likelihood ov observed readouts given model

    Parameters
    ----------
    model : SynapseIdModel
        The model to update, modified in-place
    plast_seq : PlasticitySequence
        The data to use for the update

    Returns
    -------
    nlog_like : float
        `-log(P(readouts|model))`
    """
    plast_seq = plast_seq.moveaxis(plast_seq.t_axis, 0)
    # (R,P,M,M),(R,M)
    update, initial = model.updaters()
    # (T,E,M,M) - makes a copy :(
    update = update[plast_seq.readouts[1:], plast_seq.plast_type[:-1]]
    # (E,M)
    initial = initial[plast_seq.readouts[0]]
    # (T,E,M),(T,E,M),(T,E),
    return np.log(calc_alpha_beta(update, initial)[2]).sum()
