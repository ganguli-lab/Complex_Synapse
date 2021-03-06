# -*- coding: utf-8 -*-
"""Baum-Welch update
"""
from __future__ import annotations

from typing import Dict, Tuple, Union

import numpy as np

import numpy_linalg as la
import sl_py_tools.iter_tricks as _it

import complex_synapse.synapse_base as _sb
import complex_synapse.identify.plast_seq as _ps
import complex_synapse.identify.synapse_id as _si
import complex_synapse.identify.fit_synapse as _fs
try:
    import complex_synapse.identify._baum_welch as _bw
except ImportError:
    P_OR_C = True
else:
    P_OR_C = False
# =============================================================================


def likelihood(model: _si.SynapseIdModel, data: _ps.PlasticitySequence
               ) -> float:
    """Negative log-likelihood of observed readouts given model.

    Parameters
    ----------
    model : SynapseIdModel
        The current model estimate.
    data : PlasticitySequence
        The data to use for the update.

    Returns
    -------
    nlike : float
        `-log(P(readouts|model))`.
    """
    # _,_,(E,T)->()
    return _sb.scalarise(np.log(_get_bw_abe(model, data)[2]).sum())


def state_est(model: _si.SynapseIdModel, data: _ps.PlasticitySequence
              ) -> Array:
    """Marginal probability of state occupation at each time.

    Parameters
    ----------
    model : SynapseIdModel
        The current model estimate.
    data : PlasticitySequence
        The data to use for the update.

    Returns
    -------
    state_probs : lnarray, (E,T,M)
        Current estimate of marginal occupation.
    """
    # (E,T,M),(E,T,M),_,
    alpha, beta = _get_bw_abe(model, data)[:2]
    # (E,T,M)
    return alpha * beta


def model_est(model: _si.SynapseIdModel, data: _ps.PlasticitySequence,
              steady: bool = True, combine: bool = True) -> Arrays2:
    """New, unnormalised estimate of the model.

    Parameters
    ----------
    model : SynapseIdModel
        The current model estimate.
    data : PlasticitySequence
        The data to use for the update.
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
    # (P,M,M),(M,)
    return _get_bw_abe_model(model, data, steady=steady, combine=combine)[1]


# =============================================================================
# Helpers
# =============================================================================


def _private() -> None:
    """Hold docs for private functions - decoding names.

    Suffix `p`: uses python loops, `c`: uses C-extension.
    A function called, e.g., `_abc_xyz` takes the arguments listed under `abc`
    below, and returns the tuple `xyz`. With the exception of `obj`, these sets
    of parameters should be passed as tuples, without unpacking.

    Parameters and Returns
    ======================
    obj
    ---
    model : SynapseIdModel
        The model to update, modified in-place
    data : PlasticitySequence
        The data to use for the update

    jmp
    ---
    jumps : lnarray, (R,P,M,M)
        Plasticity matrices multiplied by readout indicators of 'to' state.
    starts : lnarray, (R,M)
        Initial state distribution multiplied by readout indicators of state.
    plast_type : ArrayLike, (E,T-1), int[0:P]
        ID of plasticity type after each time-step.
    readouts : ArrayLike, (E,T), int[0:R]
        ID of readout from synapse at each time-step.

    upd
    ---
    updaters : lnarray, (E,T-1,M,M)
        Plasticity matrices multiplied by readout indicators of 'to' state,
        per experiment, per time.
    initial : lnarray, (E,M)
        Initial state distribution multiplied by readout indicators of state,
        per experiment.
    plast_type : lnarray, (E,T-1), int[0:P]
        ID of plasticity type after each time-step.

    abe
    ---
    bw_vars : Tuple[lnarray, lnarray, lnarray] (E,T,M),(E,T,M),(E,T) float
        Baum-Welch forward/backward variables.

        alpha : lnarray, (E,T,M) float[0:1]
            Normalised Baum-Welch forward variable.
        beta : lnarray, (E,T,M) float
            Scaled Baum-Welch backward variable.
        eta : lnarray, (E,T) float[1:]
            Norm of Baum-Welch forward variable.

    kwds
    ----
    steady : bool, optional
        Can we assume that `initial` is the steady-state?
        If `True`, we can use all times to estimate `initial`.
        If `False`, we can only use `t=0`. By default `True`.
    combine: bool
        Should we sum over experiments? By default `True`.
    normed: ClassVar[bool] = True
        Should we normed the result? By default `True`.
    nplast : int, optional [`_p` only]
        Number of plasticity types, P. If `None` calculate from `plast_type`.
        By default `True`.

    model
    -----
    plast : array_like, (P,M,M), float[0:1]
        New estimate of transition probability matrix.
    initial : array_like, (M,) float[0:1]
        New estimate of distribution of initial state.
    """


def _obj_jmp(model: _si.SynapseIdModel, data: _ps.PlasticitySequence) -> Arrays4:
    """convert obj -> jmp"""
    data = data.move_t_axis(-1)
    # (R,P,M,M),(R,M)
    jumps, starts = model.updaters()
    # err = la.gufuncs.make_errobj("GUfunc reported a floating point error")
    return jumps, starts, data.plast_type, data.readouts


def _jmp_upd(jmp: Arrays4) -> Arrays3:
    """convert jmp -> upd, updater matrices for each time-step"""
    jumps, starts, plast_type, readouts = jmp
    # (E,T-1,M,M)
    updaters = jumps[readouts[..., 1:], plast_type]
    # (E,M)
    initial = starts[readouts[..., 0]]
    # (E,T-1,M,M),(E,M),(E,T-1)
    return updaters, initial, plast_type


def _drop_snd(args: Arrays) -> Arrays:
    """Drop the second argument (index 1)"""
    return args[:1] + args[2:]


def _upd_calc_abe_p(upd: Arrays3) -> Arrays3:
    """Calculate BW forward/backward variables.
    """
    # drop plast_types
    updaters, initial = upd[:2]
    # (E,T,M)
    siz = initial.shape[:-1] + (updaters.shape[-3] + 1,) + updaters.shape[-1:]
    # (E,T,1,M),(E,T,M,1),(E,T,1,1),
    alpha, beta, eta = la.empty(siz).r, la.empty(siz).c, la.empty(siz[:-1]).s

    def inds(ind: int) -> Tuple[Union[int, slice]]:
        return np.s_[..., ind, :, :]

    def norm(index: int):
        ind = inds(index)
        # (E,1,1)
        eta[ind] = 1. / alpha[ind].sum(-1, keepdims=True)
        alpha[ind] *= eta[ind]

    # (T-1,E,M,M)
    updaters = updaters.moveaxis(-3, 0)
    # (E,1,M)
    alpha[inds(0)] = initial.r
    # eta(E,1,1)
    norm(0)
    for i, updater in _it.zenumerate(updaters):
        # (E,1,M) @ (E,M,M) -> (E,1,M)
        alpha[inds(i+1)] = alpha[inds(i)] @ updater
        # eta(E,1,1)
        norm(i+1)
    # (E,M,1)
    beta[inds(-1)] = 1.
    for i, updater in _it.rzenumerate(updaters):
        # (E,M,M) @ (E,M,1) * (E,1,1) -> (E,M,1)
        beta[inds(i)] = updater @ beta[inds(i+1)] * eta[inds(i+1)]
    return alpha.ur, beta.uc, eta.us


def _upd_abe_calc_model_p(upd: Arrays3, abe: Arrays3, **kwds) -> Arrays2:
    """One Baum-Welch/Rabiner-Juang update of the model.
    """
    # drop initial
    updaters, plast_type = _drop_snd(upd)
    nplast = kwds.setdefault('nplast', plast_type.max() + 1)
    if not kwds.pop('combine', True):
        axes = plast_type.shape[:-1]
        # (E,P,M,M)
        plast = np.empty(axes + updaters.shape[-3:])
        # (E,M)
        initial = np.empty(axes + updaters.shape[-1:])
        for i in np.ndindex(*axes):
            plast[i], initial[i] = _upd_abe_calc_model_p(
                [x[i] for x in upd], [y[i] for y in abe],
                combine=False, nplast=nplast, **kwds)

    axis = plast_type.ndim - 1
    # (E,T,M),(E,T,M),(E,T) -> (T,E,M),(T,E,M),(T,E)
    alpha, beta, eta = [x.moveaxis(axis, 0) for x in abe]
    # (E,T-1,M,M)
    updaters *= (alpha.c[:-1] * (beta.r * eta.s)[1:]).moveaxis(0, axis)
    # (P,M,M)
    plast = la.array([updaters[plast_type == i].sum(0) for i in range(nplast)])

    axes = tuple(range(plast_type.ndim))
    # (M,)
    if kwds.pop('steady', True):
        initial = (alpha * beta).sum(axes)
    else:
        initial = (alpha[0] * beta[0]).sum(axes[:-1])

    if kwds.pop('normed', True):
        plast /= plast.sum(axis=-1, keepdims=True)
        initial /= initial.sum(axis=-1, keepdims=True)

    return plast, initial


def _jmp_calc_abe_c(jmp: Arrays4) -> Arrays3:
    """Calculate BW forward/backward variables (using C-extension).
    """
    jumps, starts, plast_type, readouts = jmp
    if (readouts >= starts.shape[-2]).any():
        raise IndexError("readouts out of bounds")
    if (plast_type >= jumps.shape[-3]).any():
        raise IndexError("plast_type out of bounds")
    return _bw.alpha_beta(*jmp)


def _jmp_abe_calc_model_c(jmp: Arrays4, abe: Arrays3, **kwds) -> Arrays2:
    """One Baum-Welch/Rabiner-Juang update of the model (using C-extension).
    """
    # drop starts
    jumps, plast_type, readouts = _drop_snd(jmp)
    if (readouts >= jumps.shape[-4]).any():
        raise IndexError("data.readouts out of bounds")
    if (plast_type >= jumps.shape[-3]).any():
        raise IndexError("data.plast_type out of bounsds")

    nexpt = abe[2].ndim - 1  # eta

    if kwds.pop('steady', True):
        plast, initial = _bw.plast_steady(*jmp, *abe)
    else:
        plast, initial = _bw.plast_init(*jmp, *abe)

    if kwds.pop('combine', True):
        plast = plast.sum(tuple(range(nexpt)))
        initial = initial.sum(tuple(range(nexpt)))

    if kwds.pop('normed', True):
        plast /= plast.sum(axis=-1, keepdims=True)
        initial /= initial.sum(axis=-1, keepdims=True)

    return plast, initial


def _get_bw_abe(model: _si.SynapseIdModel, data: _ps.PlasticitySequence,
                ) -> Arrays3:
    """Calculate BW forward/backward variables

    Parameters
    ----------
    model : SynapseIdModel
        The model to update.
    data : PlasticitySequence
        The data to use for the update.

    Returns
    -------
    alpha : lnarray, (E,T,M)
        Normalised Baum-Welch forward variable.
    beta : lnarray, (E,T,M)
        Scaled Baum-Welch backward variable.
    eta : lnarray, (E,T)
        Norm of Baum-Welch forward variable.
    """
    jmp = _obj_jmp(model, data)
    if P_OR_C:
        upd = _jmp_upd(jmp)
        # drop plast_type
        return _upd_calc_abe_p(upd[:2])
    return _jmp_calc_abe_c(jmp)


def _get_bw_abe_model(model: _si.SynapseIdModel, data: _ps.PlasticitySequence,
                      **kwds) -> Tuple[Arrays3, Arrays2]:
    """One Baum-Welch/Rabiner-Juang update of the model,

    Parameters
    ----------
    model : SynapseIdModel
        The model to update.
    data : PlasticitySequence
        The data to use for the update.

    Keyword only:

    steady : bool, optional
        Can we assume that `initial` is the steady-state?
        If `True`, we can use all times to estimate `initial`.
        If `False`, we can only use `t=0`. By default `True`.
    combine: bool
        Should we sum over experiments? By default `True`.
    normed: ClassVar[bool] = True
        Should we normed the result? By default `True`.

    Returns
    -------
    bw_vars : Tuple[lnarray, lnarray, lnarray] (E,T,M),(E,T,M),(E,T) float
        Baum-Welch forward/backward variables.

        alpha : lnarray, (E,T,M) float[0:1]
            Normalised Baum-Welch forward variable.
        beta : lnarray, (E,T,M) float
            Scaled Baum-Welch backward variable.
        eta : lnarray, (E,T) float[1:]
            Norm of Baum-Welch forward variable.
    model : Tuple[lnarray, lnarray] (P,M,M),(M,) float
        Baum-Welch update.

        plast : array_like, (P,M,M), float[0:1]
            new estimate of transition probability matrix.
        initial : array_like, (M,) float[0:1]
            new estimate of distribution of initial state.
    """
    jmp = _obj_jmp(model, data)
    if P_OR_C:
        upd = _jmp_upd(jmp)
        abe = _upd_calc_abe_p(upd)
        model = _upd_abe_calc_model_p(upd, abe, **kwds)
    else:
        abe = _jmp_calc_abe_c(jmp)
        model = _jmp_abe_calc_model_c(jmp, abe, **kwds)
    return abe, model


# =============================================================================
# Model fitter options class
# =============================================================================


# pylint: disable=too-many-ancestors
class BaumWelchOptions(_fs.SynapseFitOptions):
    """Options for Baum-Welch synapse fitters

    The individual options can be accessed as object instance attributes
    (e.g. `obj.name`) or as dictionary items (e.g. `obj['name']`) for both
    getting and setting.

    Parameters
    ----------
    steady : bool = True
        Can we assume that `initial` is the steady-state?
        If `True`, we can use all times to estimate `initial`.
        If `False`, we can only use `t=0`. By default `True`.
    combine : bool = True
        Should we sum over experiments? By default `True`.
    normed : bool = True
        Should we normalise the result? By default `True`.
    atolx : float = 1e-5
        Absolute tolerance for `dmodel`.
    atoly : float = 1e-5
        Absolute tolerance for `dlike`.
    rtolx : float = 1e-5
        Relative tolerance for `dmodel`. Multiplies `x_scale` if given.
    rtoly : float = 1e-5
        Relative tolerance for `dlike`. Multiplies `y_scale` if given.
    max_it : int = 1e3
        Maximum number of iterations
    verbosity : int = 0
        When statistics are printed, and how verbose:
            0: do not print
            1: after iteration
            2: after iteration, detailed
            3: before iteration
            6: before iteration, detailed
            9: each iteration
            18: each iteration, detailed
        Values in different categories can be summed to produce combinations
    disp_step : int = 50
        Display progress update every `disp_step` iterations.

    All parameters are optional keywords. Any dictionary passed as positional
    parameters will be popped for the relevant items. Keyword parameters must
    be valid keys, otherwise a `KeyError` is raised.

    Properties
    ----------
    disp_before : int
        Display before starting iteration?
    disp_each : int
        Display at each iteration?
    disp_after : int
        Display after the end of iteration?
    They are interpreted as: 0=do not print, 1=print summary, 2=print detailed.
    They are related to `verbose` as `verbose = after + 3 * before + 9 * each`.
    """
    steady: bool
    combine: bool
    normed: bool

    def __init__(self, *args, **kwds) -> None:
        self.steady = True
        self.combine = True
        self.normed = True
        super().__init__(*args, **kwds)

    def bw_opts(self) -> Dict[str, bool]:
        """Those options specifically needed by Baum-Welch fitters
        """
        return {k: self[k] for k in ['steady', 'combine', 'normed']}


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
    alpha, beta, eta : lnarray (E,T,M),(E,T,M),(E,T) float
        Normalised Baum-Welch forward/backward variables and the normalisers.
    opt : BaumWelchOptions
        Optins for BW update.
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
    All of the above are stored in `BaumWelchFitter.opt`.
    See `SynapseFitter` for other options and attributes.

    See Also
    --------
    SynapseFitter.
    SynapseFitOptions, BaumWelchOptions.
    """
    # Baum-Welch forward/backward variables
    # (E,T,M)
    alpha: Array
    # (E,T,M)
    beta: Array
    # (E,T)
    eta: Array
    # BW update options
    opt: BaumWelchOptions


    def __init__(self, data: _ps.PlasticitySequence, est: _si.SynapseIdModel,
                 callback: _fs.Callback = _fs.print_callback, **kwds) -> None:
        kwds.setdefault('opt', BaumWelchOptions())
        super().__init__(data, est, callback, **kwds)
        # (E,T,M),(E,T,M),(E,T)
        self.alpha, self.beta, self.eta = _get_bw_abe(self.est, self.data)
        self.info['nlike'] = _sb.scalarise(np.log(self.eta).sum())

    def update_info(self) -> None:
        """Calculate stats for termination and display.
        """
        super().update_info()
        nlike = _sb.scalarise(np.log(self.eta).sum())
        self.info['dlike'] = self.info['nlike'] - nlike
        self.info['nlike'] = nlike

    def update_fit(self) -> None:
        """Perform a single update of the model."""
        abe, est = _get_bw_abe_model(self.est, self.data, **self.opt.bw_opts())
        # (E,T,M),(E,T,M),(E,T)
        self.alpha, self.beta, self.eta = abe
        # (P,M,M),(M,)
        self.est.plast, self.est.initial = est
        if self.opt.normed:
            self.est.sort(group=True)

    def est_occ(self, ind: _ps.Inds) -> Array:
        """Current estimate of state occupation.

        Parameters
        ----------
        ind : Tuple[Union[int, slice], ...]
            Time, experiment indices/slices to plot.

        Returns
        -------
        data : lnarray, ([E,]T,M) float[0:1]
            Estimate of state occupation.
        """
        # ([E,]T,M)
        return self.alpha[ind] * self.beta[ind]


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
    alpha, beta, eta : lnarray
        Normalised Baum-Welch forward/backward variables and the normalisers.
    opt : BaumWelchOptions
        Options for BW update.
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
    All of the above are stored in `BaumWelchFitter.opt`.
    See `SynapseFitter` for other options.

    See Also
    --------
    SynapseFitter, GroundedFitter, BaumWelchFitter.
    SynapseFitOptions, BaumWelchOptions.
    """

    def __init__(self, data: _ps.SimPlasticitySequence,
                 est: _si.SynapseIdModel, truth: _si.SynapseIdModel,
                 callback: _fs.Callback = _fs.print_callback, **kwds) -> None:
        kwds.setdefault('opt', BaumWelchOptions())
        super().__init__(data, est=est, truth=truth, callback=callback, **kwds)
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


# =============================================================================
# Aliases for hinting
# =============================================================================
Array = la.lnarray
Arrays = Tuple[Array, ...]
Arrays2 = Tuple[Array, Array]
Arrays3 = Tuple[Array, Array, Array]
Arrays4 = Tuple[Array, Array, Array, Array]
