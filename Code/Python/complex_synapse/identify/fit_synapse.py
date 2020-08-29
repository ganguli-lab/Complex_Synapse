# -*- coding: utf-8 -*-
"""Classes for updating synapse models to fit data
"""
from __future__ import annotations

import abc
from numbers import Number
import re
from typing import Callable, ClassVar, Dict, Optional, Tuple

import numpy as np

import sl_py_tools.arg_tricks as _ag
import sl_py_tools.iter_tricks as _it

from . import plast_seq as _ps
from . import synapse_id as _si
from .. import options as _opt

# =============================================================================
MESSAGES = (
    "Failure, maximum iterations reached.",
    "Success, changes in log-likelihood and model estimate below threshold.",
    "Error, invalid model generated.",
)
_FSPEC = re.compile(r'(\w*?)(\d?),(\d?)')

def _like_str(model: str = '') -> str:
    """TeX string for likelihood"""
    return f"P(\\text{{{{data}}}}|\\mathbf{{{{M}}}}{model})"


def _frmt_vars(format_spec: str) -> Tuple[Dict[str, str], str, str, str]:
    """Variables for __format__ method"""
    if format_spec == "tex":
        return {"t": "t &= {nit:d}",
                "tr": "-\\log " + _like_str("^*") + " &= {true_like:.0f}",
                "fit": "-\\log " + _like_str() + " &= {nlike:.0f}",
                "trx": ("\\lVert\\mathbf{{M}}^*-\\mathbf{{M}}\\rVert &= "
                        + "{true_dmodel:.2f}"),
                "dtr": ("\\log\\frac{{" + _like_str('^*') + "}}{{"
                        + _like_str() + "}} &= {true_dlike:.0f}"),
                "df": "\\Delta\\log " + _like_str() + " &= {dlike:.2g}",
                "dx": "\\lVert\\Delta\\mathbf{{M}}\\rVert &= {dmodel:.2g}",
                "df-": "\\Delta\\log " + _like_str() + " &= -",
                "dx-": "\\lVert\\Delta\\mathbf{{M}}\\rVert &= -",
                }, "\\begin{align*} ", " ,\\\\ ", ". \\end{align*}"
    return {"t": "nit = {nit}",
            "tr": "-log P(data|true) = {true_like:.3f}",
            "fit": "-log P(data|fit) = {nlike:.3f}",
            "trx": "\n||true-fit|| = {true_dmodel:.3g}",
            "dtr": "log[P(data|true) / P(data|fit)] = {true_dlike:.3g}",
            "df": "\nlog[P(data|fit) / P(data|prev)] = {dlike:.3g}",
            "dx": "||fit-prev|| = {dmodel:.3g}",
            "df-": "\nlog[P(data|fit) / P(data|prev)] = ---",
            "dx-": "||fit-prev|| = ---",
            }, "", ", ", ""


def _get_pos_verbosity(obj: SynapseFitter, format_spec: str):
    """Read pos & verbosity from format_spec if possible else from obj"""
    fspec = _FSPEC.fullmatch(format_spec)
    if fspec is not None:
        format_spec = fspec[1]
    if fspec is None or not fspec[2]:
        # not from a callback, so not during
        pos = 2 * np.isfinite(obj.info['dlike'])
    else:
        pos = int(fspec[2])
    if fspec is None or not fspec[3]:
        verbose = obj.opt.disp_when()[pos] == 2
    else:
        verbose = fspec[3] == '2'
    return format_spec, pos, verbose


def print_callback(obj: SynapseFitter, pos: int) -> None:
    """Callback that prints fitter state as appropriate

    Parameters
    ----------
    obj : SynapseFitter
        Object performing fit whose state we print.
    pos : int
        At what stage of the fit are we?
            0: Before first iteration.
            1: During itertions.
            2: After completion.
    """
    fmt = f"{pos},{obj.opt.disp_when()[pos]}"
    if pos == 0 and obj.opt.disp_before:
        gap = ('',) if obj.opt.disp_each else ()
        print('Before:', format(obj, fmt), *gap, sep='\n')
    elif pos == 1 and obj.opt.disp_each:
        if obj.info['nit'] % obj.opt.disp_step == 0:
            print('*', format(obj, fmt))
    elif pos == 2 and obj.opt.disp_after:
        gap = ('',) if obj.opt.disp_each or obj.opt.disp_before else ()
        print(*gap, 'After:', format(obj, fmt), MESSAGES[obj.info['result']], sep='\n')


# =============================================================================
# Options class for synapse model fitters
# =============================================================================


# pylint: disable=too-many-ancestors
class SynapseFitOptions(_opt.Options):
    """Options for synapse fitters

    Parameters
    ----------
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
    prop_attributes: ClassVar[Tuple[str, ...]] = ('disp_after', 'disp_before',
                                                  'disp_each')
    # Absolute tolerance for `dmodel`.
    atolx: float
    # Absolute tolerance for `dlike`.
    atoly: float
    # Relative tolerance for `dmodel`. Multiplies `x_scale` if given.
    rtolx: float
    # Relative tolerance for `dlike`. Multiplies `y_scale` if given.
    rtoly: float
    # Maximum number of iterations
    max_it: int
    # When statistics are printed, and how verbose:
    verbosity: int
    # Display progress update every `disp_step` iterations.
    disp_step: int

    def __init__(self, **kwds) -> None:
        self.atolx = 1e-5
        self.atoly = 1e-5
        self.rtolx = 1e-5
        self.rtoly = 1e-5
        self.max_it = 1000
        self.verbosity = 1
        self.disp_step = 50
        self.update(kwds)

    def set_disp_after(self, value: int) -> None:
        """Display after the end of iteration?

        0: do not print, 1: print summary, 2: print detailed.
        """
        if not 0 <= value < 3:
            raise ValueError(f"Allowed values: 0,1,2, not {value}.")
        change = value - self.disp_after
        if change:
            self.verbosity += change

    def set_disp_before(self, value: int) -> None:
        """Display before starting iteration?

        0: do not print, 1: print summary, 2: print detailed.
        """
        if not 0 <= value < 3:
            raise ValueError(f"Allowed values: 0,1,2, not {value}.")
        change = value - self.disp_after
        if change:
            self.verbosity += change * 3

    def set_disp_each(self, value: int) -> None:
        """Display at each iteration?

        0: do not print, 1: print summary, 2: print detailed.
        """
        if not 0 <= value < 3:
            raise ValueError(f"Allowed values: 0,1,2, not {value}.")
        change = value - self.disp_after
        if change:
            self.verbosity += change * 9

    @property
    def disp_after(self) -> int:
        """Display after the end of iteration?

        0: do not print, 1: print summary, 2: print detailed.
        """
        return self.verbosity % 3

    @property
    def disp_before(self) -> int:
        """Display before starting iteration?

        0: do not print, 1: print summary, 2: print detailed.
        """
        return (self.verbosity // 3) % 3

    @property
    def disp_each(self) -> int:
        """Display at each iteration?

        0: do not print, 1: print summary, 2: print detailed.
        """
        return (self.verbosity // 9) % 3

    def disp_when(self) -> Tuple[int, int, int]:
        """Tuple of (disp_before, disp_each, disp_after)"""
        return self.disp_before, self.disp_each, self.disp_after


# =============================================================================
# Base class for synapse model fitters
# =============================================================================


class SynapseFitter(abc.ABC):
    """Base class for synapse fitters.

    Parameters
    ----------
    data : PlasticitySequence
        The data used to fit the synapse model.
    est : SynapseIdModel
        The initial guess/current estimate of the synapse model.
    Other keywords added to `self.opts` below.

    Other Attributes
    ----------------
    prev_model : Optional[SynapseIdModel]
        During iteration or after error: Previous model estimate.
        Otherwise: `None`.
    info : Dict[str, Number]
        Statistics of current state. See below.
    opt : SynapseFitOptions
        Options for iteration and termination. See below.

    Statistics
    ----------
    nlike : float
        Negative log-likelihood of `data` given `model`.
    dlike : float
        Decrease in negative log-likelihood of `data` given `model`.
    dmodel : float
        Distance between `model` and `prev_model`.
    nit : int
        Number of iterations.
    x_thresh : float
        Threshold on `dmodel` for termination.
    y_thresh : float
        Threshold on `dlike` for termination.
    x_scale : float or None
        Typical scale for `dmodel`. By default, `None`.
    y_scale : float or None
        Typical scale for `dlike`. By default, `None`.
    result : int
        Flag for the outcome of fitting:
            -1: Error, invalid model generated.
            0: Failure, maximum iterations reached.
            1: Success, change in log-likelihood and model below threshold.
    All of the above are stored in `self.info`.

    Options
    -------
    atolx : float = 1e-5
        Absolute tolerance for `dmodel`.
    atoly : float = 1e-5
        Absolute tolerance for `dlike`.
    rtol : float = 1e-5
        Relative tolerance. Multiplies `x_scale` and `y_scale` if given.
    max_it : int = 1e3
        Maximum number of iterations
    verbosity : int = -2
        When statistics are printed, and how verbose:
            0: do not print
            1: after iteration
            2: after iteration, detailed
            3: before iteration
            6: before iteration, detailed
            9: each iteration
            18: each iteration, detailed
        Values in different categories can be summed to produce combinations
    disp_step : int = -2
        Display progress update every `disp_step` iterations.
    All of the above are stored in `SynapseFitter.opt`.

    See Also
    --------
    SynapseFitOptions.
    """
    data: _ps.PlasticitySequence
    est: _si.SynapseIdModel
    prev_est: Optional[_si.SynapseIdModel]
    # Stats of current state:
    info: Dict[str, Number]
    # options:
    opt: SynapseFitOptions

    def __init__(self, data: _ps.PlasticitySequence, est: _si.SynapseIdModel,
                 **kwds) -> None:
        """Base class for synapse fitters.

        Parameters
        ----------
        data : PlasticitySequence
            The data used to fit the synapse model.
        est : SynapseIdModel
            The initial guess/current estimate of the synapse model

        Subclass should set `self.info['nlike']`.
        """
        self.data = data
        self.est = est
        self.prev_est = None
        self.info = {'nit': 0, 'nlike': np.inf, 'dmodel': np.inf,
                     'dlike': np.inf, 'x_thresh': 0., 'y_thresh': 0.,
                     'x_scale': None, 'y_scale': None}
        self.opt = kwds.pop('opt', SynapseFitOptions())
        self.opt.update(**kwds)

    def __format__(self, format_spec: str) -> str:
        """Printing info about state of fitter

        Parameters
        ----------
        format_spec : str
            If it begins `'tex'`, produces a LaTeX equation environment.
            If it ends with `'<digit>,<digit>'`, the first digit is
            interpreted as `pos` and the second as `verbose`.

        For `pos`:
            0: Before first iteration.
            1: During itertions.
            2: After completion.
        For `verbose`:
            1: print summary,
            2: print detailed.
        """
        format_spec, pos, verbose = _get_pos_verbosity(self, format_spec)
        disp = ['t', 'fit']
        if verbose and pos:
            disp += ['df', 'dx']
        elif verbose and not format_spec:
            disp += ['df-', 'dx-']
        templates, pre, delim, post = _frmt_vars(format_spec)
        disped = [templates[k] for k in disp]
        return pre + delim.join(disped).format(**self.info) + post

    def __str__(self) -> str:
        """Printing info"""
        return self.__format__('')

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = type(self).__name__ + "(\n"
        with np.printoptions(threshold=20, precision=2):
            rpr += "    data = "
            rpr += repr(self.data).replace("\n", "\n" + " " * 11) + ",\n"
            rpr += "    model = "
            rpr += repr(self.est).replace("\n", "\n" + " " * 12) + ",\n"
            rpr += "    info = "
            rpr += repr(self.info).replace(", ", ",\n" + " " * 13) + ",\n"
        rpr += ")"
        return rpr

    def calc_thresh(self) -> None:
        """Calculate thresholds for stopping"""
        x_scale = _ag.default_eval(self.info['x_scale'], self.est.norm)
        y_scale = _ag.default(self.info['y_scale'], self.info['nlike'])
        self.info['x_thresh'] = self.opt.atolx + self.opt.rtolx * x_scale
        self.info['y_thresh'] = self.opt.atoly + self.opt.rtoly * y_scale

    def check_thresh(self) -> bool:
        """Check if last update was below threshold"""
        return (self.info['dmodel'] < self.info['x_thresh'] and
                abs(self.info['dlike']) < self.info['y_thresh'])

    @abc.abstractmethod
    def update_info(self) -> None:
        """Calculate stats for termination and display.

        Subclass must update `self.info['nlike']` and `self.info['dlike']`.
        """
        self.info['dmodel'] = (self.est - self.prev_est).norm()

    @abc.abstractmethod
    def update_fit(self) -> None:
        """Perform a single update of the model"""

    def valid(self) -> bool:
        """Check that current model estimate is valid"""
        return not self.est.nmodel and _si.valid_values(self.est)

    def run(self, callback: Callback = print_callback) -> int:
        """Run the synapse fitter until termination conditions are met

        Parameters
        ----------
        callback : Callable[[self, int]->None], optional
            Function called on every iteration, by default `print_callback`.
            Second argument:
                0: Before first iteration.
                1: During iterations.
                2: After completion.

        Returns
        -------
        result : int
            Flag for the outcome:
                -1: Error, invalid model generated.
                0: Failure, maximum iterations reached.
                1: Success, change in log-likelihood and model below threshold.
        """
        count = _it.undcount if self.opt.disp_each else _it.dcount
        callback(self, 0)
        self.info['result'] = 0
        for i in count('iteration', self.opt.max_it,
                       disp_step=self.opt.disp_step):
            self.info['nit'] = i
            self.prev_est = self.est.copy()
            self.update_fit()
            if not self.valid():
                self.info['result'] = -1
                break
            self.update_info()
            self.calc_thresh()
            callback(self, 1)
            self.prev_est = None
            if self.check_thresh():
                self.info['result'] = 1
                break
        callback(self, 2)
        return self.info['result']

    @abc.abstractmethod
    def plot_occ(self, handle: _ps.Handle, ind: _ps.Inds, **kwds) -> _ps.Plot:
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


# =============================================================================
# Fitter with ground truth
# =============================================================================


class GroundedFitter(SynapseFitter):
    """SynapseFitter where groud-truth is known.

    Parameters
    ----------
    data : SimPlasticitySequence
        The simulated data used to fit the synapse model.
    est : SynapseIdModel
        The initial guess/current estimate of the synapse model.
    truth : SynapseIdModel
        The model used to generate `data`.
    Other keywords added to `self.opt` (see `SynapseFitter`).

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

    See Also
    --------
    SynapseFitter.
    SynapseFitOptions.
    """
    data: _ps.SimPlasticitySequence
    truth: _si.SynapseIdModel

    def __init__(self, data: _ps.SimPlasticitySequence,
                 est: _si.SynapseIdModel, truth: _si.SynapseIdModel,
                 **kwds) -> None:
        """SynapseFitter where groud-truth is known.

        Parameters
        ----------
        data : SimPlasticitySequence
            The data used to fit the synapse model.
        model : SynapseIdModel
            The initial guess/current estimate of the synapse model
        truth : SynapseIdModel
            The model used to generate `data`.

        Subclass should set `self.info['true_like', 'y_scale', 'true_dlike']`.
        """
        super().__init__(data, est, **kwds)
        self.truth = truth
        self.info.update(true_dmodel=(self.truth - self.est).norm(),
                         true_like=0., x_scale=self.truth.norm())

    def __format__(self, format_spec: str) -> str:
        """Printing info about state of fitter

        Parameters
        ----------
        format_spec : str
            If it begins `'tex'`, produces a LaTeX equation environment.
            If it ends with `'<digit>,<digit>'`, the first digit is
            interpreted as `pos` and the second as `verbose`.

        For `pos`:
            0: Before first iteration.
            1: During itertions.
            2: After completion.
        For `verbose`:
            1: print summary,
            2: print detailed.
        """
        format_spec, pos, verbose = _get_pos_verbosity(self, format_spec)
        disp = ['t', 'tr', 'fit', 'trx']
        if not verbose and pos == 1 and not format_spec:
            del disp[1]
        if verbose and pos:
            disp += ['dtr', 'df', 'dx']
        elif verbose and format_spec:
            disp += ['dtr', 'df-', 'dx-']
        templates, pre, delim, post = _frmt_vars(format_spec)
        disped = [templates[k] for k in disp]
        return pre + delim.join(disped).format(**self.info) + post

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = super().__repr__()
        insert = "    truth = "
        with np.printoptions(threshold=20, precision=2):
            insert += repr(self.data).replace("\n", "\n" + " " * 12) + ",\n"
        ind = rpr.find("    info = ")
        rpr = rpr[:ind] + insert + rpr[ind:]
        return rpr

    def calc_thresh(self) -> None:
        """Calculate thresholds for stopping"""
        if self.info['x_thresh'] == 0 or self.info['y_thresh'] == 0:
            super().calc_thresh()

    @abc.abstractmethod
    def update_info(self) -> None:
        """Calculate info for termination.

        Subclass must update `self.info['nlike', 'dlike', 'true_dlike']`.
        """
        super().update_info()
        self.info['true_dmodel'] = (self.truth - self.est).norm()


# =============================================================================
Callback = Callable[[SynapseFitter, int], None]
