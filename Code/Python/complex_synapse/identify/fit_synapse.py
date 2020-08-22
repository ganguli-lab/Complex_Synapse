"""Classes for updating synapse models to fit data
"""
from __future__ import annotations

import abc
from dataclasses import dataclass
from numbers import Number
from typing import Callable, Dict, List, Optional

import numpy as np

import sl_py_tools.arg_tricks as ag
import sl_py_tools.iter_tricks as it

from . import plast_seq as ps
from . import synapse_id as si

# =============================================================================
MESSAGES = (
    "Failure, maximum iterations reached.",
    "Success, changes in log-likelihood and model estimate below threshold.",
    "Error, invalid model generated.",
)


def _like_str(model: str = '') -> str:
    """TeX string for likelihood"""
    return f"P(\\text{{{{data}}}}|\\mathbf{{{{M}}}}{model})"


def _frmt_vars(format_spec: str) -> Dict[str, str]:
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
                }
    return {"t": "nit = {nit}",
            "tr": "-log P(data|true) = {true_like:.3f}",
            "fit": "-log P(data|fit) = {nlike:.3f}",
            "trx": "\n||true-fit|| = {true_dmodel:.3g}",
            "dtr": "log[P(data|true) / P(data|fit)] = {true_dlike:.3g}",
            "df": "\nlog[P(data|fit) / P(data|prev)] = {dlike:.3g}",
            "dx": "||fit-prev|| = {dmodel:.3g}",
            }


def _frmt_out(disp: List[str], info: Dict[str, Number], frm_spec: str) -> str:
    """Variables for __format__ method"""
    if frm_spec == "tex":
        pre, delim, post = "\\begin{align*} ", " ,\\\\ ", ". \\end{align*}"
    else:
        pre, delim, post = "", ", ", ""
    return pre + delim.join(disp).format(**info) + post


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
    if pos == 0 and obj.opt.verbose >= 0:
        gap = ('',) if obj.opt.verbose >= 2 else ()
        print('Before:', obj, *gap, sep='\n')
    elif pos == 1 and obj.opt.verbose >= 2:
        if obj.info['nit'] % obj.opt.disp_step == 0:
            print('*', obj)
    elif pos == 2 and obj.opt.verbose >= 0:
        print('', 'After:', obj, MESSAGES[obj.info['result']], sep='\n')


# =============================================================================
# Options class for synapse model fitters
# =============================================================================


@dataclass
class SynapseFitOptions:
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
    verbose : int = -2
        When statistics are printed, and how verbose:
            -2: print called
            -1: print called, detailed
            0: end of iteration
            1: end of iteration, detailed
            2: each iteration
            3: each iteration, detailed
    disp_step : int = -2
        Display progress update every `disp_step` iterations.
    """
    # Absolute tolerance for `dmodel`.
    atolx: float = 1e-5
    # Absolute tolerance for `dlike`.
    atoly: float = 1e-5
    # Relative tolerance for `dmodel`. Multiplies `x_scale` if given.
    rtolx: float = 1e-5
    # Relative tolerance for `dlike`. Multiplies `y_scale` if given.
    rtoly: float = 1e-5
    # Maximum number of iterations
    max_it: int = 1000
    # When statistics are printed, and how verbose:
    verbose: int = -2
    # Display progress update every `disp_step` iterations.
    disp_step: int = -2


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
    verbose : int = -2
        When statistics are printed, and how verbose:
            -2: print called
            -1: print called, detailed
            0: end of iteration
            1: end of iteration, detailed
            2: each iteration
            3: each iteration, detailed
    disp_step : int = -2
        Display progress update every `disp_step` iterations.
    All of the above are stored in `SynapseFitter.opt`.
    """
    data: ps.PlasticitySequence
    est: si.SynapseIdModel
    prev_est: Optional[si.SynapseIdModel]
    # Stats of current state:
    info: Dict[str, Number]
    # options:
    opt: SynapseFitOptions

    def __init__(self, data: ps.PlasticitySequence, est: si.SynapseIdModel,
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
        self.opt = SynapseFitOptions(**kwds)

    def __format__(self, format_spec: str) -> str:
        """Printing info"""
        templates = _frmt_vars(format_spec)
        disp = [templates['t'], templates['fit']]
        if self.opt.verbose % 2 and np.isfinite(self.info['dlike']):
            disp += [templates['df'], templates['dx']]
        elif self.opt.verbose % 2 and not format_spec:
            disp += [templates['df'].replace('{dlike:.2g}', '-'),
                     templates['dx'].replace('{dmodel:.2g}', '-')]
        return _frmt_out(disp, self.info, format_spec)

    def __str__(self) -> str:
        """Printing info"""
        return self.__format__('')

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = type(self).__name__ + "(\n"
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
        x_scale = ag.default_eval(self.info['x_scale'], self.est.norm)
        y_scale = ag.default(self.info['y_scale'], self.info['nlike'])
        self.info['x_thresh'] = self.opt.atolx + self.opt.rtolx * x_scale
        self.info['y_thresh'] = self.opt.atoly + self.opt.rtoly * y_scale

    def check_thresh(self) -> bool:
        """Check if last update was below threshold"""
        return (self.info['dmodel'] < self.info['x_thresh'] and
                abs(self.info['dlike']) < self.info['y_thresh'])

    @abc.abstractmethod
    def update_info(self) -> None:
        """Calculate stats for termination and display.

        Subclass must update `info[nlike]` and `info[dlike]`.
        """
        self.info['dmodel'] = (self.est - self.prev_est).norm()

    @abc.abstractmethod
    def update_fit(self) -> None:
        """Perform a single update of the model"""

    def valid(self) -> bool:
        """Check that current model estimate is valid"""
        return not self.est.nmodel and si.valid_values(self.est)

    def run(self, callback: Callback = print_callback) -> int:
        """Run the synapse fitter until termination conditions are met

        Parameters
        ----------
        callback : Callable[[self, int]->None], optional
            Function called on every iteration, by default `print_callback`.
            Second argument:
                0: Before first iteration.
                1: During itertions.
                2: After completion.

        Returns
        -------
        result : int
            Flag for the outcome:
                -1: Error, invalid model generated.
                0: Failure, maximum iterations reached.
                1: Success, change in log-likelihood and model below threshold.
        """
        count = it.undcount if self.opt.verbose >= 2 else it.dcount
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
    """
    data: ps.SimPlasticitySequence
    truth: si.SynapseIdModel

    def __init__(self, data: ps.SimPlasticitySequence, est: si.SynapseIdModel,
                 truth: si.SynapseIdModel, **kwds) -> None:
        """SynapseFitter where groud-truth is known.

        Parameters
        ----------
        data : SimPlasticitySequence
            The data used to fit the synapse model.
        model : SynapseIdModel
            The initial guess/current estimate of the synapse model
        truth : SynapseIdModel
            The model used to generate `data`.

        Subclass should set `self.info` items `['true_like', 'y_scale',
        'true_dlike']`.
        """
        super().__init__(data, est, **kwds)
        self.truth = truth
        self.info.update(true_dmodel=(self.truth - self.est).norm(),
                         true_like=0., x_scale=self.truth.norm())

    def __format__(self, format_spec: str) -> str:
        """Printing info"""
        templates = _frmt_vars(format_spec)
        disp = [templates[k] for k in ['t', 'tr', 'fit', 'trx']]
        if self.prev_est is not None and not format_spec:
            del disp[1]
        if self.opt.verbose % 2 and np.isfinite(self.info['dlike']):
            disp += [templates[k] for k in ['dtr', 'df', 'dx']]
        elif self.opt.verbose % 2 and format_spec:
            disp += [templates['dtr'],
                     templates['df'].replace('{dlike:.2g}', '-'),
                     templates['dx'].replace('{dmodel:.2g}', '-')]
        return _frmt_out(disp, self.info, format_spec)

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = super().__repr__()
        insert = "    truth = "
        insert += repr(self.data).replace("\n", "\n" + " " * 12) + ",\n"
        ind = rpr.find("    info = ")
        rpr = rpr[:ind] + insert + rpr[ind:]
        return rpr

    @abc.abstractmethod
    def update_info(self) -> None:
        """Calculate info for termination.

        Subclass must update `self.info` fields `[nlike, dlike,
        true_dlike]`.
        """
        super().update_info()
        self.info['true_dmodel'] = (self.truth - self.est).norm()


# =============================================================================
Callback = Callable[[SynapseFitter, int], None]
