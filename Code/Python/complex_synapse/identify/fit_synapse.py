"""Classes for updating synapse models to fit data
"""
from __future__ import annotations

import abc
from dataclasses import dataclass
from numbers import Number
from typing import Callable, Dict, List, Optional, Tuple, Union

import numpy as np

from sl_py_tools.arg_tricks import default, default_eval
from sl_py_tools.iter_tricks import dcount, undcount

from .plast_seq import (Axes, Image, Line, PlasticitySequence,
                        SimPlasticitySequence)
from .synapse_id import SynapseIdModel, valid_values

# =============================================================================
MESSAGES = {
    -1: "Error, invalid model generated.",
    0: "Failure, maximum iterations reached.",
    1: "Success, changes in log-likelihood and model estimate below threshold."
}


def _like_str(model: str = '') -> str:
    """TeX string for likelihood"""
    return f"P(\\text{{{{data}}}}|\\mathbf{{{{M}}}}{model})"


def _frmt_vars(format_spec: str) -> Dict[str, str]:
    """Variables for __format__ method"""
    if format_spec == "tex":
        return {"t": "t &= {nit:d}",
                "fit": "-\\log " + _like_str() + " &= {nlog_like:.0f}",
                "df": "\\Delta\\log " + _like_str() + " &= {dnlog_like:.2g}",
                "dx": "\\lVert\\Delta\\mathbf{{M}}\\rVert &= {dmodel:.2g}",
                "tr": "-\\log " + _like_str("^*") + " &= {true_like:.0f}",
                "dtr": ("\\log\\frac{{" + _like_str('^*') + "}}{{"
                        + _like_str() + "}} &= {true_dlike:.0f}"),
                "trx": ("\\lVert\\mathbf{{M}}^*-\\mathbf{{M}}\\rVert &= "
                        + "{true_dmodel:.2f}")}
    return {"t": "nit = {nit}",
            "fit": "-log P(data|fit) = {nlog_like:.3f}",
            "df": "\nlog[P(data|fit) / P(data|prev)] = {dnlog_like:.3g}",
            "dx": "||fit-prev|| = {dmodel:.3g}",
            "tr": "-log P(data|true) = {true_like:.3f}",
            "dtr": "\nlog[P(data|true) / P(data|fit)] = {true_dlike:.3g}",
            "trx": "||true-fit|| = {true_dmodel:.3g}"}


def _frmt_out(disp: List[str], stats: Dict[str, Number], frm_spec: str) -> str:
    """Variables for __format__ method"""
    if frm_spec == "tex":
        pre, delim, post = "\\begin{align*} ", " ,\\\\ ", " \\end{align*}"
    else:
        pre, delim, post = "", ", ", ""
    return pre + delim.join(disp).format(**stats) + post


def print_callback(obj: SynapseFitter, pos: int) -> None:
    """Callback that prints fitter state as appropriate"""
    if pos == 0 and obj.opt.verbose >= 0:
        gap = ('',) if obj.opt.verbose >= 2 else ()
        print('Before:', obj, *gap, sep='\n')
    elif pos == 1 and obj.opt.verbose >= 2:
        if obj.stats['nit'] % obj.opt.disp_step == 0:
            print('*', obj)
    elif pos == 2 and obj.opt.verbose >= 0:
        print('', 'After:', obj, MESSAGES[obj.stats['result']], sep='\n')


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
        Absolute tolerance for `dnlog_like`.
    rtolx : float = 1e-5
        Relative tolerance for `dmodel`. Multiplies `x_scale` if given.
    rtoly : float = 1e-5
        Relative tolerance for `dnlog_like`. Multiplies `y_scale` if given.
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
    # Absolute tolerance for `dnlog_like`.
    atoly: float = 1e-5
    # Relative tolerance for `dmodel`. Multiplies `x_scale` if given.
    rtolx: float = 1e-5
    # Relative tolerance for `dnlog_like`. Multiplies `y_scale` if given.
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
    model : SynapseIdModel
        The initial guess/current estimate of the synapse model.
    Other keywords added to `self.opts` below.

    Other Attributes
    ----------------
    prev_model : Optional[SynapseIdModel]
        During itteration or after error: Previous model estimate.
        Otherwise: `None`.
    stats : Dict[str, Number]
        Statistics of current state. See below.
    opt : SynapseFitOptions
        Options for iteration and termination. See below.

    Statistics
    ----------
    nlog_like : float
        Negative log-likelihood of `data` given `model`.
    dnlog_like : float
        Decrease in negative log-likelihood of `data` given `model`.
    dmodel : float
        Distance between `model` and `prev_model`.
    nit : int
        Number of iterations.
    x_thresh : float
        Threshold on `dmodel` for termination.
    y_thresh : float
        Threshold on `dnlog_like` for termination.
    x_scale : float or None
        Typical scale for `dmodel`. By default, `None`.
    y_scale : float or None
        Typical scale for `dnlog_like`. By default, `None`.
    result : int
        Flag for the outcome of fitting:
            -1: Error, invalid model generated.
            0: Failure, maximum iterations reached.
            1: Success, change in log-likelihood and model below threshold.
    All of the above are stored in `self.stats`.

    Options
    -------
    atolx : float = 1e-5
        Absolute tolerance for `dmodel`.
    atoly : float = 1e-5
        Absolute tolerance for `dnlog_like`.
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
    data: PlasticitySequence
    model: SynapseIdModel
    prev_model: Optional[SynapseIdModel]
    # Stats of current state:
    stats: Dict[str, Number]
    # options:
    opt: SynapseFitOptions

    def __init__(self, data: PlasticitySequence, model: SynapseIdModel, **kwds
                 ) -> None:
        """Base class for synapse fitters.

        Parameters
        ----------
        data : SimPlasticitySequence
            The data used to fit the synapse model.
        model : SynapseIdModel
            The initial guess/current estimate of the synapse model

        Subclass should set `self.stats['nlog_like']`.
        """
        self.data = data
        self.model = model
        self.prev_model = None
        self.stats = {'nit': 0, 'nlog_like': np.inf, 'dmodel': np.inf,
                      'dnlog_like': np.inf, 'x_thresh': 0., 'y_thresh': 0.,
                      'x_scale': None, 'y_scale': None}
        self.opt = SynapseFitOptions(**kwds)

    def __format__(self, format_spec: str) -> str:
        """Printing stats"""
        templates = _frmt_vars(format_spec)
        disp = [templates['t'], templates['fit']]
        if self.opt.verbose % 2 and np.isfinite(self.stats['dnlog_like']):
            disp += [templates['df'], templates['dx']]
        elif self.opt.verbose % 2 and not format_spec:
            disp += [templates['df'].replace('{dnlog_like:.2g}', '-'),
                     templates['dx'].replace('{dmodel:.2g}', '-')]
        return _frmt_out(disp, self.stats, format_spec)

    def __str__(self) -> str:
        """Printing stats"""
        return self.__format__('')

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = type(self).__name__ + "(\n"
        rpr += "    data = "
        rpr += repr(self.data).replace("\n", "\n" + " " * 11) + ",\n"
        rpr += "    model = "
        rpr += repr(self.model).replace("\n", "\n" + " " * 12) + ",\n"
        rpr += "    stats = "
        rpr += repr(self.stats).replace(", ", ",\n" + " " * 13) + ",\n"
        rpr += ")"
        return rpr

    def calc_thresh(self) -> None:
        """Calculate thresholds for stopping"""
        x_scale = default_eval(self.stats['x_scale'], self.model.norm)
        y_scale = default(self.stats['y_scale'], self.stats['nlog_like'])
        self.stats['x_thresh'] = self.opt.atolx + self.opt.rtolx * x_scale
        self.stats['y_thresh'] = self.opt.atoly + self.opt.rtoly * y_scale

    def check_thresh(self) -> bool:
        """Check if last update was below threshold"""
        return (self.stats['dmodel'] < self.stats['x_thresh'] and
                abs(self.stats['dnlog_like']) < self.stats['y_thresh'])

    @abc.abstractmethod
    def update_stats(self) -> None:
        """Calculate stats for termination and display.

        Subclass must update `stats[nlog_like]` and `stats[dnlog_like]`.
        """
        self.stats['dmodel'] = (self.model - self.prev_model).norm()

    @abc.abstractmethod
    def update_model(self) -> None:
        """Perform a single update of the model"""

    def valid(self) -> bool:
        """Check that current model estimate is valid"""
        return not self.model.nmodel and valid_values(self.model)

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
        count = undcount if self.opt.verbose >= 2 else dcount
        callback(self, 0)
        self.stats['result'] = 0
        for i in count('iteration', self.opt.max_it,
                       disp_step=self.opt.disp_step):
            self.stats['nit'] = i
            self.prev_model = self.model.copy()
            self.update_model()
            if not self.valid():
                self.stats['result'] = -1
                break
            self.update_stats()
            self.calc_thresh()
            callback(self, 1)
            self.prev_model = None
            if self.check_thresh():
                self.stats['result'] = 1
                break
        callback(self, 2)
        return self.stats['result']

    @abc.abstractmethod
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


# =============================================================================
# Fitter with ground truth
# =============================================================================


class GroundedFitter(SynapseFitter):
    """SynapseFitter where groud-truth is known.

    Parameters
    ----------
    data : SimPlasticitySequence
        The simulated data used to fit the synapse model.
    model : SynapseIdModel
        The initial guess/current estimate of the synapse model.
    truth : SynapseIdModel
        The model used to generate `data`.
    Other keywords added to `self.stats` (see `SynapseFitter`).

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
    """
    data: SimPlasticitySequence
    truth: SynapseIdModel

    def __init__(self, data: SimPlasticitySequence, model: SynapseIdModel,
                 truth: SynapseIdModel, **kwds) -> None:
        """SynapseFitter where groud-truth is known.

        Parameters
        ----------
        data : SimPlasticitySequence
            The data used to fit the synapse model.
        model : SynapseIdModel
            The initial guess/current estimate of the synapse model
        truth : SynapseIdModel
            The model used to generate `data`.

        Subclass should set `self.stats` items `['true_like', 'y_scale',
        'true_dlike']`.
        """
        super().__init__(data, model, **kwds)
        self.truth = truth
        self.stats.update(true_dmodel=(self.truth - self.model).norm(),
                          true_like=0., x_scale=self.truth.norm())

    def __format__(self, format_spec: str) -> str:
        """Printing stats"""
        templates = _frmt_vars(format_spec)
        disp = [templates[k] for k in ['t', 'tr', 'fit', 'dtr', 'trx']]
        if self.prev_model is not None and not format_spec:
            del disp[1]
            del disp[-2:]
        if self.opt.verbose % 2 and np.isfinite(self.stats['dnlog_like']):
            disp += [templates['df'], templates['dx']]
        elif self.opt.verbose % 2 and format_spec:
            disp += [templates['df'].replace('{dnlog_like:.2g}', '-'),
                     templates['dx'].replace('{dmodel:.2g}', '-')]
        return _frmt_out(disp, self.stats, format_spec)

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = super().__repr__()
        insert = "    truth = "
        insert += repr(self.data).replace("\n", "\n" + " " * 12) + ",\n"
        ind = rpr.find("    stats = ")
        rpr = rpr[:ind] + insert + rpr[ind:]
        return rpr

    @abc.abstractmethod
    def update_stats(self) -> None:
        """Calculate stats for termination.

        Subclass must update `self.stats` fields `[nlog_like, dnlog_like,
        true_dlike]`.
        """
        super().update_stats()
        self.stats['true_dmodel'] = (self.truth - self.model).norm()

    # @abc.abstractmethod
    # def update_model(self) -> None:
    #     """Perform a single update of the model"""
    #     super().update_model()

    # @abc.abstractmethod
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


# =============================================================================
Callback = Callable[[SynapseFitter, int], None]
