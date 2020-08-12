"""Classes for updating synapse models to fit data
"""
from __future__ import annotations

import abc
from numbers import Number
from typing import Callable, ClassVar, Dict, Optional, Tuple, Union

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


def print_callback(obj: SynapseFitter, pos: int) -> None:
    """Callback that prints fitter state as appropriate"""
    if pos == 0 and obj.opt['verbose'] >= 0:
        gap = ('',) if obj.opt['verbose'] >= 2 else ()
        print('Before:', obj, *gap, sep='\n')
    elif pos == 1 and obj.opt['verbose'] >= 2:
        if obj.stats['nit'] % obj.opt['disp_step'] == 0:
            print('*', obj)
    elif pos == 2 and obj.opt['verbose'] >= 0:
        print('', 'After:', obj, MESSAGES[obj.stats['result']], sep='\n')


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
    Other keywords added to `self.stats` below.

    Other Attributes
    ----------------
    prev_model : Optional[SynapseIdModel]
        During itteration or after error: Previous model estimate.
        Otherwise: `None`.
    stats : Dict[str, Number]
        Statistics of current state. See below.
    opt : ClassVar[Dict[str, Number]]
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
    atolx : ClassVar[float] = 1e-5
        Absolute tolerance for `dmodel`.
    atoly : ClassVar[float] = 1e-5
        Absolute tolerance for `dnlog_like`.
    rtol : ClassVar[float] = 1e-5
        Relative tolerance. Multiplies `x_scale` and `y_scale` if given.
    max_it : ClassVar[int] = 1e3
        Maximum number of iterations
    verbose : ClassVar[int] = -2
        When statistics are printed, and how verbose:
            -2: print called
            -1: print called, detailed
            0: end of iteration
            1: end of iteration, detailed
            2: each iteration
            3: each iteration, detailed
    disp_step : ClassVar[int] = -2
        Display progress update every `disp_step` iterations.
    All of the above are stored in `SynapseFitter.opt`.
    """
    data: PlasticitySequence
    model: SynapseIdModel
    prev_model: Optional[SynapseIdModel]
    # Stats of current state:
    stats: Dict[str, Number]
    # options:
    opt: ClassVar[Dict[str, Number]] = {'atolx': 1e-5, 'atoly': 1e-5,
                                        'rtol': 1e-5, 'max_it': 1000,
                                        'verbose': -2, 'disp_step': 50}

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
        self.stats.update(kwds)

    def __str__(self) -> str:
        """Printing stats"""
        disp = f"nit = {self.stats['nit']}, "
        disp += f"-log P(data|fit) = {self.stats['nlog_like']:.3f}, "
        if self.opt['verbose'] % 2 and np.isfinite(self.stats['dnlog_like']):
            disp += f"\nlog[P(data|fit) / P(data|prev)] = "
            disp += f"{self.stats['dnlog_like']:.3g}, "
            disp += f"||fit-prev|| = {self.stats['dmodel']:.3g}, "
        return disp

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
        self.stats['x_thresh'] = self.opt['atolx'] + self.opt['rtol'] * x_scale
        self.stats['y_thresh'] = self.opt['atoly'] + self.opt['rtol'] * y_scale

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
        count = undcount if self.opt['verbose'] >= 2 else dcount
        callback(self, 0)
        self.stats['result'] = 0
        for i in count('iteration', self.opt['max_it'],
                       disp_step=self.opt['disp_step']):
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

    def plot_states(self, handle: Union[Axes, Image, Line],
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
        return NotImplemented

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
    true_nlog_like : float
        Negative log-likelihood of `data` given `truth`.
    true_dmodel : float
        Distance between `truth` and `model`.
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

        Subclass should set `self.stats` items `['true_nlog_like', 'y_scale']`.
        """
        super().__init__(data, model, **kwds)
        self.truth = truth
        self.stats.update(true_dmodel=(self.truth - self.model).norm(),
                          true_nlog_like=0., x_scale=self.truth.norm())

    def __str__(self) -> str:
        """Printing stats"""
        disp = super().__str__()
        if self.prev_model is None:
            add = f"-log P(data|truth) = {self.stats['true_nlog_like']:.3f}, "
            lin = disp.find(',') + 2
            disp = disp[:lin] + add + disp[lin:]
        if self.opt['verbose'] % 2:
            dtrue_like = self.stats['nlog_like'] - self.stats['true_nlog_like']
            disp += f"\nlog[P(data|true) / P(data|fit)] = {dtrue_like:.3g}, "
            disp += f"||true-fit|| = {self.stats['true_dmodel']:.3g},"
        return disp

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

        Subclass must update `stats[nlog_like]` and `stats[dnlog_like]`.
        """
        super().update_stats()
        self.stats['true_dmodel'] = (self.truth - self.model).norm()

    @abc.abstractmethod
    def update_model(self) -> None:
        """Perform a single update of the model"""
        super().update_model()

# =============================================================================
Callback = Callable[[SynapseFitter, int], None]
