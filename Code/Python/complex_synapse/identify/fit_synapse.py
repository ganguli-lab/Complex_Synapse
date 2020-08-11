"""Classes for updating synapse models to fit data
"""
from __future__ import annotations

import abc
from numbers import Number
from typing import Callable, ClassVar, Dict, Optional, Tuple, Union

import numpy as np

from sl_py_tools.arg_tricks import default, default_eval

from .synapse_id import SynapseIdModel, valid_values
from .plast_seq import PlasticitySequence, Axes, Image, Line
# =============================================================================
MESSAGES = {
    -1: "Error, invalid model generated.",
    0: "Failure, maximum iterations reached.",
    1: "Success, changes in log-likelihood and model estimate below threshold."
}


def print_callback(obj: SynapseFitter, pos: int) -> None:
    """Callback that prints fitter state as appropriate"""
    if pos == 2 and obj.opt['verbose'] >= 0:
        print(obj, MESSAGES[obj.stats['result']], sep='\n')
    elif pos == 1 and obj.opt['verbose'] >= 2:
        print(obj)


# =============================================================================


class SynapseFitter(abc.ABC):
    """Base class for synapse fitters.

    Parameters
    ----------
    data: PlasticitySequence
        The data used to fit the synapse model.
    model: SynapseIdModel
        The initial guess/current estimate of the synapse model.
    Other keywords added to `self.stats` below.

    Other Attributes
    ----------------
    prev_model: Optional[SynapseIdModel]
        During itteration or after error: Previous model estimate.
        Otherwise: `None`.
    stats: Dict[str, Number]
        Statistics of current state. See below.
    opt: ClassVar[Dict[str, Number]]
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
    atolx: ClassVar[float] = 1e-5
        Absolute tolerance for `dmodel`.
    atoly: ClassVar[float] = 1e-5
        Absolute tolerance for `dnlog_like`.
    rtol: ClassVar[float] = 1e-5
        Relative tolerance. Multiplies `x_scale` and `y_scale` if given.
    max_it: ClassVar[int] = 1e3
        Maximum number of iterations
    verbose: ClassVar[int] = -2
        When statistics are printed, and how verbose:
            -2: print called
            -1: print called, detailed
            0: end of iteration
            1: end of iteration, detailed
            2: each iteration
            3: each iteration, detailed
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
                                        'verbose': -2}

    def __init__(self, data: PlasticitySequence, model: SynapseIdModel, **kwds
                 ) -> None:
        self.data = data
        self.model = model
        self.prev_model = None
        self.stats = {'nit': 0, 'nlog_like': np.inf, 'dnlog_like': np.inf,
                      'dmodel': np.inf, 'x_thresh': 0., 'y_thresh': 0.,
                      'x_scale': None, 'y_scale': None}
        self.stats.update(kwds)

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
        """Calculate stats for termination.

        Subclass must update `stats[nlog_like]` and `stats[dnlog_like]`.
        """
        self.stats['dmodel'] = (self.model - self.prev_model).norm()

    @abc.abstractmethod
    def update_model(self) -> None:
        """Perform a single update of the model"""

    def valid(self) -> bool:
        """Check that current model estimate is valid"""
        return not self.model.nmodel and valid_values(self.model)

    def __str__(self) -> str:
        """Printing stats"""
        disp = f"nit = {self.stats['nit']}, "
        disp += f"-log P(data|fit) = {self.stats['nlog_like']:.3f},"
        if self.opt['verbose'] % 2:
            disp += f"\nlog[P(data|fit) / P(data|prev)] = "
            disp += f"{self.stats['dnlog_like']:.3f}, "
            disp += f"||fit-prev|| = {self.stats['dmodel']:.3f},"
        return disp

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
        self.stats['result'] = 0
        callback(self, 0)
        for i in range(self.opt['max_it']):
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
Callback = Callable[[SynapseFitter, int], None]
