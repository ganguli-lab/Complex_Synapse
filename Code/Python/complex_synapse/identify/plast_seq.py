"""Storing the results of experiments
"""
from __future__ import annotations
from typing import List, Tuple, Union
import numpy as np
import numpy_linalg as la


class PlasticitySequence:
    """The results of an experiment
    """
    # id of plasticity type after each time-step, (E,P), int[0:P]
    plast_type: la.lnarray
    # id of readout from synapse at each time-step, (E,P), int[0:R]
    readouts: la.lnarray

    def __init__(self, plast_type: np.ndarray, readouts: np.ndarray) -> None:
        self.plast_type = la.asanyarray(plast_type)
        self.readouts = la.asanyarray(readouts)

    def __getitem__(self, ind: Tuple[Union[int, slice], ...]
                    ) -> PlasticitySequence:
        return type(self)(self.plast_type[ind], self.readouts[ind])

    @property
    def nplast(self) -> int:
        """Minimum number of types of plasticity, P
        """
        return self.plast_type.max() + 1

    @property
    def nreadout(self) -> int:
        """Minimum number of readout values, R
        """
        return self.readouts.max() + 1

    @property
    def ntime(self) -> int:
        """Number of time-steps, T
        """
        return self.plast_type.shape[-1]

    @property
    def nexpt(self) -> int:
        """Number of experiment sequences stored, E
        """
        return 1 if self.readouts.ndim < 2 else self.plast_type.shape[-2]

    def readout_dwells(self) -> List[la.lnarray]:
        """Sets of dwell times for each readout value

        Returns
        -------
        dwells : List[la.lnarray] (R,)(?,)
            `dwells[i]` is an array of dwell times in readout group `i`.
        """
        return _dwells(self.readouts, self.nreadout)


class SimPlasticitySequence(PlasticitySequence):
    """The results of a simulation
    """
    # which state it is in at each time-step, (E,T), int[0:M]
    states: la.lnarray

    def __init__(self, states: np.ndarray, plast_type: np.ndarray,
                 readouts: np.ndarray) -> None:
        super().__init__(plast_type=plast_type, readouts=readouts)
        self.states = la.asanyarray(states)

    def __getitem__(self, ind: Tuple[Union[int, slice], ...]
                    ) -> SimPlasticitySequence:
        return type(self)(self.states[ind], self.plast_type[ind],
                          self.readouts[ind])

    @property
    def nstate(self) -> int:
        """Minimum number of states, M
        """
        return self.states.max() + 1

    def state_dwells(self) -> List[la.lnarray]:
        """Sets of dwell times for each state

        Returns
        -------
        dwells : List[la.lnarray] (M,)(?,)
            `dwells[i]` is an array of dwell times in state `i`.
        """
        return _dwells(self.states, self.nstate)


def _dwells(values: np.ndarray, num: int) -> List[la.lnarray]:
    """Set of dwell times

    Parameters
    ----------
    values : np.ndarray (E,T)
        value at each time-step
    num : int
        number of values

    Returns
    -------
    dwells : List[la.lnarray] (num,)(?,)
        `dwells[i` is an array of dwell times in value `i`.
    """
    changes = np.diff(values).nonzero()
    before = values[changes]
    stays = np.diff(changes[-1], prepend=0).view(la.lnarray)
    if values.ndim > 1:
        new_expt = np.diff(changes[:-1]).any(axis=0).nonzero()
        stays[new_expt] = changes[-1][new_expt]
    dwells = []
    for i in range(num):
        dwells.append(stays[before == i])
    return dwells

# TODO: plot
