"""Storing the results of experiments
"""
from __future__ import annotations

from typing import List, Tuple, Union, Sequence

import numpy as np
import matplotlib as mpl
from matplotlib.axes import Axes
from matplotlib.lines import Line2D as Line
from matplotlib.image import AxesImage as Image

import numpy_linalg as la
from sl_py_tools.containers import tuplify

from ..synapse_base import array_attrs, ArrayLike


class PlasticitySequence:
    """The results of an experiment

    Parameters (and attributes)
    ---------------------------
    plast_type : ArrayLike, (T,E), int[0:P]
        id of plasticity type after each time-step
    readouts : ArrayLike, (T,E), int[0:R]
        id of readout from synapse at each time-step

    Properties
    ----------
    nstate : int
        number of states, M.
    nplast : int
        number of plasticity types, P.
    nreadout : int
        Number of readout values, R
    ntime : int
        Number of time-steps, T
    nexpt : int
        Number of experiment sequences stored, E
    """
    # id of plasticity type after each time-step, (T,E), int[0:P]
    plast_type: la.lnarray
    # id of readout from synapse at each time-step, (T,E), int[0:R]
    readouts: la.lnarray
    # Which axis is the time-axis
    t_axis: int

    def __init__(self, plast_type: ArrayLike, readouts: ArrayLike,
                 t_axis: int = 0) -> None:
        """The results of an experiment

        Parameters
        ----------
        plast_type : np.ndarray, (T,E), int[0:P]
            id of plasticity type after each time-step
        readouts : np.ndarray, (T,E), int[0:R]
            id of readout from synapse at each time-step
        """
        self.plast_type = la.asanyarray(plast_type)
        self.readouts = la.asanyarray(readouts)
        self.t_axis = t_axis

    def __getitem__(self, ind: Tuple[Union[int, slice], ...]
                    ) -> PlasticitySequence:
        newobj = self.view()
        newobj.plast_type = self.plast_type[ind]
        newobj.readouts = self.readouts[ind]
        return newobj

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = type(self).__name__ + "(\n"
        rpr += "    plast_type = "
        rpr += repr(self.plast_type).replace("\n", "\n" + " " * 17) + ",\n"
        rpr += "    readouts = "
        rpr += repr(self.readouts).replace("\n", "\n" + " " * 15) + ",\n"
        rpr += ")"
        return rpr

    def __str__(self) -> str:
        """Short representation of object"""
        return f"{type(self).__name__} with T={self.ntime}, E={self.nexpt}"

    def view(self, **kwds) -> PlasticitySequence:
        """Copy of object, with views of array attributes

        Requires `__init__` parameter names to be the same as attribute names.
        """
        attrs = array_attrs(self)
        attrs.update(kwds)
        return type(self)(**attrs)

    def copy(self, order: str = 'C', **kwargs) -> PlasticitySequence:
        """Copy of object, with copies of array attributes

        Requires `__init__` parameter names to be the same as attribute names.
        """
        attrs = array_attrs(self)
        for k in attrs:
            if k not in kwargs.keys():
                kwargs[k] = attrs[k].copy(order)
        return type(self)(**kwargs)

    def moveaxis(self, source: Union[int, Sequence[int]],
                 destination: Union[int, Sequence[int]]) -> PlasticitySequence:
        """Change order of axes in self.readouts and self.plast_type.

        Parameters
        ----------
        source : Union[int, Tuple[int, ...]]
            Position of axis/axes to move
        destination : Union[int, Tuple[int, ...]]
            New position of axis/axes

        Returns
        -------
        newobj : PlasticitySequence
            Model with axes moved. Its attributes are views of the originals
        """
        newobj = self.view()
        newobj.readouts = np.moveaxis(self.readouts, source, destination)
        newobj.plast_type = np.moveaxis(self.plast_type, source, destination)
        source, destination = tuplify(source), tuplify(destination)
        if self.t_axis in source:
            newobj.t_axis = destination[source.index(self.t_axis)]
        return newobj

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
        return self.plast_type.shape[0]

    @property
    def nexpt(self) -> Tuple[int, ...]:
        """Number of experiment sequences stored, (E,)
        """
        return self.plast_type.shape[1:]

    @property
    def shape(self) -> Tuple[int, ...]:
        """Number of time-steps, experiment sequences stored, (T,E)
        """
        return self.plast_type.shape

    def frac(self) -> la.lnarray:
        """fraction of each plasticity type, (P,), float[0:1]
        """
        return np.histogram(self.plast_type, self.nplast, (0, self.nplast),
                            density=True)[0].view(la.lnarray)

    def readout_dwells(self) -> List[la.lnarray]:
        """Sets of dwell times for each readout value

        Returns
        -------
        dwells : List[la.lnarray] (R,)(?,)
            `dwells[i]` is an array of dwell times in readout group `i`.
        """
        return _dwells(self.readouts, self.nreadout)

    def plot(self, axs: Sequence[Union[Image, Axes]], **kwds) -> List[Image]:
        """Plot heatmaps for plast_type and readouts

        Parameters
        ----------
        axs : Sequence[Union[Image, Axes]], (2,)
            Axes to plot on, or Images to update with new data

        Returns
        -------
        imh: List[Image], (2,)
            Image objects for the plots
        """
        kwds['norm'] = mpl.colors.Normalize(0., 1. * self.nplast)
        imh = [set_plot(axs[0], self.plast_type, **kwds)]
        kwds['norm'] = mpl.colors.Normalize(0., 1. * self.nreadout)
        imh.append(set_plot(axs[1], self.readouts, **kwds))
        return imh


class SimPlasticitySequence(PlasticitySequence):
    """The results of a simulation

    Parameters (and attributes)
    ---------------------------
    plast_type : ArrayLike, (T,E), int[0:P]
        id of plasticity type after each time-step
    readouts : ArrayLike, (T,E), int[0:R]
        id of readout from synapse at each time-step
    states : ArrayLike, (T,E), int[0:M]
        which state it is in at each time-step

    Properties
    ----------
    nstate : int
        number of states, M.
    nplast : int
        number of plasticity types, P.
    nreadout : int
        Number of readout values, R
    ntime : int
        Number of time-steps, T
    nexpt : int
        Number of simulated sequences stored, E
    """
    # which state it is in at each time-step, (T,E), int[0:M]
    states: la.lnarray

    def __init__(self, plast_type: ArrayLike, readouts: ArrayLike,
                 states: ArrayLike, t_axis: int = 0) -> None:
        """The results of an simulation

        Parameters
        ----------
        plast_type : ArrayLike, (T,E), int[0:P]
            id of plasticity type after each time-step
        readouts : ArrayLike, (T,E), int[0:R]
            id of readout from synapse at each time-step
        states : ArrayLike, (T,E), int[0:M]
            which state it is in at each time-step
        """
        super().__init__(plast_type, readouts, t_axis=t_axis)
        self.states = la.asanyarray(states)

    def __getitem__(self, ind: Tuple[Union[int, slice], ...]
                    ) -> SimPlasticitySequence:
        newobj = super().__getitem__(ind)
        newobj.states = self.states[ind]
        return newobj

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = super().__repr__()
        insert = "    states = "
        insert += repr(self.states).replace("\n", "\n" + " " * 13) + ",\n"
        rpr = rpr[:-1] + insert + rpr[-1]
        return rpr

    def moveaxis(self, source: Union[int, Sequence[int]],
                 destination: Union[int, Sequence[int]]
                 ) -> SimPlasticitySequence:
        """Change order of axes in self.readouts and self.plast_type.

        Parameters
        ----------
        source : Union[int, Tuple[int, ...]]
            Position of axis/axes to move
        destination : Union[int, Tuple[int, ...]]
            New position of axis/axes

        Returns
        -------
        newobj : SimPlasticitySequence
            Model with axes moved. Its attributes are views of the originals
        """
        newobj = super().moveaxis(source, destination)
        newobj.states = np.moveaxis(self.states, source, destination)
        return newobj

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

    def plot(self, axs: Sequence[Union[Image, Line, Axes]], **kwds
             ) -> List[Union[Image, Line]]:
        """Plot heatmaps for plast_type, readouts and Line for states

        Parameters
        ----------
        axs : Sequence[Union[Image, Axes]], (3,)
            Axes to plot on, or Image/Lines to update with new data

        Returns
        -------
        imh: List[Image], (3,)
            Image/Line objects for the plots
        """
        imh = super().plot(axs, **kwds)
        kwds['line'] = True
        imh.append(set_plot(axs[2], self.states, **kwds))
        return imh


def _dwells(values: np.ndarray, num: int) -> List[la.lnarray]:
    """Set of dwell times

    Parameters
    ----------
    values : np.ndarray (T,E)
        value at each time-step
    num : int
        number of values

    Returns
    -------
    dwells : List[la.lnarray] (num,)(?,)
        `dwells[i` is an array of dwell times in value `i`.
    """
    changes = np.diff(values.T).nonzero()
    before = values[changes]
    stays = np.diff(changes[-1], prepend=0).view(la.lnarray)
    if values.ndim > 1:
        new_expt = np.diff(changes[:-1]).any(axis=0).nonzero()
        stays[new_expt] = changes[-1][new_expt]
    dwells = []
    for i in range(num):
        dwells.append(stays[before == i])
    return dwells


def set_plot(handle: Union[Axes, Image, Line], data: np.ndarray, **kwds
             ) -> Union[Image, Line]:
    """Make/update heatmaps or Lines

    Parameters
    ----------
    handle : Union[Image, Line, Axes]
        Axes to plot on, or Image/Lines to update with new data
    data : np.ndarray
        Data to use for plot

    Returns
    -------
    imh: Union[Image, Line]
        Image/Line objects for the plots
    """
    line = kwds.pop('line', False)
    if isinstance(handle, Axes):
        if line:
            kwds.setdefault('where', 'mid')
            return handle.step(np.arange(data.size), data, **kwds)[0]
        return handle.matshow(data, **kwds)
    if isinstance(handle, Line):
        handle.set_ydata(data)
    else:
        handle.set_data(data)
    return handle
