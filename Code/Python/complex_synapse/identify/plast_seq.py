"""Storing the results of experiments
"""
from __future__ import annotations

from typing import List, Tuple, Union, Sequence

import numpy as np
import matplotlib as mpl

import numpy_linalg as la
import sl_py_tools.containers as _cn
import sl_py_tools.matplotlib_tricks as _mpt

from .. import synapse_base as _sb

# =============================================================================


class PlasticitySequence:
    """The results of an experiment

    Parameters (and attributes)
    ---------------------------
    plast_type : ArrayLike, (T-1,E), int[0:P]
        id of plasticity type after each time-step
    readouts : ArrayLike, (T,E), int[0:R]
        id of readout from synapse at each time-step
    t_axis : int
        Which axis is the time axis?

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
    # id of plasticity type after each time-step, (T-1,E), int[0:P]
    plast_type: la.lnarray
    # id of readout from synapse at each time-step, (T,E), int[0:R]
    readouts: la.lnarray
    # Which axis is the time-axis
    t_axis: int

    def __init__(self, plast_type: _sb.ArrayLike, readouts: _sb.ArrayLike,
                 t_axis: int = 0) -> None:
        """The results of an experiment

        Parameters
        ----------
        plast_type : np.ndarray, (T-1,E), int[0:P]
            id of plasticity type after each time-step
        readouts : np.ndarray, (T,E), int[0:R]
            id of readout from synapse at each time-step
        """
        self.plast_type = la.asanyarray(plast_type)
        self.readouts = la.asanyarray(readouts)
        self.t_axis = t_axis % self.readouts.ndim
        if self.plast_type.shape[self.t_axis] >= self.ntime:
            p_ind = np.s_[:,]*self.t_axis + np.s_[:self.ntime-1,]
            self.plast_type = self.plast_type[p_ind]

    def __getitem__(self, ind: Inds) -> PlasticitySequence:
        return self.view(plast_type=self.plast_type[ind],
                         readouts=self.readouts[ind])

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = type(self).__name__ + "(\n"
        with np.printoptions(threshold=20):
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
        attrs = _sb.array_attrs(self)
        attrs.update(kwds)
        return type(self)(**attrs)

    def copy(self, order: str = 'C', **kwargs) -> PlasticitySequence:
        """Copy of object, with copies of array attributes

        Requires `__init__` parameter names to be the same as attribute names.
        """
        attrs = _sb.array_attrs(self)
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
        ndim = self.readouts.ndim
        source = tuple(x % ndim for x in _cn.tuplify(source))
        if self.t_axis in source:
            destination = _cn.tuplify(destination)
            newobj.t_axis = destination[source.index(self.t_axis)] % ndim
        return newobj

    def move_t_axis(self, destination: int) -> PlasticitySequence:
        """Change position of time axis in self.readouts and self.plast_type.

        Parameters
        ----------
        destination : int
            New position of time axis

        Returns
        -------
        newobj : PlasticitySequence
            Model with axes moved. Its attributes are views of the originals
        """
        return self.moveaxis(self.t_axis, destination)

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
        return self.readouts.shape[self.t_axis]

    @property
    def nexpt(self) -> Tuple[int, ...]:
        """Number of experiment sequences stored, (E,)
        """
        shape = self.readouts.shape
        return shape[:self.t_axis] + shape[self.t_axis+1:]

    @property
    def shape(self) -> Tuple[int, ...]:
        """Number of time-steps, experiment sequences stored, (T,E)
        """
        return self.readouts.shape

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
        return _dwells(self.readouts.moveaxis(self.t_axis, -1), self.nreadout)

    def plot(self, hnds: Sequence[ImHandle], **kwds) -> List[Image]:
        """Plot heatmaps for plast_type and readouts

        Parameters
        ----------
        hnds : Sequence[Union[Image, Axes]], (2,)
            Axes to plot on, or Image objects to update with new data

        Returns
        -------
        imh: List[Image], (2,)
            Image objects for the plots
        """
        if self.nexpt:
            raise ValueError("Can only plot 1 scalar experiment at a time. " +
                             f"We have nexpt={self.nexpt}")
        kwds.setdefault('cmap', 'tab20b')
        nplast = kwds.pop('nplast', self.nplast)
        nreadout = kwds.pop('nreadout', self.nreadout)
        kwds['norm'] = _int_bdry_norm(nplast, kwds['cmap'])
        plast_type = _pad(nplast - 1 - self.plast_type.r, 1)
        pth = set_plot(hnds[0], plast_type, **kwds)
        kwds['norm'] = _int_bdry_norm(nreadout, kwds['cmap'])
        roh = set_plot(hnds[1], self.readouts.r, **kwds)
        return [pth, roh]


class SimPlasticitySequence(PlasticitySequence):
    """The results of a simulation

    Parameters (and attributes)
    ---------------------------
    plast_type : ArrayLike, (T-1,E), int[0:P]
        id of plasticity type after each time-step
    readouts : ArrayLike, (T,E), int[0:R]
        id of readout from synapse at each time-step
    states : ArrayLike, (T,E), int[0:M]
        which state it is in at each time-step
    t_axis : int
        Which axis is the time axis?

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

    def __init__(self, plast_type: _sb.ArrayLike, readouts: _sb.ArrayLike,
                 states: _sb.ArrayLike, t_axis: int = 0) -> None:
        """The results of an simulation

        Parameters
        ----------
        plast_type : ArrayLike, (T-1,E), int[0:P]
            id of plasticity type after each time-step
        readouts : ArrayLike, (T,E), int[0:R]
            id of readout from synapse at each time-step
        states : ArrayLike, (T,E), int[0:M]
            which state it is in at each time-step
        """
        super().__init__(plast_type, readouts, t_axis=t_axis)
        self.states = la.asanyarray(states)

    def __getitem__(self, ind: Inds) -> SimPlasticitySequence:
        newobj = super().__getitem__(ind)
        newobj.states = self.states[ind]
        return newobj

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = super().__repr__()
        with np.printoptions(threshold=90):
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
        return _dwells(self.states.moveaxis(self.t_axis, -1), self.nstate)

    def plot(self, hnds: Sequence[Handle], **kwds) -> List[Plot]:
        """Plot heatmaps for plast_type, readouts and Line for states

        Parameters
        ----------
        hnds : Sequence[Union[Image, Axes]], (3,)
            Axes to plot on, or Image/Lines to update with new data

        Returns
        -------
        imh: List[Image], (3,)
            Image/Line objects for the plots
        """
        nplast = kwds.pop('nplast', self.nplast)
        nreadout = kwds.pop('nreadout', self.nreadout)
        label = kwds.pop('label', 'True path')
        line_opts = kwds.pop('line_opts', {})
        imh = super().plot(hnds[:-1], nplast=nplast, nreadout=nreadout, **kwds)
        kwds.update(line_opts)
        kwds['line'] = True
        kwds['label'] = label
        kwds.pop('cmap', None)
        imh.append(set_plot(hnds[-1], self.states, **kwds))
        return imh


# =============================================================================
# Helpers
# =============================================================================


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
        `dwells[i]` is an array of dwell times in value `i`.
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


def _pad(array: np.ndarray, axis: int, end: bool = True) -> np.ma.masked_array:
    """Pad with invalid along one axis"""
    width = np.zeros((array.ndim, 2), int)
    width[axis, int(end)] = 1
    padded = np.pad(array, width, constant_values=-1)
    return np.ma.masked_equal(padded, -1, copy=False)


def _int_bdry_norm(nval: int, cmap: str) -> mpl.colors.BoundaryNorm:
    """BoundaryNorm map for integer values"""
    cmap = mpl.cm.get_cmap(cmap)
    return mpl.colors.BoundaryNorm(np.arange(nval + 1.) - 0.5, cmap.N)


def set_plot(handle: Handle, data: np.ndarray, **kwds) -> Plot:
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
    trn = kwds.pop('trn', False)
    line = kwds.pop('line', False)
    if line:
        data = _mpt.stepify_data(np.arange(-0.5, data.size), data)
        data = data[::-1] if trn else data
    else:
        data = data.T if trn else data
    if isinstance(handle, mpl.axes.Axes):
        if line:
            return handle.plot(*data, **kwds)[0]
        return handle.matshow(data, **kwds)
    handle.set_data(data)
    return handle


# =============================================================================
# Hinting Aliases
# =============================================================================
Inds = Tuple[Union[int, slice], ...]
Image = mpl.image.AxesImage
Handle = Union[mpl.axes.Axes, mpl.lines.Line2D, Image]
Plot = Union[mpl.lines.Line2D, Image]
ImHandle = Union[mpl.axes.Axes, Image]
