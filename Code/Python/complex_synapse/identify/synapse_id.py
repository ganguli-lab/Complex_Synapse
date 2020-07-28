"""Model that is being fit to data
"""
from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np

import numpy_linalg as la
from numpy_linalg import convert as cvl
from sl_py_tools.arg_tricks import default_non_eval as default
from sl_py_tools.containers import tuplify
from sl_py_tools.numpy_tricks.markov import isstochastic_c

from ..builders import RNG
from ..synapse_base import ArrayLike, SynapseBase
from .plast_seq import SimPlasticitySequence


class SynapseIdModel(SynapseBase):
    """Model that is being fit to data

    Parameters (and attributes)
    ---------------------------
    plast : array_like, (P,M,M), float[0:1]
        potentiation/depression transition probability matrix.
    frac : array_like, (P,), float[0:1]
        fraction of events that are potentiating/depressing.
    readout : array_like, (M,), int[0:R].
        id of readout when in that state.
    initial : array_like, (M,) float[0:1]
        distribution of initial state

    Properties
    ----------
    nstate : int
        number of states, M.
    nplast : int
        number of plasticity types, P.
    nreadout : int
        Number of readout values, R
    """
    # id of readout when in that state, (M,), int[0:R]
    readout: la.lnarray
    # distribution of initial state, (M,) float[0:1]
    initial: la.lnarray

    def __init__(self, plast: la.lnarray,
                 frac: ArrayLike = 0.5,
                 initial: Optional[ArrayLike] = None,
                 readout: Optional[ArrayLike] = None,
                 ) -> None:
        super().__init__(plast, frac)
        nst = self.nstate
        dtype = self.plast.dtype
        self.initial = default(initial, la.asarray, la.full(nst, 1/nst, dtype))
        self.readout = default(readout, la.asarray, la.arange(nst))
        if isstochastic_c(self.plast):
            self.plast += la.identity(nst, dtype)

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = super().__repr__()
        insert = f"    initial = {self.initial!r},\n"
        insert += f"    readout = {self.readout!r},\n"
        rpr = (rpr[:-1] + insert + rpr[-1])
        return rpr

    def __array_ufunc__(self, ufunc: np.ufunc, method: str, *inputs, **kwargs):
        """Handling ufuncs with SynapseIdModels
        """
        base_result = super().__array_ufunc__(ufunc, method, *inputs, **kwargs)

        args, _ = cvl.conv_in_attr('initial', SynapseIdModel, inputs)
        conv = [True] + [False] * (ufunc.nout-1)
        outs, conv = cvl.conv_in_attr('initial', SynapseIdModel, kwargs, conv)

        results = self.initial.__array_ufunc__(ufunc, method, *args, **kwargs)
        return cvl.conv_out_attr(base_result, 'initial', results, outs, conv)

    def moveaxis(self, source: Union[int, Sequence[int]],
                 destination: Union[int, Sequence[int]]) -> SynapseIdModel:
        """Change order of axes in self.plast and self.ini6tial.

        Parameters
        ----------
        source : Union[int, Tuple[int, ...]]
            Position of axis/axes to move
        destination : Union[int, Tuple[int, ...]]
            New position of axis/axes

        Returns
        -------
        newobj : SynapseIdModel
            Model with axes moved. Its attributes are views of the originals
        """
        newobj = self.view()
        newobj.plast = np.moveaxis(self.plast, source, destination)
        newobj.initial = np.moveaxis(self.initial, source, destination)

    @property
    def nreadout(self) -> int:
        """Number of readout values, R
        """
        return self.readout.max() + 1

    def readout_indicator(self) -> la.lnarray:
        """1 in sets with each readout value, 0 elsewhere (R,M)
        """
        indicators = la.zeros((self.nreadout, self.nstate))
        for i in range(self.nreadout):
            indicators[i, self.readout == i] = 1
        return indicators

    def readout_proj(self) -> la.lnarray:
        """Projection onto sets with each readout value, (R,M,M)
        """
        return self.readout_indicator().c * la.identity(self.nstate)

    def reorder(self, inds: ArrayLike) -> None:
        """Put the states in a new order, in-place.

        Parameters
        ----------
        inds : ArrayLike, (M,)
            `inds[i] = j` means state `j` moves to position `i`.
        """
        self.plast = self.plast[(...,) + np.ix_(inds, inds)]
        self.initial = self.initial[..., inds]
        self.readout = self.readout[..., inds]

    def sort(self, group: bool = True) -> None:
        """Sort the states by decreasing eta^w

        Parameters
        ----------
        group : bool, optional
            Sort by `-eta` within groups of common readout? By default `True`
        """
        markov = (self.frac.s * self.plast).sum(-3)
        eta = - (np.ones_like(markov) - markov).inv @ self.readout
        inds = np.lexsort((-eta, self.readout)) if group else np.argsort(-eta)
        self.reorder(inds)

    def simulate(self, ntime: Optional[int] = None,
                 nexpt: Union[None, int, Sequence[int]] = None,
                 plast_seq: Optional[np.ndarray] = None,
                 rng: np.random.Generator = RNG,
                 ) -> SimPlasticitySequence:
        """Simulate markov process

        Parameters
        ----------
        ntime : Optional[int], optional
            Numper of time-steps, T, by default None
        nexpt : Union[None, int, Sequence[int]], optional
            Number of experiments, E, by default None
        plast_seq : Optional[np.ndarray], optional
            Sequence of plasticity types, by default None

        Returns
        -------
        simulation : SimPlasticitySequence, (E,T)
            The result of the simulation
        """
        # only simulate a single model
        assert self.plast.ndim == 3
        if plast_seq is None:
            nexpt = default(nexpt, tuplify, ())
            ntime = default(ntime, int, 1)
            plast_seq = rng.choice(la.arange(self.nplast),
                                   size=nexpt + (ntime,), p=self.frac)
        else:
            plast_seq = la.asarray(plast_seq)
            nexpt = plast_seq.shape[:-1]
            ntime = plast_seq.shape[-1]

        state_ind = la.arange(self.nstate)
        jump = self.plast.flattish(0, -1)
        state_ch = la.array([rng.choice(state_ind, size=plast_seq.shape, p=p)
                             for p in jump]).foldaxis(0, self.plast.shape[:-1])
        states = np.empty_like(plast_seq)
        states[..., 0] = rng.choice(state_ind, size=nexpt, p=self.initial)
        for i in range(ntime - 1):
            states[..., i+1] = state_ch[plast_seq[..., i], states[..., i],
                                        ..., i]
        readouts = self.readout[states]
        return SimPlasticitySequence(plast_seq, readouts, states)

    def _elements(self) -> la.lnarray:
        """All elements of self.plast and elf.initial, concatenated (PM**2+M,)
        """
        return np.concatenate(np.broadcast_arrays(
            self.plast.flattish(-3), self.initial, subok=True), -1)

    def norm(self, **kwargs) -> la.lnarray:
        """Norm of vector of parameters.

        Parameters
        ----------
        All passed to `np.linalg.norm`.

        Returns
        -------
        norm : la.lnarray, ()
            Norm of parameters.
        """
        return np.linalg.norm(self._elements(), axis=-1, **kwargs)

    def kl_div(self, other: SynapseIdModel) -> la.lnarray:
        """Kullback-Leibler divergence between plast & initial of self & other
        """
        mine = self._elements()
        per_param = mine * np.log((self / other)._elements())
        mine, per_param = np.broadcast_arrays(mine, per_param, subok=True)
        per_param[np.isclose(0, mine)] = 0.
        return per_param.sum(-1)


# TODO: plot, likelihood
