"""Model that is being fit to data
"""
from __future__ import annotations
from complex_synapse.synapse_mem import normalise

from typing import List, Optional, Sequence, Tuple, Union, Callable, Dict

import numpy as np
import matplotlib as mpl

import numpy_linalg as la
from numpy_linalg import convert as cvl
from sl_py_tools.arg_tricks import default_non_eval as default
from sl_py_tools.containers import tuplify
import sl_py_tools.numpy_tricks.markov as ma

from ..builders import RNG
from ..synapse_base import ArrayLike, SynapseBase
from .plast_seq import SimPlasticitySequence, Axes, Image, set_plot


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
        if ma.isstochastic_c(self.plast):
            self.plast += la.identity(nst, dtype) * max(- self.plast.min(), 1.)

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

    @property
    def nreadout(self) -> int:
        """Number of readout values, R
        """
        return self.readout.max() + 1

    def moveaxis(self, source: Union[int, Sequence[int]],
                 destination: Union[int, Sequence[int]]) -> SynapseIdModel:
        """Change order of axes in self.plast and self.initial.

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

    def normalise(self) -> None:
        """normalise plast and initial, in place"""
        ma.stochastify_d(self.plast)
        ma.stochastify_d(self.initial)

    def set_init(self, init: Optional[ArrayLike] = None) -> None:
        """Set initial to steady-state, in place"""
        if init is None:
            markov = (self.frac.s * self.plast).sum(-3)
            self.initial = ma.calc_peq_d(markov)
        else:
            self.initial = la.asarray(init)

    def reorder(self, inds: ArrayLike) -> None:
        """Put the states into a new order, in-place.

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
        fundi = np.ones_like(markov) + la.identity(self.nstate) - markov
        eta = - fundi.inv @ self.readout
        inds = np.lexsort((-eta, self.readout)) if group else np.argsort(-eta)
        self.reorder(inds)

    def readout_indicator(self) -> la.lnarray:
        """1 in sets with each readout value, 0 elsewhere. (R,M)
        """
        indicators = la.zeros((self.nreadout, self.nstate))
        for i in range(self.nreadout):
            indicators[i, self.readout == i] = 1
        return indicators

    def updaters(self) -> Tuple[la.lnarray, la.lnarray]:
        """Projected plasticity matrices, (R,P,M,M), (R,M)
        """
        indic = self.readout_indicator()
        updaters = self.plast.expand_dims(-4) * indic.expand_dims((-3, -2))
        initial = self.initial * indic
        return updaters, initial

    def _elements(self) -> la.lnarray:
        """All elements of self.plast and self.initial, concatenated (PM**2+M,)
        """
        return np.concatenate((self.plast.flattish(-3), self.initial), -1)

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
        plast_seq : Optional[np.ndarray], optional, (T,E)
            Sequence of plasticity types, by default None

        Returns
        -------
        simulation : SimPlasticitySequence, (T,E)
            The result of the simulation
        """
        # only simulate a single model
        assert self.plast.ndim == 3
        if plast_seq is None:
            nexpt = default(nexpt, tuplify, ())
            ntime = default(ntime, int, 1)
            # (T,E)
            plast_seq = rng.choice(la.arange(self.nplast), p=self.frac,
                                   size=(ntime,) + nexpt)
        else:
            plast_seq = la.asarray(plast_seq)
            nexpt = plast_seq.shape[1:]
            ntime = plast_seq.shape[0]

        state_ind = la.arange(self.nstate)
        expt_ind = tuple(la.arange(expt) for expt in nexpt)
        # (PM,M)
        jump = self.plast.flattish(0, -1)
        # (P,M,T,E)
        state_ch = la.array([rng.choice(state_ind, size=plast_seq.shape, p=p)
                             for p in jump]).foldaxis(0, self.plast.shape[:-1])
        # (T,E)
        states = np.empty_like(plast_seq)
        # (E,)
        states[0] = rng.choice(state_ind, size=nexpt, p=self.initial)
        for i in range(ntime - 1):
            # (E,)
            states[i+1] = state_ch[(plast_seq[i], states[i], i) + expt_ind]
        # (T,E)
        readouts = self.readout[states]
        return SimPlasticitySequence(plast_seq, readouts, states, t_axis=0)

    def plot(self, axs: Sequence[Union[Axes, Image]], **kwds) -> List[Image]:
        """Plot heatmaps for initial and plast

        Parameters
        ----------
        axs : Sequence[Union[Image, Axes]], (P+1,)
            Axes to plot on, or Images to update with new data

        Returns
        -------
        imh: List[Image], (P+1,)
            Image objects for the plots
        """
        kwds['norm'] = mpl.colors.Normalize(0., 1.)
        imh = [set_plot(axs[0], self.initial, **kwds)]
        for axh, mat in zip(axs[1:], self.plast):
            imh.append(set_plot(axh, mat, **kwds))
        return imh

    @classmethod
    def build(cls, builder: Callable[..., Dict[str, la.lnarray]],
              nst: int, frac: ArrayLike = 0.5,
              extra_args=(), **kwargs) -> SynapseBase:
        """Build model from function.

        Parameters
        ----------
        builder : function
            function object with parameters (n,...) that returns dictionary
            {plast, weight, etc}, ignoring any extra elements.
        nst: int
            total number of states, passed to `builder`.
        frac : float
            fraction of events that are potentiating, default=0.5.
        extra_args, **kwargs
            *extra_args, **kwargs: passed to `builder`.

        Returns
        -------
        synobj
            SynapseBase instance
        """
        keys = builder(nst, *extra_args, **kwargs)
        keys['frac'] = frac
        keys.pop('signal')
        keys.setdefault('readout', _weight_readout(keys.pop('weight')))
        obj = cls(**keys)
        obj.normalise()
        obj.set_init()
        return obj


# =============================================================================


def _weight_readout(weight: ArrayLike) -> la.lnarray:
    """Convert weight vector to readout vector"""
    return None if weight is None else np.unique(weight, return_inverse=True
                                                 )[1].astype(int)
