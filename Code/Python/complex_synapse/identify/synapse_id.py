"""Model that is being fit to data
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib as mpl
import numpy as np

import numpy_linalg as la
import numpy_linalg.convert as cvl
import sl_py_tools.arg_tricks as ag
import sl_py_tools.containers as cn
import sl_py_tools.numpy_tricks.markov as ma

from .. import builders as _bld
from .. import synapse_base as _sb
from . import plast_seq as _ps

# =============================================================================


class SynapseIdModel(_sb.SynapseBase):
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
    nmodel : Tuple[int]
        Number and shape of models being broadcast.
    """
    # id of readout when in that state, (M,), int[0:R]
    readout: la.lnarray
    # distribution of initial state, (M,) float[0:1]
    initial: la.lnarray

    def __init__(self, plast: la.lnarray,
                 frac: _sb.ArrayLike = 0.5,
                 initial: Optional[_sb.ArrayLike] = None,
                 readout: Optional[_sb.ArrayLike] = None,
                 ) -> None:
        super().__init__(plast, frac)
        nst = self.nstate
        self.initial = ag.default_non_eval(initial, la.asarray, la.ones(nst)/nst)
        self.readout = ag.default_non_eval(readout, la.asarray, la.arange(nst))
        if ma.isstochastic_c(self.plast):
            self.plast += la.identity(nst) * max(- self.plast.min(), 1.)

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = super().__repr__()
        insert = f"    initial = {self.initial!r},\n"
        insert += f"    readout = {self.readout!r},\n"
        rpr = rpr[:-1] + insert + rpr[-1]
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

    @property
    def nmodel(self) -> Tuple[int, ...]:
        """Number and shape of models being broadcast."""
        return la.gufuncs.array_return_shape('(p,m,m),(m)->()',
                                             self.plast, self.initial)

    def moveaxis(self, source: Union[int, Sequence[int]],
                 destination: Union[int, Sequence[int]]) -> SynapseIdModel:
        """Change order of axes in `self.plast` and `self.initial`.

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
        """normalise `plast` and `initial`, in place"""
        ma.stochastify_d(self.plast)
        ma.stochastify_d(self.initial)

    def set_init(self, init: Optional[_sb.ArrayLike] = None) -> None:
        """Set initial to steady-state, in place

        Parameters
        ----------
        init : None or ArrayLike (M,), optional
            Vector to use as initial probability distribution.
            By default `None` -> use steady-state distribution.
        """
        if init is None:
            markov = (self.frac.s * self.plast).sum(-3)
            self.initial = ma.calc_peq_d(markov)
        else:
            self.initial = la.asarray(init)

    def reorder(self, inds: _sb.ArrayLike) -> None:
        """Put the states into a new order, in-place.

        Parameters
        ----------
        inds : ArrayLike[int] (M,)
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

    def readout_prob(self) -> la.lnarray:
        """Readout probabilities (R,M)

        rp[r,i] = Prob(readout=r|state=i).
        """
        indicators = la.zeros((self.nreadout, self.nstate))
        for i in range(self.nreadout):
            indicators[i, self.readout == i] = 1
        return indicators

    def updaters(self) -> Tuple[la.lnarray, la.lnarray]:
        """Projected plasticity matrices & initial, (R,P,M,M), (R,M)

        update[r,p,i,j] = Prob(i(t+1)=j, w(t+1)=r|i(t)=i, mu(t)=p),
        init[r,i] = Prob(w(0)=r, i(0)=i)
        """
        indic = self.readout_prob()
        updaters = self.plast.expand_dims(-4) * indic.expand_dims((-3, -2))
        initial = self.initial * indic
        return updaters, initial

    def elements(self) -> la.lnarray:
        """All elements of self.plast and self.initial, concatenated.

        Returns
        -------
        elems : la.lnarray (PM**2+M,)
            When `self.nmodel = ()`, concatenation of ravelled `self.plast` and
            `self.initial`. Otherwise, broadcasts over `nmodel` axes.
        """
        if not self.nmodel:
            return np.concatenate((self.plast.ravel(), self.initial))
        vectors = (self.plast.ravelaxes(-3), self.initial)
        bcast_vectors = la.gufuncs.broadcast_matrices('(a),(b)', *vectors)
        return np.concatenate(bcast_vectors, -1)

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
        return np.linalg.norm(self.elements(), axis=-1, **kwargs)

    def kl_div(self, other: SynapseIdModel) -> la.lnarray:
        """Kullback-Leibler divergence between plast & initial of self & other
        """
        mine, theirs = self.elements(), other.elements()
        per_param = mine * np.log(mine / theirs)
        mine, per_param = np.broadcast_arrays(mine, per_param, subok=True)
        per_param[np.isclose(0, mine)] = 0.
        return per_param.sum(-1)

    def simulate(self, ntime: Optional[int] = None,
                 nexpt: Union[None, int, Sequence[int]] = None,
                 plast_seq: Optional[np.ndarray] = None,
                 rng: np.random.Generator = _bld.RNG,
                 ) -> _ps.SimPlasticitySequence:
        """Simulate markov process

        If `plast_seq` is omitted, both `ntime` and `nexpt` are required.

        Parameters
        ----------
        ntime : Optional[int], optional
            Numper of time-steps, T, by default None
        nexpt : Union[None, int, Sequence[int]], optional
            Number of experiments, E, by default None
        plast_seq : Optional[np.ndarray], optional, (T,E)
            Sequence of plasticity types, by default None
        rng : np.random.Generator, optional
            random number generator, by default: builders.RNG

        Returns
        -------
        simulation : SimPlasticitySequence, (T,E)
            The result of the simulation
        """
        # only simulate a single model
        assert not self.nmodel
        if plast_seq is None:
            nexpt = ag.default_non_eval(nexpt, cn.tuplify, ())
            ntime = ag.default_non_eval(ntime, int, 1)
            # (T-1,E)
            plast_seq = rng.choice(la.arange(self.nplast), p=self.frac,
                                   size=(ntime - 1,) + nexpt)
        else:
            # (T-1,E)
            plast_seq = la.asarray(plast_seq)
            nexpt = plast_seq.shape[1:]
            ntime = plast_seq.shape[0] + 1

        state_ind = la.arange(self.nstate)
        expt_ind = np.ix_(*(la.arange(expt) for expt in nexpt))
        expt_ind = (expt_ind,) if len(nexpt) == 1 else expt_ind
        # (PM,M)
        jump = self.plast.ravelaxes(0, -1)
        # (P,M,T-1,E)
        state_ch = la.array([rng.choice(state_ind, size=plast_seq.shape, p=p)
                             for p in jump]).unravelaxis(0,
                                                         self.plast.shape[:-1])
        # (T,E)
        states = la.empty((ntime,) + nexpt, int)
        # (E,)
        states[0] = rng.choice(state_ind, size=nexpt, p=self.initial)
        for i, plt in enumerate(plast_seq):
            # (E,)
            states[i+1] = state_ch[(plt, states[i], i) + expt_ind]
        # (T,E)
        readouts = self.readout[states]
        return _ps.SimPlasticitySequence(plast_seq, readouts, states, t_axis=0)

    def plot(self, axs: Sequence[_ps.ImHandle], **kwds) -> List[_ps.Image]:
        """Plot heatmaps for initial and plast

        Parameters
        ----------
        axs : Sequence[Image or Axes], (P+1,)
            Axes to plot on, or Images to update with new data

        Returns
        -------
        imh: List[Image], (P+1,)
            Image objects for the plots
        """
        if self.nmodel:
            raise ValueError("Can only plot 1 scalar model at a time. " +
                             f"We have nmodel={self.nmodel}")
        trn = kwds.pop('trn', False)
        kwds.setdefault('norm', mpl.colors.Normalize(0., 1.))
        initial = self.initial.r if trn else self.initial.c
        imh = [_ps.set_plot(axs[0], initial, **kwds)]
        for axh, mat in zip(axs[1:], self.plast):
            imh.append(_ps.set_plot(axh, mat, **kwds))
        return imh

    @classmethod
    def build(cls, builder: Callable[..., Dict[str, la.lnarray]],
              nst: int, frac: _sb.ArrayLike = 0.5,
              extra_args=(), **kwargs) -> SynapseIdModel:
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
        *extra_args, **kwargs
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

    @classmethod
    def rand(cls, nst, *args, **kwargs) -> SynapseIdModel:
        """Random model

        Synapse model with random transition matrices

        Parameters
        ----------
            All passed to `cls.build`, `builders.build_rand`
            or `sl_py_tools.numpy_tricks.markov_param`:
        nst: int
            total number of states
        npl: int
            total number of plasticity types
        frac : float
            fraction of events that are potentiating, default=0.5.
        binary : bool
            is the weight vector binary? Otherwise it's linear. Default: False
        sp : float
            sparsity, default: 1.
        rng : np.random.Generator, optional
            random number generator, by default: builders.RNG
        ...
            extra arguments passed to `cls.build`, `builders.build_rand`
            or `sl_py_tools.numpy_tricks.markov...`.

        Returns
        -------
        synobj
            SynapseIdModel instance
        """
        return cls.build(build_rand, nst, *args, **kwargs)


# =============================================================================


def _weight_readout(weight: _sb.ArrayLike) -> la.lnarray:
    """Convert weight vector to readout vector"""
    if weight is None:
        return None
    return np.unique(weight, return_inverse=True)[1].astype(int)


def build_rand(nst: int, npl: int = 2, binary: bool = False,
               rng: np.random.Generator = _bld.RNG,
               **kwds) -> Dict[str, la.lnarray]:
    """Make a random model.

    Make a random model, i.e. transition rates are random numbers.
    For use with `SynapseMemoryModel.Build`.

    Parameters
    ----------
    n : int
        total number of states
    binary : bool
        is the weight vector binary? Otherwise it's linear. Default: False
    sparsity : float
        sparsity, default: 1.
    rng : np.random.Generator, optional
        random number generator, by default: builders.RNG
    optional arguments
        passed to `sl_py_tools.numpy_tricks.markov.rand_trans`.

    Returns
    -------
    dictionary : Dict[str, la.lnarray]
        plast : la.lnarray
            potentiation/depression transition matrices
        weight : la.lnarray
            synaptic weights (linear/binary)
        signal : la.lnarray
            desired signal contribution from each plasticity type
    """
    kwds.setdefault('rng', rng)
    return _bld.build_generic(ma.rand_trans_d, nst, npl, binary, **kwds)


def valid_shapes(model: SynapseIdModel) -> bool:
    """Do attributes (plast, weight, frac) have correct shapes?
    """
    vld = model.plast.shape[-2] == model.nstate
    vld &= model.initial.shape[-1] == model.nstate
    vld &= model.frac.shape[-1] == model.nplast
    vld &= model.readout.shape[-1] == model.nstate
    return vld


def valid_values(model: SynapseIdModel) -> bool:
    """Do attributes (plast, initial, frac) have valid values?
    """
    vld = np.isfinite(model.plast).all() and np.isfinite(model.initial).all()
    vld &= ma.isstochastic_d(model.plast, model.StochThresh)
    vld &= ma.isstochastic_d(model.initial, model.StochThresh)
    vld &= ma.isstochastic_d(model.frac, model.StochThresh)
    return vld


def well_behaved(model: SynapseIdModel, cond: bool = False) -> bool:
    """Are attributes (plast, initial) finite, and is Zinv well conditioned?

    Parameters
    ----------
    cond : bool, optional
        Do we check if Zinv is well conditioned? By default `False`.
    """
    vld = np.isfinite(model.plast).all() and np.isfinite(model.initial).all()
    if cond:
        markov = (model.frac.s * model.plast).sum(-3)
        fundi = np.ones_like(markov) + la.identity(model.nstate) - markov
        vld &= np.linalg.cond(fundi).max() < model.CondThresh
    return vld
