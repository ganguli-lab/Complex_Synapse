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

# TODO: readout_proj, simulate, sort, plot
