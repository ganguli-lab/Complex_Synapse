"""Model that is being fit to data
"""
from __future__ import annotations
from typing import Dict, List, Tuple, Union, Optional
import numpy as np
import numpy_linalg as la
from numpy_linalg import convert as cvl
from ..synapse_base import SynapseBase
from ..synapse_mem import ArrayLike
from sl_py_tools.numpy_tricks import markov as ma
from sl_py_tools.arg_tricks import default_non_eval as default


class SynapseIdModel(SynapseBase):
    """Model that is being fi to data
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
            self.plast += la.identity(nst, dtype)

    def dict_copy(self, keys=(), order='C', **kwds) -> Dict[str, la.lnarray]:
        """Dictionary with copies of data attributes.
        """
        keys += ('initial', 'readout')
        return super().dict_copy(keys=keys, order=order, **kwds)

    def __repr__(self) -> str:
        """Accurate representation of object"""
        rpr = super().__repr__()
        insert = f"    initial = {self.initial!r},\n"
        insert += f"    readout = {self.readout!r},\n"
        rpr = (rpr[:-1] + insert + rpr[-1])
        return rpr

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
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
        """Minimum number of readout values, R
        """
        return self.readout.max() + 1

# TODO: readout_proj, simulate, sort, plot
