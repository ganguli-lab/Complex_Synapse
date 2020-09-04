"""Recording and playing back synapse fittinf
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

# import numpy as np

import numpy_linalg as la
import sl_py_tools.arg_tricks as _ag
# import sl_py_tools.iter_tricks as _it

from . import plast_seq as _ps
from . import synapse_id as _si
from . import fit_synapse as _fs
from .. import synapse_base as _sb

# =============================================================================


class RecordingCallback:
    """A wrapper for other callbacks that records the progress of a fitter.

    Parameters
    ----------
    ind : Tuole[int|slice, ...]|None
        The snippet of data to store for describing path estimates.
    callback : Callable[[SynapseFitter, int] -> None]|None, optional
        A different callback being wrapped by this one.

    Attributes
    ----------
    info : Dict[str, lnarray]
        Dictionary carrying all of the recorded data.
    """
    ind: _ps.Inds
    callback: Optional[_fs.Callback] = None
    step_num: int
    info: Dict[str, la.lnarray]
    est: Optional[la.lnarray]
    truth: Optional[la.lnarray]
    occ: Optional[la.lnarray]

    def __init__(self, ind: Optional[_ps.Inds] = None,
                 callback: Optional[_fs.Callback] = None
                 ) -> None:
        self.ind, self.callback, self.step_num = ind, callback, 0
        self.info, self.est, self.truth, self.occ = {}, None, None, None

    def __call__(self, fitter: _fs.SynapseFitter, pos: int) -> List:
        if pos == 0:
            self.setup(fitter)
            self.update(fitter)
        elif pos == 1 and fitter.info['nit'] % fitter.opt.disp_step == 0:
            self.update(fitter)
        elif pos == 2:
            if fitter.info['nit'] % fitter.opt.disp_step:
                self.update(fitter)
            self.cleanup(fitter)
        if self.callback is not None:
            return self.callback(fitter, pos)
        return []

    def setup(self, fitter: _fs.SynapseFitter) -> None:
        """Setup the callback to record data"""
        self.step_num = 0
        maxit = fitter.opt.max_it // fitter.opt.disp_step + 2
        nelm = fitter.est.nplast * fitter.est.nstate**2 + fitter.est.nstate
        self.est = la.empty((maxit, nelm), fitter.est.plast.dtype)
        if self.ind is not None:
            occ = fitter.est_occ(self.ind)
            self.occ = la.empty((maxit,) + occ.shape, occ.dtype)
        if isinstance(fitter, _fs.GroundedFitter):
            self.truth = _si.elements(fitter.truth)
        for key, val in fitter.info.items():
            self.info[key] = la.empty((maxit,), type(val))

    def update(self, fitter: _fs.SynapseFitter) -> None:
        """Record one data point"""
        self.est[self.step_num] = _si.elements(fitter.est)
        if self.ind is not None:
            self.occ[self.step_num] = fitter.est_occ(self.ind)
        for key, val in fitter.info.items():
            if key != 'result':
                self.info[key][self.step_num] = val
        self.step_num += 1

    def cleanup(self, fitter: _fs.SynapseFitter) -> None:
        """Clean up recorded data"""
        # store final entry:
        if self.info['nit'][self.step_num-1] < fitter.info['nit']:
            self.update(fitter)
        # trim arrays
        for key, val in self.info.items():
            if key != 'result':
                self.info[key] = val[:self.step_num]
        # store other data in `info`
        self.info['occ'] = self.occ[:self.step_num]
        self.info['est'] = self.est[:self.step_num]
        if self.truth is not None:
            self.info['truth'] = self.truth
        self.info['frc'] = fitter.est.frac
        self.info['rdo'] = fitter.est.readout
        self.info.update(_sb.array_attrs(fitter.data))
        self.est, self.truth, self.occ = None, None, None


# =============================================================================
# Replaying fitters
# =============================================================================


class FitterReplay(_fs.SynapseFitter):
    """Class to replay the saved progress of another fitter

    Parameters
    ----------
    data : PlasticitySequence
        The data being fit to.
    saved : Dict[str, lnarray|Number]
        The saved data, the `info` attribute of a `RecordingCallback`.
    callback : Callable[[self, int]->None], optional
        Function called on every iteration, by default `print_callback`.
        Second argument:
            0: Before first iteration.
            1: During iterations.
            2: After completion.

    Attributes
    ----------
    saved_est : la.lnarray (S, PM**2 + M)
        The sequence of estimated models.
    saved_occ : la.lnarray (S, M, T) or (S, T)
        The sequence of estimated paths.
    """
    saved: Optional[Dict[str, la.lnarray]] = None
    saved_est: _si.SynapseIdModel
    saved_occ: la.lnarray
    _ind: int

    def __init__(self, saved: Dict[str, la.lnarray],
                 callback: _fs.Callback = _fs.print_callback, **kwds) -> None:
        self.saved = saved.copy()

        self.saved_est = _si.from_elements(self.saved.pop('est'),
                                           self.saved.pop('frc'),
                                           self.saved.pop('rdo'))
        self.saved_occ = self.saved.pop('occ')
        data = _ps.PlasticitySequence(self.saved.pop('plast_type'),
                                      self.saved.pop('readouts'))
        self._ind = 0
        kwds.setdefault('data', data)
        kwds.setdefault('est', self.saved_est[0])
        kwds.setdefault('max_it', self.saved_est.nmodel[0])
        kwds.setdefault('disp_step', 1)
        super().__init__(callback=callback, **kwds)

    def update_info(self) -> None:
        """Calculate stats for termination and display.
        """
        for key, val in self.saved.items():
            self.info[key] = val[self._ind]

    def update_fit(self) -> None:
        """Perform a single update of the model"""
        self.est = self.saved_est[self._ind]

    def est_occ(self, ind: _ps.Inds) -> la.lnarray:
        """Current estimate of state occupation

        Parameters
        ----------
        ind : Tuple[Union[int, slice], ...]
            Time, experiment indices/slices to plot

        Returns
        -------
        data : lnarray,  (T,M) float[0:1] or (T,) int[0:M]
            Estimate of state occupation
        """
        # (T.M)
        return self.saved_occ[self._ind]

    def init(self) -> List[Any]:
        """Prepare for first iteration.
        """
        self._ind = 0
        self.info['result'] = 0
        self.update_fit()
        self.update_info()
        return self.callback(self, 0)

    def step(self, step_num: int) -> List[Any]:
        """One update step

        Parameters
        ----------
        step_num : int
            Number of steps completed
        """
        self._ind = step_num
        self.update_fit()
        self.update_info()
        if self._ind == self.opt.max_it - 1:
            self.info['result'] = int(self.check_thresh())
            return self.callback(self, 2)
        return self.callback(self, 1)


class GroundedFitterReplay(FitterReplay, _fs.GroundedFitter):
    """Class to replay the saved progress of another fitter with ground truth

    Parameters
    ----------
    data : PlasticitySequence
        The data being fit to.
    saved : Dict[str, lnarray|Number]
        The saved data, the `info` attribute of a `RecordingCallback`.
    callback : Callable[[self, int]->None], optional
        Function called on every iteration, by default `print_callback`.
        Second argument:
            0: Before first iteration.
            1: During iterations.
            2: After completion.
    """

    def __init__(self, saved: Dict[str, la.lnarray],
                 callback: _fs.Callback = _fs.print_callback, **kwds) -> None:
        saved = saved.copy()
        data = _ps.SimPlasticitySequence(saved['plast_type'],
                                         saved['readouts'],
                                         saved.pop('states'))
        truth = _si.from_elements(saved.pop('truth'),
                                  saved['frc'], saved['rdo'])
        kwds.setdefault('data', data)
        kwds.setdefault('truth', truth)
        super().__init__(saved=saved, callback=callback, **kwds)
