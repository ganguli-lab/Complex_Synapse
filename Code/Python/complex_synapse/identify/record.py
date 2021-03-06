"""Recording and playing back synapse fittinf
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

# import numpy as np

import numpy_linalg as la
# import sl_py_tools.arg_tricks as _ag
# import sl_py_tools.iter_tricks as _it
import sl_py_tools.containers as _cn

import complex_synapse.identify.plast_seq as _ps
import complex_synapse.identify.synapse_id as _si
import complex_synapse.identify.fit_synapse as _fs
import complex_synapse.synapse_base as _sb

_log = logging.getLogger(__name__)
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
    callback: Optional[_fs.Callback] = None
    step_num: int
    info: Dict[str, la.lnarray]
    est: Optional[la.lnarray]
    truth: Optional[la.lnarray]
    occ: Optional[la.lnarray]

    def __init__(self, callback: Optional[_fs.Callback] = None) -> None:
        self.callback, self.step_num = callback, 0
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
        nelm = _si.num_elements(fitter.est)
        self.est = la.empty((maxit, nelm), fitter.est.plast.dtype)
        # (T[,E],M) or (T[,E])
        occ = fitter.est_occ(slice(None))
        # (S,T[,E],M) or (S,T[,E])
        self.occ = la.empty((maxit,) + occ.shape, occ.dtype)
        if isinstance(fitter, _fs.GroundedFitter):
            self.truth = _si.elements(fitter.truth)
        for key, val in fitter.info.items():
            self.info[key] = la.empty((maxit,), type(val))

    def update(self, fitter: _fs.SynapseFitter) -> None:
        """Record one data point"""
        self.est[self.step_num] = _si.elements(fitter.est)
        self.occ[self.step_num] = fitter.est_occ(slice(None))
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
        self.update_info()

    def update_info(self) -> None:
        """Calculate stats for termination and display.
        """
        for key, val in self.saved.items():
            self.info[key] = val[self._ind]

    def update_fit(self) -> None:
        """Perform a single update of the model"""
        self.est = self.saved_est[self._ind]

    def est_occ(self, ind: _ps.Inds = slice(None)) -> la.lnarray:
        """Current estimate of state occupation

        Parameters
        ----------
        ind : Tuple[Union[int, slice], ...]
            Time, experiment indices/slices to plot

        Returns
        -------
        data : lnarray, (T,M) float[0:1] or (T,) int[0:M]
            Estimate of state occupation
        """
        # (T,M)
        return self.saved_occ[(self._ind,) + _cn.tuplify(ind)]

    def init(self) -> List[Any]:
        """Prepare for first iteration.

        Returns
        -------
        output : Any
            Whatever the callback returns.
        """
        _log.debug("Calling fitter.init")
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

        Returns
        -------
        output : Any
            Whatever the callback returns.
        """
        self._ind = step_num
        self.update_fit()
        self.update_info()
        if self._ind == self.opt.max_it - 1:
            self.info['result'] = int(self.check_thresh())
            _log.debug("Calling last fitter.step(%d)", step_num)
            return self.callback(self, 2)
        _log.debug("Calling fitter.step(%d)", step_num)
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
        self.saved = saved.copy()
        data = _ps.SimPlasticitySequence(self.saved['plast_type'],
                                         self.saved['readouts'],
                                         self.saved.pop('states'))
        truth = _si.from_elements(self.saved.pop('truth'),
                                  self.saved['frc'], self.saved['rdo'])
        kwds.setdefault('data', data)
        kwds.setdefault('truth', truth)
        super().__init__(saved=self.saved, callback=callback, **kwds)
