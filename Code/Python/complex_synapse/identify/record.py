"""Recording and playing back synapse fittinf
"""
from __future__ import annotations

from typing import Dict, Optional

# import numpy as np

import numpy_linalg as la
# import sl_py_tools.arg_tricks as _ag
import sl_py_tools.iter_tricks as _it

from . import plast_seq as _ps
from . import synapse_id as _si
from . import fit_synapse as _fs

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

    def __call__(self, fitter: _fs.SynapseFitter, pos: int) -> None:
        if pos == 0:
            self.setup(fitter)
        elif pos == 1 and fitter.info['nit'] % fitter.opt.disp_step == 0:
            self.update(fitter)
        elif pos == 2:
            self.cleanup(fitter)
        if self.callback is not None:
            self.callback(fitter, pos)

    def setup(self, fitter: _fs.SynapseFitter) -> None:
        """Setup the callback to record data"""
        self.step_num = 0
        maxit = fitter.opt.max_it // fitter.opt.disp_step
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

    Attributes
    ----------
    saved_est : la.lnarray (S, PM**2 + M)
        The sequence of estimated models.
    saved_occ : la.lnarray (S, M, T) or (S, T)
        The sequence of estimated paths.
    """
    saved: Dict[str, la.lnarray]
    saved_est: _si.SynapseIdModel
    saved_occ: la.lnarray
    _ind: int

    def __init__(self, data: _ps.PlasticitySequence,
                 saved: Dict[str, la.lnarray], **kwds) -> None:
        if not isinstance(self, GroundedFitterReplay):
            saved = saved.copy()

        self.saved_est = _si.from_elements(self.saved.pop('est'),
                                           self.saved.pop('frc'),
                                           self.saved.pop('rdo'))
        self.saved_occ = self.saved['occ']
        self._ind = 0
        super().__init__(data, self.saved_est[0], **kwds)

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

    def run(self, callback: _fs.Callback = _fs.print_callback) -> int:
        """Rerun the synapse fitter until it runs out

        Parameters
        ----------
        callback : Callable[[self, int]->None], optional
            Function called on every iteration, by default `print_callback`.
            Second argument:
                0: Before first iteration.
                1: During iterations.
                2: After completion.

        Returns
        -------
        result : int
            Flag for the outcome:
                -1: Error, invalid model generated.
                0: Failure, maximum iterations reached.
                1: Success, change in log-likelihood and model below threshold.
        """
        count = _it.undcount if self.opt.disp_each else _it.dcount
        callback(self, 0)
        self.info['result'] = 0
        for i in count('iteration', self.saved_est.nmodel[0]):
            self._ind = i
            self.update_fit()
            self.update_info()
            callback(self, 1)
        self.info['result'] = int(self.check_thresh())
        callback(self, 2)
        return self.info['result']


class GroundedFitterReplay(FitterReplay, _fs.GroundedFitter):
    """Class to replay the saved progress of another fitter with ground truth
    """

    def __init__(self, data: _ps.PlasticitySequence,
                 saved: Dict[str, la.lnarray], **kwds) -> None:
        self.saved = saved.copy()
        truth = _si.from_elements(self.saved.pop('truth'), self.saved['frc'],
                                  self.saved['rdo'])
        super().__init__(data=data, saved={}, truth=truth, **kwds)
