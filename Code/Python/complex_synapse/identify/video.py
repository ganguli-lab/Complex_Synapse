"""Creating video frames for synapse fitting
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# for type hints
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colors import Normalize
from matplotlib.colorbar import Colorbar
from matplotlib.gridspec import GridSpec

from sl_py_tools.arg_tricks import default
import sl_py_tools.matplotlib_tricks as mpt

from .fit_synapse import SynapseFitter, GroundedFitter
from .plast_seq import Image, Line

# =============================================================================
# Helpers
# =============================================================================
# width: cbar=1, init=2, plast=12, plast_seq=*
# height: model=12, plast_type=2, readout=2, states=6
WIDTHS = {'c': 1, 'i': 2, 'p': 12}
HEIGHTS = {'pr': 2, 'st': 6, 'md': WIDTHS['p']}
SCL = 0.4


def _whr(npl: int = 2, row: int = 0) -> Dict[str, List[int]]:
    """width, height ratios"""
    return {'width_ratios': [WIDTHS[k] for k in ['i'] + ['p'] * npl + ['c']],
            'height_ratios': [HEIGHTS[k] for k in
                              ['md', 'pr', 'pr', 'st', 'md'][row:]]}


def _cbar(obj: FitterVideo, model: str) -> Colorbar:
    """Add a colourbar for a model"""
    return obj.fig.colorbar(obj.imh[model][0], cax=obj.axh[model][0],
                            label="Probability")


def fig_expt(npl: int = 2) -> Tuple[Figure, List[Axes], List[Axes]]:
    """Create a figure with axes for an experimental fit

    Parameters
    ----------
    npl : int, optional
        Number of types of plasticity, by default 2

    Returns
    -------
    fig : Figure
        The figure object for video frames.
    ps_ax : List[Axes]
        The axes objects for the `PlasticitySequence` data
    fit_ax : List[Axes]
        The axes objects for the `SynapseIdModel` fit.
    """
    # number of grid-ratio units
    l_x = (WIDTHS['c'] + WIDTHS['i']) + WIDTHS['p'] * npl
    l_y = HEIGHTS['md'] + 2 * HEIGHTS['pr'] + HEIGHTS['st']
    # number of grid units
    n_x, n_y = 2 + npl, 4
    # width: cbar=1, init=2, plast=12, plast_seq=*
    # height: model=12, plast_type=2, readout=2, states=6
    fig = plt.figure(figsize=(SCL * l_x, SCL * l_y), frameon=False,
                     constrained_layout=True)
    gsp = fig.add_gridspec(n_y, n_x, **_whr(npl, 1))
    # plast_seq, ht=4
    ps_ax = data_axes(fig, gsp, top=0)
    # fit_model, ht=12
    fit_ax = model_axes(fig, gsp, npl, row=-1)
    return fig, ps_ax, fit_ax


def fig_sim(npl: int = 2) -> Tuple[Figure, List[Axes], List[Axes], List[Axes]]:
    """Create a figure with axes for an simulation fit.

    Parameters
    ----------
    npl : int, optional
        Number of types of plasticity, by default 2

    Returns
    -------
    fig : Figure
        The figure object for video frames.
    ps_ax : List[Axes]
        The axes objects for the `PlasticitySequence` data
    fit_ax : List[Axes]
        The axes objects for the `SynapseIdModel` fit.
    true_ax : List[Axes]
        The axes objects for the ground-truth `SynapseIdModel`.
    """
    # number of grid-ratio units
    l_x = (WIDTHS['c'] + WIDTHS['i']) + WIDTHS['p'] * npl
    l_y = 2 * HEIGHTS['md'] + 2 * HEIGHTS['pr'] + HEIGHTS['st']
    # number of grid units
    n_x, n_y = 2 + npl, 5
    # width: cbar=1, init=2, plast=12, plast_seq=*
    # height: model=12, plast_type=2, readout=2, states=6
    fig = plt.figure(figsize=(SCL * l_x, SCL * l_y), frameon=False,
                     constrained_layout=True)
    gsp = fig.add_gridspec(n_y, n_x, **_whr(npl, 0))
    # true_model, ht=12
    true_ax = model_axes(fig, gsp, npl, row=0)
    # plast_seq, ht=10
    ps_ax = data_axes(fig, gsp, top=1)
    # fit_model, ht=12
    fit_ax = model_axes(fig, gsp, npl, row=-1)
    return fig, ps_ax, fit_ax, true_ax


def model_axes(fig: Figure, gsp: GridSpec, npl: int, row: int) -> List[Axes]:
    """Create axes for a synapse model.
    """
    cax = fig.add_subplot(gsp[row, -1])
    iax = fig.add_subplot(gsp[row, 0])
    pax = [fig.add_subplot(gsp[row, i+1], sharey=iax) for i in range(npl)]
    return [cax, iax] + pax


def data_axes(fig: Figure, gsp: GridSpec, top: int) -> List[Axes]:
    """Create axes for an experiment/simulation.
    """
    pax = fig.add_subplot(gsp[top, :])
    rax = fig.add_subplot(gsp[top+1, :], sharex=pax)
    sax = fig.add_subplot(gsp[top+2, :], sharex=pax)
    return [pax, rax, sax]


def label_model(axs: List[Axes], labels: List[str], **kwds):
    mpt.clean_axes(axs[0], **kwds)
    """Label axes for model"""
    axs[1].set_xticks([])
    axs[1].set_ylabel('Iniitial state')
    title = axs[1].set_title(' ', loc='left', in_layout=False, pad=20)
    mpt.clean_axes(axs[1], **kwds)
    for axh, lab in zip(axs[2:], labels[1:]):
        axh.set_ylabel("From state")
        axh.set_xlabel("To state")
        axh.set_title(lab)
        axh.xaxis.set_ticks_position('bottom')
        mpt.clean_axes(axh, **kwds)
    t_x, t_y = title.get_position()
    axs[1].text(t_x - 0.9, t_y + 0.08, labels[0], in_layout=False,
                fontsize=kwds.get('fontsize', 24), transform=axs[1].transAxes)


def label_data(axs: List[Axes], **kwds):
    """Label axes for plasticity sequence"""
    axs[0].set_yticks([])
    axs[1].set_yticks([])
    axs[2].set_ylabel("State")
    axs[2].set_xlabel("Time")
    for axh, lab in zip(axs, ["Plasticity type", "Synaptic efficacy",
                              "Occupation probability"]):
        axh.set_title(lab)
        axh.xaxis.set_ticks_position('bottom')
        axh.set_aspect('auto')
        mpt.clean_axes(axh, **kwds)


# =============================================================================
# Fitter video class
# =============================================================================


class FitterVideo:
    """Class to produce video frames showing fitter in action.
    """
    fitter: SynapseFitter
    ind: Inds
    norm: Normalize
    fig: Figure
    axh: Dict[str, List[Axes]]
    imh: Dict[str, List[Union[Image, Line, Colorbar]]]

    def __init__(self, fitter: SynapseFitter, inds: Inds,
                 norm: Optional[Normalize] = None, **kwds) -> None:
        self.fitter = fitter
        self.ind = inds
        self.norm = default(norm, Normalize(vmin=0., vmax=1.))

        self.axh = {}
        self.imh = {}
        npl = self.fitter.model.nplast
        if isinstance(self.fitter, GroundedFitter):
            (self.fig, self.axh['ps'],
             self.axh['fit'], self.axh['tr']) = fig_sim(npl)
        else:
            self.fig, self.axh['ps'], self.axh['fit'] = fig_expt(npl)
            self.axh['tr'] = []
        self.axh['st'] = self.axh['ps'][2]

        self.create_plots(**kwds)

    def __call__(self, fitter: SynapseFitter, pos: int) -> None:
        """Execute callback for fitter"""
        if pos == 0:
            self.create_plots()
        elif pos == 1:
            self.update_plots()

    def create_plots(self, **kwds) -> None:
        """Create initial plots"""
        plast_labs = kwds.pop('plast_labs', ["Potentiation", "Depression"])
        model_labs = kwds.pop('model_labs', ["Fit model", "True model"])
        opt = {'norm': self.norm, 'zorder': 0}
        fmt = {'box': False, 'tight': False}

        self.imh['st'] = self.fitter.plot_occ(self.axh['st'], self.ind, **opt)
        opt['zorder'] = 10
        self.imh['ps'] = self.fitter.data[self.ind].plot(self.axh['ps'])
        label_data(self.axh['ps'], **fmt)

        self.imh['fit'] = self.fitter.model.plot(self.axh['fit'][1:], **opt)
        self.imh['cfit'] = _cbar(self, 'fit')
        label_model(self.axh['fit'], model_labs[:1] + plast_labs, **fmt)

        if self.axh['tr']:
            self.imh['tr'] = self.fitter.truth.plot(self.axh['tr'][1:],  **opt)
            self.imh['ctr'] = _cbar(self, 'tr')
            label_model(self.axh['tr'], model_labs[1:] + plast_labs, **fmt)

        plt.draw()

    def update_plots(self) -> None:
        """Update plots after iteration"""
        self.imh['st'] = self.fitter.plot_occ(self.imh['st'], self.ind)
        self.imh['fit'] = self.fitter.model.plot(self.imh['fit'])
        plt.draw()


# =============================================================================
# Hint valiases
# =============================================================================
Inds = Tuple[Union[int, slice], ...]
