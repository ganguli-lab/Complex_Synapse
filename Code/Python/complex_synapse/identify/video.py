"""Creating video frames for synapse fitting
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from sl_py_tools.arg_tricks import default
import sl_py_tools.matplotlib_tricks as mpt

from .fit_synapse import SynapseFitter, GroundedFitter, print_callback
from .plast_seq import Image, Line

# =============================================================================
# Helpers
# =============================================================================
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# width: cbar=1, init=2, plast=12, plast_seq=*
# height: model=12, plast_type=2, readout=2, states=6
WIDTHS = {'c': 1, 'i': 2, 'p': 12}
HEIGHTS = {'pr': 2, 'st': 6, 'md': WIDTHS['p']}
SCL = 0.4


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
    fig, gsp = _make_fig(npl, row=1)
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
    fig, gsp = _make_fig(npl, row=0)
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
    """Label axes for model"""
    mpt.clean_axes(axs[0], **kwds)
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
    size = kwds.get('titlefontsize', kwds.get('fontsize', 24))
    axs[1].text(t_x - 0.9, t_y + 0.06, f"\\textit{{{labels[0]}}}",
                in_layout=False, transform=axs[1].transAxes, fontsize=size)
                # bbox={'boxstyle': 'Square', 'facecolor': 'none'})


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

    Parameters
    ----------
    fitter : SynapseFitter
        Object that performs model fit
    inds : Inds
        Indices for `fitter.data` for snippet to plot. Must result in
        `fitter.data[ind].nexpt = ()`.
    fname : Optional[str], optional
        Strings for filenames - produced by `fname.format(frams_number)`,
        by default `None`: don't save files.
    norm : Optional[Normalize], optional
        Data range for probability heatmaps, by default `Normalize(0, 1)`.
    """
    fitter: SynapseFitter
    ind: Inds
    fname: Optional[str]
    norm: Normalize
    fig: Figure
    axh: Dict[str, List[Axes]]
    imh: Dict[str, List[Union[Image, Line, Colorbar]]]

    def __init__(self, fitter: SynapseFitter, inds: Inds,
                 fname: Optional[str] = None,
                 norm: Optional[Normalize] = None, **kwds) -> None:
        """Class to produce video frames showing fitter in action.

        Parameters
        ----------
        fitter : SynapseFitter
            Object that performs model fit
        inds : Inds
            Indices for `fitter.data` for snippet to plot. Must result in
            `fitter.data[ind].nexpt = ()`.
        fname : Optional[str], optional
            Strings for filenames - produced by `fname.format(frams_number)`,
            by default `None`: don't save files.
        norm : Optional[Normalize], optional
            Data range for probability heatmaps, by default `Normalize(0, 1)`.
        """
        self.fitter = fitter
        self.ind = inds
        self.norm = default(norm, Normalize(vmin=0., vmax=1.))
        self.fname = fname

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

        # self.create_plots(**kwds)

    def __call__(self, fitter: SynapseFitter, pos: int) -> None:
        """Execute callback for fitter"""
        if pos == 0:
            self.create_plots()
            print_callback(fitter, 0)
        elif pos == 1:
            if fitter.stats['nit'] % fitter.opt['disp_step'] == 0:
                self.update_plots()
                self.savefig(self.fitter.stats['nit'])
        elif pos == 2:
            print_callback(fitter, 2)

    def create_plots(self, **kwds) -> None:
        """Create initial plots"""
        pty_lab = kwds.pop('pty_labels', ["dep.", "pot."])
        rdo_lab = kwds.pop('rdo_labels', ["weak", "strong"])
        model_lab = kwds.pop('model_labels', ["Fit model", "True model"])
        plast_lab = kwds.pop('plast_labels', ["Potentiation", "Depression"])
        opt = {'norm': self.norm, 'zorder': 0}
        pso = {'zorder': 10, 'cmap': 'YlOrBr'}
        fmt = {'box': False, 'tight': False, **kwds}

        self.imh['st'] = self.fitter.plot_occ(self.axh['st'], self.ind, **opt)
        self.imh['ps'] = self.fitter.data[self.ind].plot(self.axh['ps'], **pso)
        label_data(self.axh['ps'], **fmt)
        self.imh['cpt'] = [_cbarp(self, 0, pty_lab), _cbarp(self, 1, rdo_lab)]
        opt.pop('zorder', None)

        self.imh['fit'] = self.fitter.model.plot(self.axh['fit'][1:], **opt)
        self.imh['cfit'] = _cbarm(self, 'fit')
        label_model(self.axh['fit'], model_lab[:1] + plast_lab, **fmt)

        if self.axh['tr']:
            self.imh['tr'] = self.fitter.truth.plot(self.axh['tr'][1:],  **opt)
            self.imh['ctr'] = _cbarm(self, 'tr')
            label_model(self.axh['tr'], model_lab[1:] + plast_lab, **fmt)

        # self.fig.canvas.draw_idle()
        plt.draw()

    def update_plots(self) -> None:
        """Update plots after iteration"""
        self.imh['st'] = self.fitter.plot_occ(self.imh['st'], self.ind)
        self.imh['fit'] = self.fitter.model.plot(self.imh['fit'])
        plt.draw()
        # self.fig.canvas.draw_idle()

    def savefig(self, fileno: Union[None, int, str]) -> None:
        """Save current figure as a file"""
        if self.fname:
            self.fig.savefig(self.fname.format(fileno))


# =============================================================================
# Helpers
# =============================================================================


def _gratios(npl: int, row: int) -> Tuple[List[int], List[int]]:
    """width, height ratios"""
    # width: cbar=1, init=2, plast=12, plast_seq=*
    # height: model=12, plast_type=2, readout=2, states=6
    return ([HEIGHTS[k] for k in ['md', 'pr', 'pr', 'st', 'md'][row:]],
            [WIDTHS[k] for k in ['i'] + ['p'] * npl + ['c']])


def _gsarg(npl: int, row: int) -> Tuple[Tuple[int, int], Dict[str, List[int]]]:
    """grid spec arguments"""
    # width: cbar=1, init=2, plast=12, plast_seq=*
    # height: model=12, plast_type=2, readout=2, states=6
    grs = _gratios(npl, row)
    keys = ['height_ratios', 'width_ratios']
    return tuple(map(len, grs)), dict(zip(keys, grs))


def _fsiz(npl: int, row: int) -> Tuple[float, float]:
    """width, height of figure"""
    # width: cbar=1, init=2, plast=12, plast_seq=*
    # height: model=12, plast_type=2, readout=2, states=6
    grs = _gratios(npl, row)
    return tuple(SCL * sum(g) for g in reversed(grs))


def _make_fig(npl: int = 2, row: int = 0) -> Tuple[Figure, GridSpec]:
    """Create a figure and grid spec"""
    fig = plt.figure(figsize=_fsiz(npl, row), frameon=True,
                     constrained_layout=True)
    args, kwargs = _gsarg(npl, row)
    gsp = fig.add_gridspec(*args, **kwargs)
    return fig, gsp


def _cbarm(obj: FitterVideo, model: str) -> Colorbar:
    """Add a colourbar for a model"""
    return obj.fig.colorbar(obj.imh[model][0], cax=obj.axh[model][0],
                            label="Probability")


def _cbarp(obj: FitterVideo, ind: int, keys: List[str]) -> Colorbar:
    """Add a colourbar for a plasticity sequence"""
    cbh = obj.fig.colorbar(obj.imh['ps'][ind], ax=obj.axh['ps'][ind],
                           drawedges=True, fraction=0.3, aspect=5)
    cbh.set_ticks(np.arange(len(keys)))
    cbh.set_ticklabels(keys)
    return cbh


# =============================================================================
# Hint valiases
# =============================================================================
Inds = Tuple[Union[int, slice], ...]
Figure = mpl.figure.Figure
Axes = mpl.axes.Axes
Normalize = mpl.colors.Normalize
Colorbar = mpl.colorbar.Colorbar
GridSpec = mpl.gridspec.GridSpec
