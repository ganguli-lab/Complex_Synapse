"""Creating video frames for synapse fitting
"""
from __future__ import annotations

from typing import ChainMap, Dict, List, Tuple, Union
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import sl_py_tools.containers as cn
import sl_py_tools.matplotlib_tricks as mpt

from .fit_synapse import SynapseFitter, GroundedFitter, print_callback
from .plast_seq import Image, Line
Normalize = mpl.colors.Normalize

# =============================================================================
# Helpers
# =============================================================================
mpt.rc_colours()
mpt.rc_fonts('sans-serif')


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
    pax = fig.add_subplot(gsp[top, :-1])
    rax = fig.add_subplot(gsp[top+1, :-1], sharex=pax)
    sax = fig.add_subplot(gsp[top+2, :-1], sharex=pax)
    cpax = fig.add_subplot(gsp[top, -1])
    crax = fig.add_subplot(gsp[top+1, -1])
    csax = fig.add_subplot(gsp[top+2, -1])
    return [pax, rax, sax, cpax, crax, csax]


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


def label_data(axs: List[Axes], cax: List[Axes], **kwds):
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
    for axh in cax:
        mpt.clean_axes(axh, **kwds)


# =============================================================================
# Fitter video options class
# =============================================================================


class VideoOptions:
    """Visual options for FitterVideo

    Parameters
    ----------
    plast_type : List[str] = ["dep.", "pot."]
        Colourbar labels for `fitter.data.plasticity_type`.
    readout : List[str] = ["weak", "strong"]
        Colourbar labels for `fitter.data.readout`.
    model : List[str] = ["Fit, $\\mathbf{M}$", "Truth, $\\mathbf{M}^*$"]
        Overall labels for `fitter.model` and `fitter.truth`..
    plast : List[str] = ["Potentiation", "Depression"]
        Transition matrix labels.
    txt_opt : Dict[str, bool] = {'box': False, 'tight': False}
         Options for text. See `sl_py_tools.matplotlib_tricks.clean_axes`.
    im_opt : Dict[str, str] = {'cmap': 'YlOrBr'}
        Options for image plots.
    All of them are keyword only. Unknown keywords are passed to `txt_opt` if
    `sl_py_tools.matplotlib_tricks.clean_axes` takes them, or `im_opt` if not.
    """
    # Text for labels
    plast_type: List[str]
    readout: List[str]
    model: List[str]
    plast: List[str]
    # keyword options
    txt_opt: Dict[str, bool]
    im_opt: Dict[str, str]

    def __init__(self, **kwds) -> None:
        self.plast_type = kwds.pop('plast_type', ["dep.", "pot."])
        self.readout = kwds.pop('readout', ["weak", "strong"])
        self.model = kwds.pop('model',
                              ["Fit, $\\mathbf{M}$", "Truth, $\\mathbf{M}^*$"])
        self.plast = kwds.pop('plast', ["Potentiation", "Depression"])
        self.txt_opt = kwds.pop('txt_opt', {'box': False, 'tight': False})
        self.im_opt = kwds.pop('im_opt', {'cmap': 'YlOrBr'})
        for key, val in ChainMap(self.txt_opt, self.im_opt).items():
            kwds.setdefault(key, val)
        self.txt_opt.update(mpt.clean_axes_keys(kwds))
        self.im_opt.update(kwds)

    def model_labs(self, ind: Inds) -> List[str]:
        """Labels for a model's heatmaps"""
        return cn.listify(self.model[ind]) + self.plast

    def __repr__(self) -> str:
        return (type(self).__name__ + "plast_type={self.plast_type}, "
                + "readout={self.readout}, model={self.model}, "
                + "plast={self.plast}, txt_opt={self.txt_opt}, "
                + "im_opt={self.im_opt}, ")


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
    fname : str, optional
        Strings for filenames - produced by `fname.format(frams_number)`,
        by default `""`: don't save files.
    norm : Optional[Normalize], optional
        Data range for probability heatmaps, by default `Normalize(0, 1)`.
    """
    ind: Inds
    fname: str
    norm: Normalize
    fig: Figure
    axh: Dict[str, List[Axes]]
    imh: Dict[str, List[Union[Image, Line, Colorbar]]]
    opt: VideoOptions

    def __init__(self, fitter: SynapseFitter, inds: Inds, fname: str = "",
                 norm: Normalize = Normalize(0., 1.), **kwds) -> None:
        """Video frame producing callback for SynapseFitter.

        Parameters
        ----------
        fitter : SynapseFitter
            Object that performs model fit
        inds : Inds
            Indices for `fitter.data` for snippet to plot. Must result in
            `fitter.data[ind].nexpt = ()`.
        fname : str, optional
            Strings for filenames - produced by `fname.format(frams_number)`,
            by default `""`: don't save files.
        norm : Optional[Normalize], optional
            Data range for probability heatmaps, by default `Normalize(0, 1)`.
        Other keywords passed to `self.opt`
        """
        self.ind = inds
        self.fname = fname
        self.norm = norm
        self.opt = VideoOptions(**kwds)

        self.axh = {}
        self.imh = {}
        npl = fitter.model.nplast
        if isinstance(fitter, GroundedFitter):
            (self.fig, self.axh['ps'],
             self.axh['fit'], self.axh['tr']) = fig_sim(npl)
        else:
            self.fig, self.axh['ps'], self.axh['fit'] = fig_expt(npl)
            self.axh['tr'] = []
        self.axh['cps'] = self.axh['ps'][3:]
        self.axh['ps'] = self.axh['ps'][:3]
        self.axh['st'] = self.axh['ps'][2]
        self.axh['cps'][2].set(frame_on=False, xticks=[], yticks=[])

    def __call__(self, fitter: SynapseFitter, pos: int) -> None:
        """Execute callback for fitter"""
        if pos == 0:
            self.create_plots(fitter)
            print_callback(fitter, 0)
        elif pos == 1:
            if fitter.stats['nit'] % fitter.opt.disp_step == 0:
                self.update_plots(fitter)
                self.savefig(fitter.stats['nit'])
        elif pos == 2:
            print_callback(fitter, 2)

    def create_plots(self, fitter: SynapseFitter) -> None:
        """Create initial plots"""
        fit_lab, tru_lab = self.opt.model_labs(0), self.opt.model_labs(1)
        mdo = {**self.opt.im_opt, 'zorder': 0, 'norm': self.norm}
        pso = {**self.opt.im_opt, 'zorder': 10, 'nplast': fitter.model.nplast,
               'nreadout': fitter.model.nreadout}

        self.imh['st'] = fitter.plot_occ(self.axh['st'], self.ind, **mdo)
        self.imh['ps'] = fitter.data[self.ind].plot(self.axh['ps'], **pso)
        self.imh['cpt'] = [_cbarp(self, 0, self.opt.plast_type),
                           _cbarp(self, 1, self.opt.readout)]
        label_data(self.axh['ps'], self.axh['cps'][:2], **self.opt.txt_opt)

        self.imh['fit'] = fitter.model.plot(self.axh['fit'][1:], **mdo)
        self.imh['cfit'] = _cbarm(self, 'fit')
        label_model(self.axh['fit'], fit_lab, **self.opt.txt_opt)

        if self.axh['tr']:
            self.imh['tr'] = fitter.truth.plot(self.axh['tr'][1:], **mdo)
            self.imh['ctr'] = _cbarm(self, 'tr')
            label_model(self.axh['tr'], tru_lab, **self.opt.txt_opt)

        self.write_stats(fitter)
        # self.fig.canvas.draw_idle()
        plt.draw()

    def update_plots(self, fitter: SynapseFitter) -> None:
        """Update plots after iteration"""
        self.imh['st'] = fitter.plot_occ(self.imh['st'], self.ind)
        self.imh['fit'] = fitter.model.plot(self.imh['fit'])
        self.write_stats(fitter)
        plt.draw()
        # self.fig.canvas.draw_idle()

    def savefig(self, fileno: Union[None, int, str]) -> None:
        """Save current figure as a file"""
        if self.fname:
            self.fig.savefig(self.fname.format(fileno))

    def write_stats(self, fitter: SynapseFitter) -> None:
        """Write stats next to states' axes"""
        txt = format(fitter, 'tex')
        if 'stats' in self.imh:
            self.imh['stats'].set_text(txt)
        else:
            axs = self.axh['cps'][2]
            self.imh['stats'] = axs.text(5, 1, txt, va='top', ha='right',
                                         transform=axs.transAxes)


# =============================================================================
# Helpers
# =============================================================================


def _gratios(npl: int, row: int) -> Tuple[List[int], List[int]]:
    """width, height ratios"""
    # width: cbar=1, init=2, plast=12, plast_seq=*
    # height: model=12, plast_type=2, readout=2, states=6
    widths = {'c': 1, 'i': 2, 'p': 12}
    heights = {'pr': 2, 'st': 6, 'md': widths['p']}
    return ([heights[k] for k in ['md', 'pr', 'pr', 'st', 'md'][row:]],
            [widths[k] for k in ['i'] + ['p'] * npl + ['c']])


def _make_fig(npl: int = 2, row: int = 0) -> Tuple[Figure, GridSpec]:
    """Create a figure and grid spec"""
    scl = 0.3
    keys = ['height_ratios', 'width_ratios']

    grs = _gratios(npl, row)
    fsiz = tuple(scl * sum(g) for g in reversed(grs))
    args, kwargs = tuple(map(len, grs)), dict(zip(keys, grs))

    fig = plt.figure(figsize=fsiz, frameon=True, constrained_layout=True)
    gsp = fig.add_gridspec(*args, **kwargs)
    return fig, gsp


def _cbarm(obj: FitterVideo, model: str) -> Colorbar:
    """Add a colourbar for a model"""
    return obj.fig.colorbar(obj.imh[model][0], cax=obj.axh[model][0],
                            label="Probability")


def _cbarp(obj: FitterVideo, ind: int, keys: List[str]) -> Colorbar:
    """Add a colourbar for a plasticity sequence variable"""
    cbh = obj.fig.colorbar(obj.imh['ps'][ind], cax=obj.axh['cps'][ind],
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
Colorbar = mpl.colorbar.Colorbar
GridSpec = mpl.gridspec.GridSpec
