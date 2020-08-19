"""Creating video frames for synapse fitting
"""
from __future__ import annotations

from typing import Any, ChainMap, Dict, List, Tuple, Union
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import sl_py_tools.containers as cn
import sl_py_tools.matplotlib_tricks as mpt

from .fit_synapse import SynapseFitter, GroundedFitter, print_callback
from .plast_seq import Image, Line

Normalize = mpl.colors.Normalize
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# =============================================================================
# Fitter video options class
# =============================================================================


class VideoOptions:
    """Visual options for FitterVideo

    Parameters
    ----------
    plast_type : List[str] = ["Plasticity type", "dep.", "pot."]
        Title and colourbar labels for `fitter.data.plasticity_type`.
    readout : List[str] = ["Synaptic efficacy", "weak", "strong"]
        Title and colourbar labels for `fitter.data.readout`.
    state : List[str] = ["Occupation probability", "Time", "State"]
        Title and axes labels for `fitter.data.state`, `fitter.model.plot_occ`.
    model : List[str] = ["Fit, $\\mathbf{M}$", "Truth, $\\mathbf{M}^*$"]
        Overall labels for `fitter.model` and `fitter.truth`..
    plast : List[str] = ["Potentiation", "Depression"]
        Transition matrix labels.
    transpose : bool = False
        Transpose the layout of the video? Stored in `txt_opt` and `im_opt`
        under the key `'trn'`
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
    state: List[str]
    model: List[str]
    plast: List[str]
    # keyword options
    txt_opt: Dict[str, bool]
    im_opt: Dict[str, str]

    def __init__(self, **kwds) -> None:
        self.plast_type = kwds.pop('plast_type', ["Plasticity type",
                                                  "dep.", "pot."])
        self.readout = kwds.pop('readout', ["Synaptic efficacy",
                                            "weak", "strong"])
        self.state = kwds.pop('state', ["Occupation probability",
                                        "Time", "State"])
        self.model = kwds.pop('model',
                              ["Fit, $\\mathbf{M}$", "Truth, $\\mathbf{M}^*$"])
        self.plast = kwds.pop('plast', ["Potentiation", "Depression"])

        transpose = kwds.pop('transpose', False)
        self.txt_opt = kwds.pop('txt_opt', {'box': False, 'tight': False})
        self.im_opt = kwds.pop('im_opt', {'cmap': 'YlOrBr'})
        for key, val in ChainMap(self.txt_opt, self.im_opt).items():
            kwds.setdefault(key, val)
        self.txt_opt.update(mpt.clean_axes_keys(kwds))
        self.im_opt.update(kwds)
        if transpose:
            self.txt_opt['trn'] = True
            self.im_opt['trn'] = True

    def model_labs(self, ind: Inds) -> List[str]:
        """Labels for a model's heatmaps"""
        return cn.listify(self.model[ind]) + self.plast

    def update(self, opt: cn.Dictable[str, Any], **kwds) -> None:
        """Update options."""
        opt = dict(opt, **kwds)
        for key, val in opt.items():
            if key in {'txt_opt', 'im_opt'}:
                getattr(self, key).update(val)
            else:
                setattr(self, key, val)

    @property
    def transpose(self) -> bool:
        """Transpose the layout of the video?"""
        return self.txt_opt.get('trn', False)

    @transpose.setter
    def transpose(self, value: bool):
        """Transpose the layout of the video?"""
        self.txt_opt['trn'] = value
        self.im_opt['trn'] = value

    @transpose.deleter
    def transpose(self):
        """Transpose the layout of the video?"""
        self.txt_opt.pop('trn', None)
        self.im_opt.pop('trn', None)

    def __repr__(self) -> str:
        return type(self).__name__ + f"(**{self.__dict__})"


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
    imh: Dict[str, List[Union[Image, Line, mpl.text.Text]]]
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
        self.opt = kwds.pop('opt', VideoOptions())
        self.opt.update(kwds)

        self.axh = {}
        self.imh = {}
        fopt = (fitter.model.nplast, self.opt.transpose)
        if isinstance(fitter, GroundedFitter):
            (self.fig, psax,
             self.axh['fit'], self.axh['tr']) = fig_sim(*fopt)
        else:
            self.fig, psax, self.axh['fit'] = fig_expt(*fopt)
            self.axh['tr'] = []
        self.axh['ps'] = psax[1::2]
        self.axh['pt'] = psax[:2]
        self.axh['ro'] = psax[2:4]
        self.axh['st'] = psax[4:]
        self.axh['info'] = self.axh['tr'][0] if self.axh['tr'] else psax[4]
        # self.axh['info'].set(frame_on=False, xticks=[], yticks=[])

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
        legend = bool(self.axh['tr'])
        lbo = self.opt.txt_opt
        mdo = {**self.opt.im_opt, 'zorder': 0, 'norm': self.norm}
        pso = {**self.opt.im_opt, 'zorder': 10, 'nplast': fitter.model.nplast,
               'nreadout': fitter.model.nreadout}

        self.imh['st'] = fitter.plot_occ(self.axh['st'][1], self.ind, **mdo)
        self.imh['ps'] = fitter.data[self.ind].plot(self.axh['ps'], **pso)
        label_data(*self._fai('pt'), self.opt.plast_type, **lbo)
        label_data(*self._fai('ro'), self.opt.readout, **lbo)
        label_sim(*self._fai('st')[1:], self.opt.state, legend=legend, **lbo)

        self.imh['fit'] = fitter.model.plot(self.axh['fit'][1:], **mdo)
        label_model(*self._fai('fit'), fit_lab, **lbo)

        if self.axh['tr']:
            self.imh['tr'] = fitter.truth.plot(self.axh['tr'][1:], **mdo)
            label_model(*self._fai('tr'), tru_lab, cbar=False, **lbo)

        self.write_stats(fitter)
        # self.fig.canvas.draw_idle()
        plt.draw()

    def update_plots(self, fitter: SynapseFitter) -> None:
        """Update plots after iteration"""
        trn = self.opt.transpose
        self.imh['st'] = fitter.plot_occ(self.imh['st'], self.ind, trn=trn)
        self.imh['fit'] = fitter.model.plot(self.imh['fit'], trn=trn)
        self.write_stats(fitter)
        plt.draw()
        # self.fig.canvas.draw_idle()

    def savefig(self, fileno: Union[None, int, str]) -> None:
        """Save current figure as a file"""
        if self.fname:
            self.fig.savefig(self.fname.format(fileno))

    def write_stats(self, fitter: SynapseFitter) -> None:
        """Write stats next to states' axes"""
        self.imh['stats'] = write_stats(
            format(fitter, 'tex'), self.imh.get('stats', self.axh['info']),
            self.opt.transpose)

    def _fai(self, key: str) -> Tuple[Figure, List[Axes], Image]:
        """Get figure, axes, image handles"""
        ps_keys = ('pt', 'ro', 'st')
        if key in ps_keys:
            return self.fig, self.axh[key], self.imh['ps'][ps_keys.index(key)]
        return self.fig, self.axh[key], self.imh[key][0]


# =============================================================================
# Helpers
# =============================================================================


def fig_expt(npl: int = 2, trn: bool = False) -> Tuple[Figure, List[Axes],
                                                       List[Axes]]:
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
    fig, gsp = _make_fig(npl, row=1, trn=trn)
    # plast_seq, ht=4
    ps_ax = _data_axes(fig, gsp, top=0, trn=trn)
    # fit_model, ht=12
    fit_ax = _model_axes(fig, gsp, npl, row=-1, trn=trn)
    return fig, ps_ax, fit_ax


def fig_sim(npl: int = 2, trn: bool = False) -> Tuple[Figure, List[Axes],
                                                      List[Axes], List[Axes]]:
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
    fig, gsp = _make_fig(npl, row=0, trn=trn)
    # true_model, ht=12
    true_ax = _model_axes(fig, gsp, npl, row=0, trn=trn)
    # plast_seq, ht=10
    ps_ax = _data_axes(fig, gsp, top=1, trn=trn)
    # fit_model, ht=12
    fit_ax = _model_axes(fig, gsp, npl, row=-1, trn=trn)
    return fig, ps_ax, fit_ax, true_ax


def label_model(fig: Figure, axs: List[Axes], imh: Image, labels: List[str],
                **kwds):
    """Label axes for model"""
    trn = kwds.pop('trn', False)
    if kwds.pop('cbar', True):
        cbopt = {'orientation': 'horizontal'} if trn else {}
        fig.colorbar(imh, cax=axs[0], label="Probability", **cbopt)
        mpt.clean_axes(axs[0], **kwds)
    else:
        axs[0].set(frame_on=False, xticks=[], yticks=[])

    if trn:
        axs[1].xaxis.set_ticks_position('bottom')
    to_set = ('yticks', 'xlabel') if trn else ('xticks', 'ylabel')
    axs[1].set(**{to_set[0]: [], to_set[1]: 'Initial state'})
    size = kwds.get('titlefontsize', kwds.get('fontsize', 24))
    tlab = f"\\textit{{{labels[0]}}}" if trn else ' '
    title = axs[1].set_title(tlab, loc='left', pad=20, fontsize=size)
    mpt.clean_axes(axs[1], **kwds)

    for axh, lab in zip(axs[2:], labels[1:]):
        axh.set_ylabel("From state")
        axh.set_xlabel("To state")
        axh.set_title(lab)
        axh.xaxis.set_ticks_position('bottom')
        mpt.clean_axes(axh, **kwds)

    if not trn:
        t_x, t_y = title.get_position()
        axs[1].text(t_x - 0.9, t_y + 0.06, f"\\textit{{{labels[0]}}}",
                    in_layout=False, transform=axs[1].transAxes, fontsize=size)
                    # bbox={'boxstyle': 'Square', 'facecolor': 'none'})


def label_data(fig: Figure, axs: List[Axes], imh: Image, labels: List[str],
               **kwds):
    """Label axes for plasticity sequence"""
    trn = kwds.pop('trn', False)
    cbopt = {'orientation': 'horizontal'} if trn else {}
    cbh = fig.colorbar(imh, cax=axs[0], drawedges=True, aspect=5, **cbopt)
    cbh.set_ticks(np.arange(len(labels)-1))
    cbh.set_ticklabels(labels[1:])
    if trn:
        axs[0].tick_params(labelrotation=45)
    mpt.clean_axes(axs[0], **kwds)

    title = labels[0].replace(" ", "\n") if trn else labels[0]
    axs[1].set_title(title)
    if trn:
        axs[1].set_xticks([])
        axs[1].yaxis.set_ticks_position('right')
    else:
        axs[1].set_yticks([])
        axs[1].xaxis.set_ticks_position('bottom')
    axs[1].set_aspect('auto')
    mpt.clean_axes(axs[1], **kwds)


def label_sim(axs: List[Axes], lnh: Line, labels: List[str], **kwds):
    """Label axes for plasticity sequence"""
    trn = kwds.pop('trn', False)
    if kwds.pop('legend', False):
        axs[0].legend(handles=[lnh], loc=10)
        mpt.clean_axes(axs[0], **kwds)
    axs[0].set(frame_on=False, xticks=[], yticks=[])
    axlabels = labels[:0:-1] if trn else labels[1:]
    axs[1].set_xlabel(axlabels[0])
    axs[1].set_ylabel(axlabels[1])
    title = labels[0].replace(" ", "\n") if trn else labels[0]
    axs[1].set_title(title)
    axs[1].xaxis.set_ticks_position('bottom')
    if trn:
        axs[1].yaxis.set_label_position('right')
        axs[1].yaxis.set_ticks_position('right')
    axs[1].set_aspect('auto')
    mpt.clean_axes(axs[1], **kwds)


def write_stats(txt: str, handle: Union[Axes, Text], trn: bool = False
                ) -> Text:
    """Write stats next to states' axes"""
    # if trn:
    #     txt = _rep_every_other(txt, '\\\\', '&')
    if isinstance(handle, Text):
        handle.set_text(txt)
        return handle
    pos, yal = (3, 'top') if trn else (0, 'bottom')
    return handle.text(0.5, pos, txt, va=yal, ha='center',
                       transform=handle.transAxes)

# =============================================================================
# Private helpers
# =============================================================================


def _make_fig(npl: int, row: int, trn: bool) -> Tuple[Figure, GridSpec]:
    """Create a figure and grid spec"""
    scl = 0.3
    keys = ['height_ratios', 'width_ratios']

    grs = _gratios(npl, row, trn)
    fsiz = tuple(scl * sum(g) for g in reversed(grs))
    args, kwargs = tuple(map(len, grs)), dict(zip(keys, grs))

    fig = plt.figure(figsize=fsiz, frameon=True, constrained_layout=True)
    gsp = fig.add_gridspec(*args, **kwargs)
    gsp = TransposeGridSpec(gsp) if trn else gsp
    return fig, gsp


def _gratios(npl: int, row: int, trn: bool) -> Tuple[List[int], List[int]]:
    """width, height ratios"""
    # width: cbar=1, init=2, plast=12, plast_seq=*
    # height: model=12, plast_type=2, readout=2, states=6
    widths = {'c': 1, 'i': 2, 'p': 12}
    heights = {'pr': 2, 'st': 6, 'md': widths['p']}
    rats = ([heights[k] for k in ['md', 'pr', 'pr', 'st', 'md'][row:]],
            [widths[k] for k in ['i'] + ['p'] * npl + ['c']])
    return rats[::-1] if trn else rats


def _model_axes(fig: Figure, gsp: GridSpec, npl: int, row: int, trn: bool
                ) -> List[Axes]:
    """Create axes for a synapse model.
    """
    cax = fig.add_subplot(gsp[row, -1])
    iax = fig.add_subplot(gsp[row, 0])
    share = {'sharex': iax} if trn else  {'sharey': iax}
    pax = [fig.add_subplot(gsp[row, i+1], **share) for i in range(npl)]
    return [cax, iax] + pax


def _data_axes(fig: Figure, gsp: GridSpec, top: int, trn: bool) -> List[Axes]:
    """Create axes for an experiment/simulation.
    """
    pax = fig.add_subplot(gsp[top, :-1])
    share = {'sharey': pax} if trn else  {'sharex': pax}
    rax = fig.add_subplot(gsp[top+1, :-1], **share)
    sax = fig.add_subplot(gsp[top+2, :-1], **share)
    cpax = fig.add_subplot(gsp[top, -1])
    crax = fig.add_subplot(gsp[top+1, -1])
    csax = fig.add_subplot(gsp[top+2, -1])
    return [cpax, pax, crax, rax, csax, sax]


def _rep_every_other(txt: str, old: str, new: str) -> str:
    """Replace everyother occurence"""
    ind = -len(old)
    while old in txt[ind+len(old):]:
        ind = txt.find(old, ind+len(old))
        txt = txt[:ind] + new + txt[ind+len(old):]
        ind = txt.find(old, ind+len(new))
    return txt


# =============================================================================
# Helper class
# =============================================================================


class TransposeGridSpec:
    """A grid spec that transposes subscripts
    """
    gsp: GridSpec

    def __init__(self, gsp: GridSpec) -> None:
        self.gsp = gsp

    def __getitem__(self, ind: Inds) -> mpl.gridspec.SubplotSpec:
        return self.gsp[ind[::-1]]

    def __getattr__(self, attr: str) -> Any:
        return getattr(self.gsp, attr)


# =============================================================================
# Hint valiases
# =============================================================================
Inds = Tuple[Union[int, slice], ...]
Figure = mpl.figure.Figure
Axes = mpl.axes.Axes
Colorbar = mpl.colorbar.Colorbar
GridSpec = mpl.gridspec.GridSpec
Text = mpl.text.Text
