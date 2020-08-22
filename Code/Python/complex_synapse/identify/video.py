"""Creating video frames for synapse fitting
"""
from __future__ import annotations

from collections import ChainMap
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import sl_py_tools.arg_tricks as ag
import sl_py_tools.containers as cn
import sl_py_tools.matplotlib_tricks as mpt

from . import plast_seq as _ps
from . import fit_synapse as _fs

mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# =============================================================================
# Fitter video options class
# =============================================================================


class VideoOptions:
    """Visual options for FitterVideo

    Parameters
    ----------
    plast_type : List[str]
        Title and colourbar labels for `fitter.data.plasticity_type`.
        By default: `["Plasticity type", "dep.", "pot."]`.
    readout : List[str]
        Title and colourbar labels for `fitter.data.readout`. By default:
        `["Synaptic efficacy", "weak", "strong"]`.
    state : List[str]
        Title, axes, line labels for `fitter.data.state` and
        `fitter.model.plot_occ`. By default:
        `["Occupation probability", "Time", "State", "True path"]`.
    model : List[str]
        Overall labels for `fitter.model` and `fitter.truth`. By default:
        `["Fit, $\\mathbf{M}$", "Truth, $\\mathbf{M}^*$"]`.
    plast : List[str]
        Transition matrix titles, labels for [colourbar, initial state,
        from state, to state]. By default:
        `[["Potentiation", "Depression"],
          ["Probability", "Initial state", "From state", "To state"]]`.
    transpose : bool
        Transpose the layout of the video? Stored in `txt_opt` and `im_opt`
        under the key `'trn'`. By default: `False`
    txt_opt : Dict[str, bool]
         Options for text. See `sl_py_tools.matplotlib_tricks.clean_axes`.
         By default: `{'box': False, 'tight': False}`.
    im_opt : Dict[str, str]
        Options for image plots. By default: `{'cmap': 'YlOrBr'}`.
    All of them are keyword only. Unknown keywords are passed to `txt_opt` if
    `sl_py_tools.matplotlib_tricks.clean_axes` takes them, or `im_opt` if not.
    """
    # Text for labels
    plast_type: List[str]
    readout: List[str]
    state: List[str]
    model: List[List[str]]
    plast: List[str]
    # keyword options
    txt_opt: Dict[str, bool]
    im_opt: Dict[str, str]

    def __init__(self, **kwds) -> None:
        self.plast_type = kwds.pop('plast_type',
                                   ["Plasticity type", "dep.", "pot."])
        self.readout = kwds.pop('readout',
                                ["Synaptic efficacy", "weak", "strong"])
        self.state = kwds.pop('state', ["Occupation probability", "Time",
                                        "State", "True path"])
        self.model = kwds.pop('model',
                              ["Fit, $\\mathbf{M}$", "Truth, $\\mathbf{M}^*$"])
        self.plast = kwds.pop('plast', [["Potentiation", "Depression"],
                                        ["Probability", "Initial state",
                                         "From state", "To state"]])

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

    def model_labs(self, ind: _ps.Inds) -> Tuple[List[str], List[str]]:
        """Labels for a model's heatmaps"""
        return cn.listify(self.model[ind]) + self.plast[0], self.plast[1]

    def update(self, opt: cn.Dictable[str, Any], **kwds) -> None:
        """Update options.

        If `txt_opt` and `im_opt` are specified, they are updated with the new
        value, rather than being replaced like other variables.
        """
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
    norm : Optional[Normalize], optional keyword
        Data range for probability heatmaps, by default `Normalize(0, 1)`.
    """
    ind: _ps.Inds
    fname: str
    norm: mpl.colors.Normalize
    fig: Figure
    axh: Dict[str, AxList]
    imh: Dict[str, List[Disp]]
    opt: VideoOptions

    def __init__(self, npl: int, ground: bool, ind: _ps.Inds, fname: str = "",
                 opt: Optional[VideoOptions] = None, **kwds) -> None:
        """Video frame producing callback for SynapseFitter.

        Parameters
        ----------
        fitter : SynapseFitter
            Object that performs model fit
        ind : Tuple[Union[int, slice], ...]
            Indices for `fitter.data` for snippet to plot. Must result in
            `fitter.data[ind].nexpt = ()`.
        fname : str, optional
            Strings for filenames - produced by `fname.format(frams_number)`,
            by default `""` -> don't save files.
        opt : VideoOptions or None, optional
            Options for video display, by default `None` -> `VideoOptions()`.
        vmin, vmax : float, optional keyword
            Lower and Upper bounds for probability heatmaps, ignored if `norm`
            is given, by default `(0., 1.)`
        norm : mpl.colors.Normalize, optional keyword
            Range for probability heatmaps, by default `Normalize(vmin, vmax)`.
        Other keywords passed to `self.opt.update`
        """
        self.ind = ind
        self.fname = fname
        vmin, vmax = kwds.pop('vmin', 0.), kwds.pop('vmax', 1.)
        self.norm = kwds.pop('norm', mpl.colors.Normalize(vmin, vmax))
        self.opt = ag.default(opt, VideoOptions())
        self.opt.update(kwds)

        self.axh = {}
        self.imh = {}
        figax = vid_fig(npl, ground, self.opt.transpose)
        self.fig, psax, self.axh['fit'], self.axh['tr'] = figax
        self.axh['ps'] = psax[1::2]
        self.axh['pt'] = psax[:2]
        self.axh['ro'] = psax[2:4]
        self.axh['st'] = psax[4:]
        self.axh['info'] = self.axh['tr' if self.ground else 'st'][0]

    def __call__(self, fitter: _fs.SynapseFitter, pos: int) -> None:
        """Callback that displays fitter state as appropriate

        Parameters
        ----------
        obj : SynapseFitter
            Object performing fit whose state we display.
        pos : int
            At what stage of the fit are we?
                0: Before first iteration.
                1: During itertions.
                2: After completion.
        """
        if pos == 0:
            self.create_plots(fitter)
            _fs.print_callback(fitter, 0)
        elif pos == 1:
            if fitter.info['nit'] % fitter.opt.disp_step == 0:
                self.update_plots(fitter)
                self.savefig(fitter.info['nit'])
        elif pos == 2:
            _fs.print_callback(fitter, 2)

    @property
    def ground(self) -> bool:
        """Do we have ground truth?"""
        return bool(self.axh['tr'])

    def create_plots(self, fitter: _fs.SynapseFitter) -> None:
        """Create initial plots

        Parameters
        ----------
        obj : SynapseFitter
            Object performing fit whose state we display.
        """
        fit_lab, tru_lab = self.opt.model_labs(0), self.opt.model_labs(1)
        lbo = self.opt.txt_opt
        mdo = {**self.opt.im_opt, 'zorder': 0, 'norm': self.norm}
        pso = {**self.opt.im_opt, 'zorder': 10, 'nplast': fitter.est.nplast,
               'nreadout': fitter.est.nreadout}

        self.imh['st'] = fitter.plot_occ(self.axh['st'][1], self.ind, **mdo)
        self.imh['ps'] = fitter.data[self.ind].plot(self.axh['ps'], **pso)
        label_data(self.axh['pt'], self.opt.plast_type, **lbo)
        label_data(self.axh['ro'], self.opt.readout, **lbo)
        label_sim(self.axh['st'], self.opt.state, leg=self.ground, **lbo)

        self.imh['fit'] = fitter.est.plot(self.axh['fit'][1:], **mdo)
        label_model(self.axh['fit'], *fit_lab, cbar=True, **lbo)

        if self.ground:
            self.imh['tr'] = fitter.truth.plot(self.axh['tr'][1:], **mdo)
            label_model(self.axh['tr'], *tru_lab, cbar=False, **lbo)

        self.imh['info'] = write_info(format(fitter, 'tex'), self.axh['info'],
                                      trn=self.opt.transpose,
                                      size=lbo.get('tickfontsize', 10))
        # self.fig.canvas.draw_idle()
        plt.draw()
        plt.show()

    def update_plots(self, fitter: _fs.SynapseFitter) -> None:
        """Update plots after iteration

        Parameters
        ----------
        obj : SynapseFitter
            Object performing fit whose state we update.
        """
        trn = self.opt.transpose
        fitter.plot_occ(self.imh['st'], self.ind, trn=trn)
        fitter.est.plot(self.imh['fit'], trn=trn)
        self.imh['info'].set_text(format(fitter, 'tex'))
        plt.draw()
        # self.fig.canvas.draw_idle()

    def savefig(self, fileno: Union[None, int, str]) -> None:
        """Save current figure as a file

        Parameters
        ----------
        fileno : int, str or None
            Identifier for video frame to save. Passed to `self.fname.format`.
        """
        if self.fname:
            self.fig.savefig(self.fname.format(fileno))


# =============================================================================
# Helpers
# =============================================================================


def vid_fig(npl: int = 2, ground: bool = False, trn: bool = False
            ) -> Tuple[Figure, AxList, AxList, AxList]:
    """Create a figure with axes for an experiment/simulation fit.

    Parameters
    ----------
    npl : int, optional
        Number of types of plasticity, by default 2
    ground : bool, optional
        Do we have ground truth? By default `False`
    trn : bool, optional
        Transpose the layout of the video? By default `False`

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
    fig, gsp = _make_fig(npl, row=int(not ground), trn=trn)
    # true_model
    true_ax = _model_axes(fig, gsp, npl, row=0, trn=trn) if ground else []
    # plast_seq
    ps_ax = _data_axes(fig, gsp, top=int(ground), trn=trn)
    # fit_model
    fit_ax = _model_axes(fig, gsp, npl, row=-1, trn=trn)
    return fig, ps_ax, fit_ax, true_ax


def label_model(axs: AxList, titles: List[str], labels: List[str], cbar: bool,
                **kwds) -> None:
    """Label axes and create colourbar for model

    Parameters
    ----------
    axs : List[Axes] (P+2,)
        Axes for the [colourbar, initial state, each plasticity matrix ...]
    titles : List[str] (P+1,)
        Titles for [the full model, each plasticty matrix ...]
    labels : List[str] (4,)
        Labels for [colourbar, initial state, from state, to state]
    cbar :  bool
        Create a colourbar?
    trn : bool, optional keyword
        Transpose the layout of the video? By default `False`
    Other keywords passed to `sl_py_tools.matplotlib_tricks.clean_axes`.
    """
    trn = kwds.pop('trn', False)
    if cbar:
        fig, imh = axs[1].get_figure(), axs[1].get_images()[0]
        cbopt = {'orientation': 'horizontal'} if trn else {}
        fig.colorbar(imh, cax=axs[0], label=labels[0], **cbopt)
        mpt.clean_axes(axs[0], **kwds)
    else:
        axs[0].set(frame_on=False, xticks=[], yticks=[])

    axs[1].set_title(f"\\textit{{{titles[0]}}}", pad=20)
    if trn:
        axs[1].xaxis.set_ticks_position('bottom')
        axs[1].set(yticks=[], xlabel=labels[1])
    else:
        axs[1].set(xticks=[], ylabel=labels[1])
    mpt.clean_axes(axs[1], **kwds)

    for axh, lab in zip(axs[2:], titles[1:]):
        axh.set_ylabel(labels[2])
        axh.set_xlabel(labels[3])
        axh.set_title(lab)
        axh.xaxis.set_ticks_position('bottom')
        mpt.clean_axes(axh, **kwds)


def label_data(axs: AxList, labels: List[str], **kwds) -> None:
    """Label axes and create colourbar for plasticity-type/readout sequence

    Parameters
    ----------
    axs : List[Axes] (2,)
        The Axes for the [colourbar/key, the image plot]
    labels : List[str] (P+1,) or (R+1)
        Text for [Axes title, each value on the colourbar/key]
    trn : bool, optional keyword
        Transpose the layout of the video? By default `False`
    Other keywords passed to `sl_py_tools.matplotlib_tricks.clean_axes`.
    """
    trn = kwds.pop('trn', False)
    cbopt = {'orientation': 'horizontal'} if trn else {}
    fig, imh = axs[1].get_figure(), axs[1].get_images()[0]
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


def label_sim(axs: AxList, labels: List[str], leg: bool, **kwds) -> None:
    """Label axes and create legend for state sequence

    Parameters
    ----------
    axs : List[Axes] (2,)
        The axes for the [legend, the true path plot]
    labels : List[str] (3,)
        Text for the [title, time-axis, state-axis, true state line labels]
    leg :  bool
        Create a legend?
    trn : bool, optional keyword
        Transpose the layout of the video? By default `False`
    Other keywords passed to `sl_py_tools.matplotlib_tricks.clean_axes`.
    """
    trn = kwds.pop('trn', False)
    if leg:
        axs[0].legend(handles=axs[1].get_lines(), labels=labels[3:], loc=10)
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


def write_info(txt: str, axs: mpl.axes.Axes, **kwds) -> mpl.text.Text:
    """Write info next to states/truth axes"""
    trn = kwds.pop('trn', False)
    kwds.update(transform=axs.transAxes, ha='center')
    pos, kwds['va'] = ((.5, 3), 'top') if trn else ((.5, 0), 'bottom')
    return axs.text(*pos, txt, **kwds)

# =============================================================================
# Private helpers
# =============================================================================


def _make_fig(npl: int, row: int, trn: bool) -> Tuple[Figure, GridSpec]:
    """Create a figure and grid spec"""
    # width: cbar=1, init=2, plast=12, {plast_seq=sum-cbar}
    # height: model=12, plast_type=2, readout=2, states=6
    sizes = {'c': 1, 'i': 2, 'p': 12, 'pr': 2, 'st': 6, 'md': 12}
    scl = 0.3
    keys = ['height_ratios', 'width_ratios']

    ratios = ([sizes[k] for k in ['md', 'pr', 'pr', 'st', 'md'][row:]],
              [sizes[k] for k in ['i'] + ['p'] * npl + ['c']])
    ratios = ratios[::-1] if trn else ratios

    fsiz = tuple(scl * sum(g) for g in reversed(ratios))
    gsiz, kwargs = tuple(map(len, ratios)), dict(zip(keys, ratios))

    fig = plt.figure(figsize=fsiz, frameon=True, constrained_layout=True)
    gsp = fig.add_gridspec(*gsiz, **kwargs)
    gsp = TransposeGridSpec(gsp) if trn else gsp
    return fig, gsp


def _model_axes(fig: Figure, gsp: GridSpec, npl: int, row: int, trn: bool
                ) -> AxList:
    """Create axes for a synapse model.
    """
    cax = fig.add_subplot(gsp[row, -1])
    iax = fig.add_subplot(gsp[row, 0])
    share = {'sharex' if trn else 'sharey': iax}
    pax = [fig.add_subplot(gsp[row, i+1], **share) for i in range(npl)]
    return [cax, iax] + pax


def _data_axes(fig: Figure, gsp: GridSpec, top: int, trn: bool) -> AxList:
    """Create axes for an experiment/simulation.
    """
    pax = fig.add_subplot(gsp[top, :-1])
    share = {'sharey' if trn else 'sharex': pax}
    rax = fig.add_subplot(gsp[top+1, :-1], **share)
    sax = fig.add_subplot(gsp[top+2, :-1], **share)
    cpax = fig.add_subplot(gsp[top, -1])
    crax = fig.add_subplot(gsp[top+1, -1])
    csax = fig.add_subplot(gsp[top+2, -1])
    return [cpax, pax, crax, rax, csax, sax]


def _rep_every_other(txt: str, old: str, new: str) -> str:
    """Replace every-other occurence"""
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

    def __getitem__(self, ind: _ps.Inds) -> mpl.gridspec.SubplotSpec:
        return self.gsp[ind[::-1]]

    def __getattr__(self, attr: str) -> Any:
        return getattr(self.gsp, attr)


# =============================================================================
# Hint valiases
# =============================================================================
Figure = mpl.figure.Figure
GridSpec = mpl.gridspec.GridSpec
Disp = Union[mpl.lines.Line2D, _ps.Image, mpl.text.Text]
TxHandle = Union[mpl.axes.Axes, mpl.text.Text]
AxList = List[mpl.axes.Axes]
