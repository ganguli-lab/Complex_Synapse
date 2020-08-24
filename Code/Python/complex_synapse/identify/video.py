"""Creating video frames for synapse fitting
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import sl_py_tools.arg_tricks as ag
import sl_py_tools.containers as cn
import sl_py_tools.matplotlib_tricks as mpt

from .. import synapse_base as _sb
from . import fit_synapse as _fs
from . import plast_seq as _ps

mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# =============================================================================
# Fitter video options class
# =============================================================================


class _Options:
    """Base class for video options"""

    def __contains__(self, key: str) -> bool:
        """Is key valid for self.update?"""
        return hasattr(self, 'set_' + key) or key in _sb.instance_attrs(self)

    def update(self, opt: cn.Dictable[str, Any] = (), **kwds) -> None:
        """Update options.
        """
        opt = dict(opt, **kwds)
        for key, val in opt.items():
            self.set(key, val)

    def set(self, key: str, val: Any) -> None:
        """Set an option.
        """
        if hasattr(self, 'set_' + key):
            getattr(self, 'set_' + key)(val)
        else:
            setattr(self, key, val)

    def __repr__(self) -> str:
        return type(self).__name__ + f"(**{self.__dict__})"


class VideoLabels(_Options):
    """Tiles, axes labels, etc. for FitterVideo

    Parameters
    ----------
    plast_type : List[str]
        Title, axis and colourbar labels for `fitter.data.plasticity_type`.
        By default: `["Plasticity type", "", "dep.", "pot."]`.
    readout : List[str]
        Title, axis and colourbar labels for `fitter.data.readout`. By default:
        `["Synaptic efficacy", "", "weak", "strong"]`.
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
    All of them are keyword only.
    """
    # Text for labels
    plast_type: List[str]
    readout: List[str]
    state: List[str]
    model: List[List[str]]
    plast: List[str]

    def __init__(self, **kwds) -> None:
        self.plast_type = ["Plasticity type", "", "dep.", "pot."]
        self.readout = ["Synaptic efficacy", "", "weak", "strong"]
        self.state = ["Occupation probability", "Time", "State", "True path"]
        self.model = ["Fit, $\\mathbf{M}$", "Truth, $\\mathbf{M}^*$"]
        self.plast = [["Potentiation", "Depression"],
                      ["Probability", "Initial state", "From state",
                       "To state"]]
        self.keys = {'transpose'} | _sb.instance_attrs(self).keys()
        self.update(kwds)

    def set_transpose(self, transpose: Optional[bool]) -> None:
        """Put the time-axis label on the lowest/leftmost of the data axes.

        Also adds/removes linebreaks from axes titles.
        Does nothing if `transpose` is `None`.
        """
        if transpose is None:
            return
        labels = (self.plast_type[1], self.readout[1], self.state[1])
        label = labels[list(map(bool, labels)).index(True)]
        labels = (label, "", "") if transpose else ("", "", label)
        self.plast_type[1], self.readout[1], self.state[1] = labels

        swap = (" ", "\n") if transpose else ("\n", " ")
        self.plast_type[0] = self.plast_type[0].replace(*swap)
        self.readout[0] = self.readout[0].replace(*swap)
        self.state[0] = self.state[0].replace(*swap)

    def model_labs(self, ind: _ps.Inds) -> Tuple[List[str], List[str]]:
        """Labels for a model's heatmaps"""
        return cn.listify(self.model[ind]) + self.plast[0], self.plast[1]


class VideoLayout(_Options):
    """Tiles, axes labels, etc. for FitterVideo

    Parameters
    ----------
    mrows : Sequence[int] (2,)
        which row for [estimated model, true model]
    mcols : Sequence[int] (4,)
        which column for [colourbar, initial, first,last+1 plasticity matrices]
    drows : Sequence[int] (3,)
        rows for (plast_type, readout, states),
    dcol : Sequence[int or slice] (2,)
        cols for (legend, plot)
    sizes : Dict[str, float]
        Grid size ratios for various axes.
        Widths (for model) - c: colourbar, i: `initial`, p: `plast`,
        Heights (for model) - md: model,
        Heights (for data) - pr: `plast_type` and `readout`, st: states,
        scale: Controls figure size: inches per grid-ratio unit.
        By default `None -> dict(c=1, i=2, p=12, pr=2, st=6, md=12, scale=0.3)`
    All of them are keyword only. Unknown keywords are passed to `sizes` if it
    is an existing key, or stored as an attribute if not.
    """
    # row/column assignments for models/data
    mrows: _ps.Inds
    mcols: _ps.Inds
    drows: _ps.Inds
    dcols: _ps.Inds
    # height/width ratios for rows/columns
    sizes: Dict[str, float]

    def __init__(self, **kwds) -> None:
        self.mrows = (-1, 0)
        self.mcols = (-1, 0, 1)
        self.drows = (0, 1, 2)
        self.dcols = (-1, slice(-1))
        # width: [c]bar=1, [i]nit=2, [p]last=12, {plast_seq=sum-cbar}
        # height: [m]o[d]el=12, [p]last_type/[r]eadout=2, [st]ates=6
        self.sizes = dict(c=1, i=2, p=12, pr=2, st=6, md=12, scale=0.3)
        self.update(kwds)

    def __contains__(self, key: str) -> bool:
        """Is key valid for self.update?"""
        return super().__contains__(key) or key in self.sizes

    def set(self, key: str, val: Any) -> None:
        """Set an option.

        If `sizes` is specified, it is updated with the new value, rather than
        being replaced like other variables.
        """
        if key == 'sizes':
            self.sizes.update(val)
        elif key in self.sizes:
            self.sizes[key] = val
        else:
            super().set(key, val)
        if key in {'mcols', 'mrows', 'dcols', 'drows'}:
            setattr(self, key, cn.tuplify(getattr(self, key)))

    def set_ground(self, ground: Optional[bool]) -> None:
        """Put the time-axis label on the lowest/leftmost of the data axes.

        Also adds/removes linebreaks from axes titles.
        Does nothing if `transpose` is `None`.
        """
        if ground is None:
            return
        change = int(ground) - self.drows[0]
        if change:
            self.drows = tuple(x + change for x in self.drows)

    def set_npl(self, npl: Optional[int]) -> None:
        """Put the time-axis label on the lowest/leftmost of the data axes.

        Also adds/removes linebreaks from axes titles.
        Does nothing if `transpose` is `None`.
        """
        if npl is None:
            return
        if len(self.mcols) - 2 != npl:
            self.mcols = self.mcols[:2] + tuple(self.mcols[2] + np.arange(npl))

    def set_fitter(self, fitter: Optional[_fs.SynapseFitter]) -> None:
        """Set `npl` and `ground` from a `SynapseFitter`"""
        self.set_npl(fitter.est.nplast)
        self.set_ground(isinstance(fitter, _fs.GroundedFitter))

    def grid_ratios(self) -> Tuple[List[int], List[int]]:
        """Ratios of row heights and column widths"""
        cols = ['p'] * len(self.mcols)
        rows = ['md'] * (4 + self.drows[0])
        for i, col in enumerate(('c', 'i')):
            cols[self.mcols[i]] = col
        for i, row in enumerate(('pr', 'pr', 'st')):
            rows[self.drows[i]] = row
        return [self.sizes[k] for k in rows], [self.sizes[k] for k in cols]



class VideoOptions(_Options):
    """Visual options for FitterVideo

    Parameters
    ----------
    txt : VideoLabels
        Title, axis, legend and colourbar labels.
    layout : VideoLayout
        Options for arrangement and sizes of video frames.
    transpose : bool
        Transpose the layout of the video? Stored in `txt_opt` and `im_opt`
        under the key `'trn'`. By default: `False`
    txt_opt : Dict[str, bool]
         Options for text. See `sl_py_tools.matplotlib_tricks.clean_axes`.
         By default: `{'box': False, 'tight': False}`.
    im_opt : Dict[str, str]
        Options for image plots. By default: `{'cmap': 'YlOrBr'}`.
    All of them are keyword only. Unknown keywords are passed to `txt_opt` if
    `sl_py_tools.matplotlib_tricks.clean_axes` takes them, `txt` or `layout`
    if it is an existing attribute, or `im_opt` if not.
    """
    # Text for labels
    txt: VideoLabels
    # Layout
    layout: VideoLayout
    # keyword options
    txt_opt: Dict[str, bool]
    im_opt: Dict[str, str]

    def __init__(self, **kwds) -> None:

        transpose = kwds.pop('transpose', False)
        self.txt = kwds.pop('txt', VideoLabels())
        self.layout = kwds.pop('sizes', VideoLayout())
        self.txt_opt = kwds.pop('txt_opt', {'box': False, 'tight': False})
        self.im_opt = kwds.pop('im_opt', {'cmap': 'YlOrBr'})

        for key, val in self.txt_opt.items():
            kwds.setdefault(key, val)
        self.txt_opt.update(mpt.clean_axes_keys(kwds))
        self.update(kwds)
        if transpose:
            self.txt_opt['trn'] = True
            self.im_opt['trn'] = True
        self.txt.set_transpose(transpose)

    def update(self, opt: cn.Dictable[str, Any] = (), **kwds) -> None:
        """Update options.

        Attributes are updated with the new value, rather than being replaced.
        """
        super().update(opt, **kwds)
        self.txt.set_transpose(self.transpose)

    def set(self, key: str, val: Any) -> None:
        """Set an option.
        """
        if key == 'transpose':
            self.transpose = val
            return
        myattrs = _sb.instance_attrs(self)
        for attr in myattrs.values():
            if key in attr:
                attr.update(**{key: val})
                return
        if key in myattrs:
            myattrs[key].update(val)
        else:
            self.im_opt[key] = val

    @property
    def transpose(self) -> bool:
        """Transpose the layout of the video?"""
        return self.txt_opt.get('trn', False)

    @transpose.setter
    def transpose(self, value: bool) -> None:
        """Transpose the layout of the video?"""
        self.txt_opt['trn'] = value
        self.im_opt['trn'] = value
        self.txt.set_transpose(value)

    @transpose.deleter
    def transpose(self) -> None:
        """Transpose the layout of the video?"""
        self.txt_opt.pop('trn', None)
        self.im_opt.pop('trn', None)
        self.txt.set_transpose(False)


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
    opt : VideoOptions, optional
        Choices for Axes labels, plot styles, etc
    norm : Optional[Normalize], optional keyword
        Data range for probability heatmaps, by default `Normalize(0, 1)`.

    Attributes
    ----------
    fig: Figure
        The figure containing video frames
    axh: Dict[str, List[Axes]]
        fit: Axes for estimated model [colourbar, initial, plasticity matrices]
        tr: for true model, [colourbar, initial, plasticity matrices] or []
        ps: for expt/sim heatmaps/plot, [plast_type, readouts, true path]
        pt: for plasticity types used [legend, heatmap]
        ro: for readouts types observed [legend, heatmap]
        st: for true/estimated path [legend, heatmap & plot]
        info: for display of the state of the fitter
    imh: Dict[str, List[Image or Line or Text]]
        fit: List[Image] - estimated model [initial, plasticity matrices]
        tr: List[Image] - true model [initial, plasticity matrices]
        ps: List[Image, Image, Line] - expt/sim [plast_type, readouts, path]
        st: Image - heatmap for estimated path
        info: Text - the display of the fitter's state
    """
    ind: _ps.Inds
    fname: str
    norm: mpl.colors.Normalize
    fig: Figure
    axh: Dict[str, AxList]
    imh: Dict[str, List[Disp]]
    opt: VideoOptions

    def __init__(self, fitter: _fs.SynapseFitter, ind: _ps.Inds, fname: str = "",
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
        self.opt.layout.set_fitter(fitter)
        figax = vid_fig(self.opt.transpose, self.opt.layout)
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
        ft_lab, tr_lab = self.opt.txt.model_labs(0), self.opt.txt.model_labs(1)
        lbo = self.opt.txt_opt
        mdo = {**self.opt.im_opt, 'zorder': 0, 'norm': self.norm}
        pso = {**self.opt.im_opt, 'zorder': 10, 'nplast': fitter.est.nplast,
               'nreadout': fitter.est.nreadout}

        self.imh['st'] = fitter.plot_occ(self.axh['st'][1], self.ind, **mdo)
        self.imh['ps'] = fitter.data[self.ind].plot(self.axh['ps'], **pso)
        label_data(self.axh['pt'], self.opt.txt.plast_type, leg=True, **lbo)
        label_data(self.axh['ro'], self.opt.txt.readout, leg=True, **lbo)
        label_sim(self.axh['st'], self.opt.txt.state, leg=self.ground, **lbo)

        self.imh['fit'] = fitter.est.plot(self.axh['fit'][1:], **mdo)
        label_model(self.axh['fit'], *ft_lab, cbar=True, **lbo)

        if self.ground:
            self.imh['tr'] = fitter.truth.plot(self.axh['tr'][1:], **mdo)
            label_model(self.axh['tr'], *tr_lab, cbar=False, **lbo)

        self.imh['info'] = write_info(format(fitter, 'tex'), self.axh['info'],
                                      trn=self.opt.transpose,
                                      size=lbo.get('tickfontsize', 10))
        # self.fig.canvas.draw_idle()
        plt.draw()
        # plt.show()

    def update_plots(self, fitter: _fs.SynapseFitter) -> None:
        """Update plots after iteration

        Parameters
        ----------
        obj : SynapseFitter
            Object performing fit whose state we update.
        """
        fitter.plot_occ(self.imh['st'], self.ind, trn=self.opt.transpose)
        fitter.est.plot(self.imh['fit'], trn=self.opt.transpose)
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


def vid_fig(trn: bool = False, layout: Optional[VideoLayout] = None,
            ) -> Tuple[Figure, AxList, AxList, AxList]:
    """Create a figure with axes for an experiment/simulation fit.

    Parameters
    ----------
    trn : bool, optional
        Transpose the layout of the video? By default `False`
    layout : VideoLayout, optional
        Options for arrangement and sizes of video frames.

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
    # width/height ratios, etc
    layout = ag.default(layout, VideoLayout())
    gnd = layout.drows[0]

    # figure, grid-spec
    fig, gsp = _make_fig(layout.grid_ratios(), trn, layout.sizes['scale'])

    # plast_seq
    ps_ax = _data_axes(fig, gsp, layout.drows, layout.dcols)
    # fit_model
    fit_ax = _model_axes(fig, gsp, layout.mrows[0], layout.mcols)
    # true_model
    tr_ax = _model_axes(fig, gsp, layout.mrows[1], layout.mcols) if gnd else []

    return fig, ps_ax, fit_ax, tr_ax


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


def label_data(axs: AxList, labels: List[str], leg: bool, **kwds) -> None:
    """Label axes and create colourbar for plasticity-type/readout sequence

    Parameters
    ----------
    axs : List[Axes] (2,)
        The Axes for the [colourbar/key, the image plot]
    labels : List[str] (P+1,) or (R+1)
        Text for [Axes title, each value on the colourbar/key]
    leg :  bool
        Create a legend?
    trn : bool, optional keyword
        Transpose the layout of the video? By default `False`
    Other keywords passed to `sl_py_tools.matplotlib_tricks.clean_axes`.
    """
    trn = kwds.pop('trn', False)
    if leg:
        cblabels = labels[2:]
        cbopt = {'orientation': 'horizontal'} if trn else {}
        cbopt.update(drawedges=True, aspect=5, ticks=np.arange(len(cblabels)))
        fig, imh = axs[1].get_figure(), axs[1].get_images()[0]
        cbh = fig.colorbar(imh, cax=axs[0], **cbopt)
        cbh.set_ticklabels(cblabels)
        if trn:
            axs[0].tick_params(labelrotation=45)
        mpt.clean_axes(axs[0], **kwds)
    else:
        axs[0].set(frame_on=False, xticks=[], yticks=[])

    axs[1].set_title(labels[0])
    if labels[1]:
        axs[1].set_ylabel(labels[1])
    if trn:
        axs[1].set_xticks([])
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

    axs[1].set_title(labels[0])
    axs[1].xaxis.set_ticks_position('bottom')
    axlabels = labels[2:0:-1] if trn else labels[1:3]
    axs[1].set_xlabel(axlabels[0])
    if axlabels[1]:
        axs[1].set_ylabel(axlabels[1])
    axs[1].set_aspect('auto')
    mpt.clean_axes(axs[1], **kwds)


def write_info(txt: str, axs: mpl.axes.Axes, **kwds) -> mpl.text.Text:
    """Write info in an axes

    Parameters
    ----------
    txt : str
        The text to write.
    axs : matplotlib.axes.Axes
        The `Axes` object on which we write.
    trn : bool, optional keyword
        Use a transposed layout? By default `False`.
    Other keywords passed to `matplotlib.text.Text`.

    Returns
    -------
    txh : matplotlib.text.Text
        `Text` object that was written.
    """
    trn = kwds.pop('trn', False)
    kwds.update(transform=axs.transAxes, ha='center')
    pos, kwds['va'] = ((.5, 3), 'top') if trn else ((.5, 0), 'bottom')
    return axs.text(*pos, txt, **kwds)

# =============================================================================
# Private helpers
# =============================================================================


def _make_fig(ratios: Tuple[List[int], List[int]], trn: bool, scale: float
              ) -> Tuple[Figure, GridSpec]:
    """Create a figure and grid spec"""
    keys = ['height_ratios', 'width_ratios']
    ratios = ratios[::-1] if trn else ratios

    fsiz = tuple(scale * sum(g) for g in reversed(ratios))
    gsiz, kwargs = tuple(map(len, ratios)), dict(zip(keys, ratios))

    fig = plt.figure(figsize=fsiz, frameon=True, constrained_layout=True)
    gsp = fig.add_gridspec(*gsiz, **kwargs)
    gsp = TransposeGridSpec(gsp) if trn else gsp
    return fig, gsp


def _model_axes(fig: Figure, gsp: GridSpec, row: int, cols: Sequence[int]
                ) -> AxList:
    """Create axes for a synapse model.

    cols : Sequence[int] (4,)
        which column for (colourbar, initial, first,last+1 plasticity matrices)
    """
    shared = 'x' if isinstance(gsp, TransposeGridSpec) else 'y'
    cax = fig.add_subplot(gsp[row, cols[0]])
    iax = fig.add_subplot(gsp[row, cols[1]])
    share = {'share' + shared: iax}
    pax = [fig.add_subplot(gsp[row, i], **share) for i in cols[2:]]
    return [cax, iax] + pax


def _data_axes(fig: Figure, gsp: GridSpec, rows: _ps.Inds, cols: _ps.Inds
               ) -> AxList:
    """Create axes for an experiment/simulation.
    """
    shared = 'y' if isinstance(gsp, TransposeGridSpec) else 'x'
    cpax = fig.add_subplot(gsp[rows[0], cols[0]])
    crax = fig.add_subplot(gsp[rows[1], cols[0]])
    csax = fig.add_subplot(gsp[rows[2], cols[0]])
    pax = fig.add_subplot(gsp[rows[0], cols[1]])
    share = {'share' + shared: pax}
    rax = fig.add_subplot(gsp[rows[1], cols[1]], **share)
    sax = fig.add_subplot(gsp[rows[2], cols[1]], **share)
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
