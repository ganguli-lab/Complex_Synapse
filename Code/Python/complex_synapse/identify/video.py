# -*- coding: utf-8 -*-
"""Creating video frames for synapse fitting
"""
from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import tempfile
from pathlib import Path

import matplotlib as mpl
import matplotlib.animation as mpla
import matplotlib.pyplot as plt
import numpy as np

import sl_py_tools.arg_tricks as ag
import sl_py_tools.containers as cn
import sl_py_tools.graph_plots as gp
import sl_py_tools.matplotlib_tricks as mpt
import sl_py_tools.options_classes as op
from sl_py_tools.graph_tricks import MultiDiGraph

import complex_synapse.identify.fit_synapse as fs
import complex_synapse.identify.plast_seq as ps

mpt.rc_colours()
mpt.rc_fonts('sans-serif')
_log = logging.getLogger(__name__)
# =============================================================================
# Fitter video options classes
# =============================================================================


# pylint: disable=too-many-ancestors
class VideoLabels(op.Options):
    """Tiles, axes labels, etc. for `FitterPlots`.

    The individual options can be accessed as object instance attributes
    (e.g. `obj.name`) or as dictionary items (e.g. `obj['name']`) for both
    getting and setting.

    Parameters
    ----------
    plast_type : List[str]
        Title, axis and colourbar labels for `fitter.data.plasticity_type`.
        By default: `["Plasticity type", "", "dep.", "pot."]`.
    readout : List[str]
        Title, axis and colourbar labels for `fitter.data.readout`. By default:
        `["Synaptic efficacy", "", "weak", "strong"]`.
    state : List[str]
        Title, axes and line labels for `fitter.data.state` and
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
        Transpose the layout of the video, swapping rows and columns?

    All parameters are optional keywords. Any dictionary passed as positional
    parameters will be popped for the relevant items. Keyword parameters must
    be valid keys, otherwise a `KeyError` is raised.
    """
    prop_attributes: op.Attrs = ('transpose',)
    key_last: op.Attrs = ('transpose',)
    # Text for labels
    plast_type: List[str]
    readout: List[str]
    state: List[str]
    model: List[List[str]]
    plast: List[str]

    def __init__(self, *args, **kwds) -> None:
        self.plast_type = ["Plasticity type", "", "dep.", "pot."]
        self.readout = ["Synapse strength", "", "weak", "strong"]
        self.state = ["Occupation probability", "Time", "State", "True path"]
        self.model = ["Fit, $\\mathbf{M}$", "Truth, $\\mathbf{M}^*$"]
        self.plast = [["Potentiation", "Depression"],
                      ["Probability", "Initial state", "From state",
                       "To state"]]
        # args = op.sort_dicts(args, ('transpose',), -1)
        # kwds = op.sort_dict(kwds, ('transpose',), -1)
        super().__init__(*args, **kwds)

    def model_labs(self, ind: ps.Inds) -> Tuple[List[str], List[str]]:
        """Labels for a model's heatmaps."""
        return cn.listify(self.model[ind]) + self.plast[0], self.plast[1]

    @property
    def transpose(self) -> bool:
        """Do we swap rows and columns?"""
        return bool(self.plast_type[1])

    @transpose.setter
    def transpose(self, value: Optional[bool]) -> None:
        """Put the time-axis label on the lowest/leftmost of the data axes.

        Also adds/removes linebreaks from axes titles.
        Does nothing if `transpose` is `None`.

        It is better to set `VideoOptions.transpose` in the parent
        object instead of setting this property directly.
        """
        if value is None or value == self.transpose:
            return
        # label = self.plast_type[1] + self.readout[1] + self.state[1]
        # labels = (label, "", "") if transpose else ("", "", label)
        # self.plast_type[1], self.readout[1], self.state[1] = labels
        label = self.state[1] if value else ""
        self.plast_type[1], self.readout[1] = label, label

        swap = (" ", "\n") if value else ("\n", " ")
        self.plast_type[0] = self.plast_type[0].replace(*swap)
        self.readout[0] = self.readout[0].replace(*swap)
        self.state[0] = self.state[0].replace(*swap)


# pylint: disable=too-many-ancestors
class VideoLayout(op.Options):
    """Axes sizes, positions etc. for `FitterPlots`.

    The individual options can be accessed as object instance attributes
    (e.g. `obj.name`) or as dictionary items (e.g. `obj['name']`) for both
    getting and setting.

    Parameters
    ----------
    mrows : Sequence[int] (2,)
        which row for [estimated model, true model].
    mcols : Sequence[int] (4,)
        which column for [colourbar, initial, first,last+1 plasticity matrices]
    drows : Sequence[int] (3,)
        rows for (plast_type, readout, states).
    dcol : Sequence[int or slice] (2,)
        cols for (legend, plot).
    sizes : Dict[str, float]
        Grid size ratios for various axes.
        Widths (for model) - c: colourbar, i: `initial`, p: `plast`,
        Heights (for model) - md: model,
        Heights (for data) - pr: `plast_type` and `readout`, st: `states`,
        scale: Controls figure size: inches per grid-ratio unit.
        By default `None -> dict(c=1, i=2, p=12, pr=2, st=6, md=12, scale=0.3)`
    transpose : bool
        Transpose the layout of the video, swapping rows and columns?
    ground : bool
        Do we have axes for a ground truth model?
    npl : int
        How many plasticity types are there?
    verbosity : int
        Display verbose information?
    fitter : SynapseFitter
        Get `npl`, `ground` and `verbose` from `fitter.est`, `truth` and `opt`.
    When constructing a `FitterVideo`, `set_fitter` is called. Therefore,
    there is usually no need to provide `npl`, `ground` or `verbose`.

    All parameters are optional keywords. Any dictionary passed as positional
    parameters will be popped for the relevant items. Keyword parameters must
    be valid keys, otherwise a `KeyError` is raised.

    Notes
    -----
    Attributes can also be referenced and modified by subscripting with the
    attribute name. If the name is not found, it will search `sizes`.
    """
    map_attributes: op.Attrs = ('sizes',)
    prop_attributes: op.Attrs = ('transpose', 'ground', 'npl', 'verbosity')
    key_last: op.Attrs = ('transpose', 'ground', 'npl', 'verbosity', 'fitter')
    # row/column assignments for models/data
    mrows: ps.Inds
    mcols: ps.Inds
    drows: ps.Inds
    dcols: ps.Inds
    # height/width ratios for rows/columns
    sizes: Dict[str, float]
    # swap rows and columns?
    _transpose: bool
    # display verbose information?
    _verbosity: int

    def __init__(self, *args, **kwds) -> None:
        self.mrows = (-1, 0)
        self.mcols = (-1, -2, 0, 1)
        self.drows = (0, 1, 2)
        self.dcols = (-1, slice(-1))
        # width: [c]olourbar, [g]raph, [i]nitial, [p]last, {data=sum-cbar}
        # height: [m]o[d]el, [p]last_type/[r]eadout, [st]ates
        self.sizes = dict(c=1, g=8, i=2, p=12, pr=2, st=6, md=12, scale=0.3)
        self._transpose = False
        self._verbosity = 1
        # order = ('transpose', 'ground', 'npl', 'verbosity', 'fitter')
        # args = op.sort_dicts(args, order, -1)
        # kwds = op.sort_dict(kwds, order, -1)
        super().__init__(*args, **kwds)

    def __setitem__(self, key: str, val: Any) -> None:
        """Set an option.

        If `sizes` is specified, it is updated with the new value, rather than
        being replaced like other variables.
        """
        super().__setitem__(key, val)
        if key in {'mcols', 'mrows', 'dcols', 'drows'}:
            setattr(self, key, cn.tuplify(getattr(self, key)))

    def grid_ratios(self) -> Tuple[List[int], List[int]]:
        """Ratios of row heights and column widths.

        Returns
        -------
        height_ratios : List[int]
            Ratios of row heights
        width_ratios : List[int]
            Ratios of column widths
        """
        cols = ['p'] * len(self.mcols)
        rows = ['md'] * (4 + self.drows[0])
        for i, col in enumerate(('c', 'g', 'i')):
            cols[self.mcols[i]] = col
        for i, row in enumerate(('pr', 'pr', 'st')):
            rows[self.drows[i]] = row
        ratios = [self.sizes[k] for k in rows], [self.sizes[k] for k in cols]
        return ratios[::-1] if self._transpose else ratios

    def gspec_opts(self) -> Tuple[Tuple[float, float], Tuple[int, int],
                                  Dict[str, List[int]]]:
        """Options for creating a `Figure` and `GridSpec`."""
        keys = ['height_ratios', 'width_ratios']
        ratios = self.grid_ratios()

        fsiz = tuple(self.sizes['scale'] * sum(g) for g in reversed(ratios))
        gsiz, kwargs = tuple(map(len, ratios)), dict(zip(keys, ratios))
        return fsiz, gsiz, kwargs

    def set_fitter(self, fitter: Optional[fs.SynapseFitter]) -> None:
        """Set `npl`, `ground` and `verbose` from a `SynapseFitter`.

        Does nothing if `fitter` is `None`.
        """
        if fitter is None:
            return
        self.npl = fitter.est.nplast
        self.ground = isinstance(fitter, fs.GroundedFitter)
        self.verbosity = fitter.opt.disp_each

    @property
    def transpose(self) -> bool:
        """Do we swap rows and columns?"""
        return self._transpose

    @transpose.setter
    def transpose(self, transpose: Optional[bool]) -> None:
        """Choose whether to swap rows and columns.

        Does nothing if `transpose` is `None`.

        It is better to set `VideoOptions.transpose` in the parent
        object instead of setting this property directly.
        """
        if transpose is None:
            return
        self._transpose = transpose

    @property
    def ground(self) -> bool:
        """Do we have axes for a ground truth model?"""
        return bool(self.drows[0])

    @ground.setter
    def ground(self, ground: Optional[bool]) -> None:
        """Include axes for ground truth model?

        Does nothing if `ground` is `None`.
        """
        if ground is None:
            return
        change = int(ground) - self.drows[0]
        if change:
            self.drows = tuple(x + change for x in self.drows)

    @property
    def npl(self) -> int:
        """How many plasticity types are there?"""
        return len(self.mcols) - 3

    @npl.setter
    def npl(self, npl: Optional[int]) -> None:
        """Set the number of plasticity types.

        Does nothing if `npl` is `None`.
        """
        if npl is None:
            return
        if len(self.mcols) - 3 != npl:
            self.mcols = self.mcols[:3] + tuple(self.mcols[3] + np.arange(npl))

    @property
    def verbosity(self) -> int:
        """Do we display verbose information?
        """
        return self._verbosity

    @verbosity.setter
    def verbosity(self, verbose: Optional[int]) -> None:
        """Choose whether to display verbose information.

        Does nothing if `verbose` is `None`.
        """
        if verbose is None:
            return
        self._verbosity = verbose


# pylint: disable=too-many-ancestors
class VideoOptions(op.MasterOptions, fallback='im_opt'):
    """Visual options for `FitterPlots`.

    The individual options can be accessed as object instance attributes
    (e.g. `obj.name`) or as dictionary items (e.g. `obj['name']`) for both
    getting and setting.

    Parameters
    ----------
    txt : VideoLabels
        Title, axis, legend and colourbar labels.
    layout : VideoLayout
        Options for arrangement and sizes of video frames.
    transpose : bool
        Transpose the layout of the video? Stored in `txt_opt` and `im_opt`
        under the key `'trn'`. By default: `False`
    ax_opt : Dict[str, Any]
         Options for `Axes`. See `sl_py_tools.matplotlib_tricks.clean_axes`.
         By default: `{'box': False, 'tight': False}`.
    im_opt : ImageOptions
        Options for image plots. By default: `ImageOptions()`.
    ln_opt : Dict[str, Any]
        Options for line plots. By default: `{}`.
    an_opt : AnimationOptions
        Options for animating the video.
    gr_opt : GraphOptions
        Options for drawing the graph of the models.

    All parameters are optional keywords. Any dictionary passed as positional
    parameters will be popped for the relevant items.

    Notes
    -----
    All parameters are keyword only. Unknown keywords are passed to `ax_opt`
    if `sl_py_tools.matplotlib_tricks.clean_axes` takes them, `txt`, `layout`
    or `ln_opt` if it is an existing attribute/key, or `im_opt` if not.

    Attributes can also be referenced and modified by subscripting with the
    attribute name. If the name is not found, it will search `txt`, `layout`,
    `ax_opt`, `ln_opt` and `im_opt`. Setting a new key adds it to `im_opt`.
    To add a new item to a `VideoOptions` instance, set it as an attribute.
    """
    map_attributes: op.Attrs = ('txt', 'layout', 'image', 'line', 'axes',
                                'animate', 'graph')
    prop_attributes: op.Attrs = ('transpose',)
    # Text for labels
    txt: VideoLabels
    # Layout
    layout: VideoLayout
    # Heatmaps
    image: mpt.ImageOptions
    # Line plots
    line: Dict[str, Any]
    # Fonts for labels
    axes: mpt.AxesOptions
    # Animation
    animate: mpt.AnimationOptions
    # Graph
    graph: gp.GraphOptions

    def __init__(self, *args, **kwds) -> None:
        self.txt = VideoLabels()
        self.layout = VideoLayout()
        self.line = {}
        self.image = mpt.ImageOptions(trn=False)
        self.axes = mpt.AxesOptions(box=False, tight=False)
        # self.ax_opt.trn = False
        self.animate = mpt.AnimationOptions()
        self.graph = gp.GraphOptions({'nodes.mult': 600, 'edges.mult': 2,
                                       'edges.mut_scale': 3})
        self.graph.set_layout(gp.linear_layout, sep=(0.0, 1.0))

        super().__init__(*args, **kwds)

    @property
    def transpose(self) -> bool:
        """Transpose the layout of the video?"""
        return self.axes.get('trn', False)

    @transpose.setter
    def transpose(self, value: bool) -> None:
        """Transpose the layout of the video?"""
        if value is None:
            return
        # self.ax_opt['trn'] = transpose
        self.image['trn'] = value
        self.txt.transpose = value
        self.layout.transpose = value
        sep = (1.0, 0.0) if value else (0.0, 1.0)
        self.graph.set_layout(gp.linear_layout, sep=sep)
# pylint: enable=too-many-ancestors


# =============================================================================
# Animation
# =============================================================================


def animate(fitter: fs.SynapseFitter, **kwargs) -> mpla.FuncAnimation:
    """Animate a fitter video.

    Parameters
    ----------
    fitter : SynapseFitter
        The synapse fitter we will animate. An instance of `FitterReplay` or
        one of its subclasses is recommended. Its `callback` must be an
        instance of `FitterPlots`.
    Other keywords
        Passed to `FuncAnimation`.

    Returns
    -------
    ani : FuncAnimation
        The animation.

    To view the animation call `plt.show()`.
    To save a video, call `ani.save(...)`.
    See <https://matplotlib.org/api/animation_api.html>.
    """
    if not isinstance(fitter.callback, FitterPlots):
        raise TypeError
    _log.debug("Calling build before creating animation")
    fitter.callback.build(fitter)
    kwargs.setdefault('init_func', fitter.init)
    kwargs.setdefault('frames', fitter.opt.max_it)
    opt = fitter.callback.opt.an_opt.copy()
    opt.update(kwargs)
    _log.debug("Creating animation")
    return mpla.FuncAnimation(fitter.callback.fig, fitter.step, **opt)


# =============================================================================
# Fitter video class
# =============================================================================


class FitterPlots:
    """Class to produce video frames showing fitter in action.

    Parameters
    ----------
    inds : Inds
        Indices for `fitter.data` for snippet to plot. Must result in
        `fitter.data[ind].nexpt = ()`.
    opt : VideoOptions, optional
        Choices for Axes labels, plot styles, etc.
    norm : Optional[Normalize], optional keyword
        Data range for probability heatmaps, by default `Normalize(0, 1)`.
    Other keywords passed to `self.opt.update`.

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
    imh: Dict[str, List[Image|Line|Text]]
        fit: List[Image] - estimated model [initial, plasticity matrices]
        tr: List[Image] - true model [initial, plasticity matrices]
        ps: List[Plot] - experiment/simulation data [plast_type, readout, path]
        st: List[Image] - heatmap for estimated path
        info: List[Text] - the display of the fitter's state
    """
    ind: ps.Inds
    norm: mpl.colors.Normalize
    fig: Optional[Figure]
    axh: Dict[str, AxList]
    imh: Dict[str, List[Disp]]
    grf: Dict[str, MultiDiGraph]
    opt: VideoOptions

    def __init__(self, ind: ps.Inds, opt: Optional[VideoOptions] = None,
                 **kwds) -> None:
        """Video frame producing callback for `SynapseFitter`.

        Parameters
        ----------
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
        self.opt = ag.default(opt, VideoOptions())
        self.opt.update(kwds)
        self.norm = self.opt.image.pop('norm')
        self.fig = None
        self.axh = {}
        self.imh = {}
        self.grf = {}

    def __call__(self, fitter: fs.SynapseFitter, pos: int) -> List[Disp]:
        """Callback that displays fitter state as appropriate

        Parameters
        ----------
        obj : SynapseFitter
            Object performing fit whose state we display.
        pos : int
            At what stage of the fit are we?
                0: Before first iteration.
                1: During iterations.
                2: After completion.
        """
        if pos == 0:
            self.build(fitter)
            fs.print_callback(fitter, 0)
        elif pos == 1:
            self.update_plots(fitter)
        elif pos == 2:
            if fitter.info['nit'] % fitter.opt.disp_step:
                self.update_plots(fitter)
            fs.print_callback(fitter, 2)
            _log.debug("Final callback")

        # self.fig.canvas.draw_idle()
        # plt.draw()
        # plt.show()
        if self.opt.animate.blit:
            return self.updated()
        return None

    @property
    def ground(self) -> bool:
        """Do we have ground truth?"""
        return bool(self.axh['tr'])

    def _make_fig(self) -> None:
        """Create the figure and axes for the video frames"""
        _log.debug("Create figure and axes")
        figax = vid_fig(self.opt.layout)
        self.fig, psax, self.axh['fit'], self.axh['tr'] = figax
        # All PlasticitySequence main plots
        self.axh['ps'] = psax[1::2]
        # Individual PlasticitySequence plots & their legends
        self.axh['pt'] = psax[:2]
        self.axh['ro'] = psax[2:4]
        self.axh['st'] = psax[4:]
        # Which axes to use to display `fitter.info`
        self.axh['info'] = self.axh['tr' if self.ground else 'st'][:1]

    def _create_plots(self, fitter: fs.SynapseFitter) -> None:
        """Create initial plots

        Should only be called after `make_fig`.

        Parameters
        ----------
        fitter : SynapseFitter
            Object performing fit whose state we display.
        """
        _log.debug("Create plots")
        mdo = {**self.opt.image, 'zorder': 0, 'norm': self.norm}
        pso = {**self.opt.image, 'zorder': 10, 'nplast': fitter.est.nplast,
               'nreadout': fitter.est.nreadout, 'line_opts': self.opt.line}
        txo = {'size': self.opt.axes.get('tickfontsize', 10),
               'clip_on': False}
        verbose = self.opt.layout.verbosity

        self.imh['st'] = [fitter.plot_occ(self.axh['st'][1], self.ind, **mdo)]
        self.imh['ps'] = fitter.data[self.ind].plot(self.axh['ps'], **pso)

        mdo['gopts'] = self.opt.graph
        (self.imh['fit'],
         self.grf['fit']) = fitter.est.plot(self.axh['fit'][1:], **mdo)

        if self.ground:
            (self.imh['tr'],
             self.grf['tr']) = fitter.truth.plot(self.axh['tr'][1:], **mdo)

        if verbose:
            self.imh['info'] = [write_info(format(fitter, f'tex0,{verbose}'),
                                           self.axh['info'][0], **txo)]

    def _format_axes(self, fitter: fs.SynapseFitter) -> None:
        """Format axes labels, legends and colourbars

        Should only be called after `create_plots`.
        """
        _log.debug("Format axes, except info")
        ft_lab, tr_lab = self.opt.txt.model_labs(0), self.opt.txt.model_labs(1)
        lbo = {'trn': self.opt.transpose, **self.opt.axes}
        nst = fitter.est.nstate

        label_data(self.axh['pt'], self.opt.txt.plast_type, leg=True, **lbo)
        label_data(self.axh['ro'], self.opt.txt.readout, leg=True, **lbo)
        label_sim(self.axh['st'], self.opt.txt.state, leg=self.ground, **lbo)

        lbo['rad'] = self.opt.graph.rad
        label_model(self.axh['fit'], *ft_lab, cbar=True, nst=nst, **lbo)
        self.axh['info'][0].set_clip_on(False)

        if self.ground:
            label_model(self.axh['tr'], *tr_lab, cbar=False, nst=nst, **lbo)

    def _info_axes(self, **kwargs) -> None:
        """Fix the size for the Text's Axes

        Should be called after first `draw`.
        """
        _log.debug('resizing info axes')
        if self.opt.layout.verbosity:
            # make sure we have a renderer
            frmt = kwargs.pop('format', 'pdf')
            with tempfile.TemporaryDirectory() as tmpdirname:
                self.savefig(Path(tmpdirname, 'blah'), format=frmt)
            # get bbox of text and transform to figure coordinates
            txt_bbox = self.imh['info'][0].get_window_extent().frozen()
            transf = self.fig.transFigure.inverted()
            # patch = self.imh['info'][0].get_bbox_patch()
            # txt_bbox = patch.get_bbox()
            # transf = patch.get_transform() - self.fig.transFigure
            # set bbox of axes in figure coordinates
            ax_bbox = transf.transform_bbox(txt_bbox)
            self.axh['info'][0].set_position(ax_bbox)
            self.imh['info'][0].set_in_layout(False)
            self.imh['info'][0].set_position((0, 0.5))
            self.imh['info'][0].set_ha('left')

    def build(self, fitter: fs.SynapseFitter, **kwds) -> None:
        """Create the figure and axes with the first plot on them

        Parameters
        ----------
        fitter : SynapseFitter
            The object that performs the fit, whose state we display.
        """
        if self.fig is None:
            self.opt.layout.set_fitter(fitter)
            self._make_fig()
            self._create_plots(fitter)
            self._format_axes(fitter)
            # self.fig.canvas.draw_idle()
            # plt.draw()
            self._info_axes(**kwds)
        # self.fig.set_constrained_layout(False)

    def update_plots(self, fitter: fs.SynapseFitter) -> None:
        """Update plots after iteration

        Parameters
        ----------
        obj : SynapseFitter
            Object performing fit whose state we update.
        """
        _log.debug("update plots")
        # if self.imh['info'][0].get_in_layout():
        #     self._info_axes()
        trn, verbose = self.opt.transpose, self.opt.layout.verbosity
        fitter.plot_occ(self.imh['st'][0], self.ind, trn=trn)
        fitter.est.plot(self.imh['fit'], self.grf['fit'], trn=trn)
        if verbose:
            self.imh['info'][0].set_text(format(fitter, f'tex1,{verbose}'))

    def updated(self) -> List[Disp]:
        """The artists that are updated at each step
        """
        out = self.imh['st'] + self.imh['fit'][1:] + self.imh['ps'][2:]
        out += self.imh['fit'][0].collection + self.imh.get('info', [])
        return out

    def savefig(self, filename: str, *args, **kwds) -> None:
        """Save current figure as a file

        Parameters
        ----------
        filename : int, str or None
            Name of file to save.
        """
        self.fig.savefig(filename, *args, **kwds)


# =============================================================================
# Helpers
# =============================================================================


def vid_fig(layout: Optional[VideoLayout] = None,
            ) -> Tuple[Figure, AxList, AxList, AxList]:
    """Create a figure with axes for an experiment/simulation fit.

    Parameters
    ----------
    trn : bool, optional
        Transpose the layout of the video? By default `False`.
    layout : VideoLayout, optional
        Options for arrangement and sizes of video frames.

    Returns
    -------
    fig : Figure
        The figure object for video frames.
    ps_ax : List[Axes]
        The axes objects for the `PlasticitySequence` data.
    fit_ax : List[Axes]
        The axes objects for the `SynapseIdModel` fit.
    true_ax : List[Axes]
        The axes objects for the ground-truth `SynapseIdModel`.
    """
    # width/height ratios, etc
    layout = ag.default(layout, VideoLayout())
    gnd = layout.ground
    fsiz, gsiz, kwargs = layout.gspec_opts()

    # figure, grid-spec
    fig = plt.figure(figsize=fsiz, frameon=True, constrained_layout=True)
    gsp = fig.add_gridspec(*gsiz, **kwargs)
    gsp = TransposeGridSpec(gsp) if layout.transpose else gsp

    # plast_seq
    ps_ax = _data_axes(fig, gsp, layout.drows, layout.dcols)
    # fit_model
    fit_ax = _model_axes(fig, gsp, layout.mrows[0], layout.mcols)
    # true_model
    tr_ax = _model_axes(fig, gsp, layout.mrows[1], layout.mcols) if gnd else []

    return fig, ps_ax, fit_ax, tr_ax


def label_model(axs: AxList, titles: List[str], labels: List[str], cbar: bool,
                nst: int, **kwds) -> None:
    """Label axes and create colourbar for model.

    Parameters
    ----------
    axs : List[Axes] (P+2,)
        Axes for the [colourbar, initial state, each plasticity matrix ...].
    titles : List[str] (P+1,)
        Titles for [the full model, each plasticty matrix ...].
    labels : List[str] (4,)
        Labels for [colourbar, initial state, from state, to state].
    cbar :  bool
        Create a colourbar?
    trn : bool, optional keyword
        Transpose the layout of the video? By default `False`.
    Other keywords passed to `sl_py_tools.matplotlib_tricks.clean_axes`.
    """
    trn = kwds.pop('trn', False)
    rad = np.max(np.abs((kwds.pop('rad', [-0.7, 0.35])))) / 1.95
    if cbar:
        fig, imh = axs[2].get_figure(), axs[2].get_images()[0]
        cbopt = {'orientation': 'horizontal'} if trn else {}
        fig.colorbar(imh, cax=axs[0], label=labels[0], **cbopt)
        mpt.clean_axes(axs[0], **kwds)
    else:
        axs[0].set(frame_on=False, xticks=[], yticks=[])

    axs[1].set(frame_on=False, xticks=[], yticks=[], aspect=1)
    para, perp = (-0.5, nst - 0.5), (-rad * nst, rad * nst)
    if trn:
        axs[1].set(xlim=para, ylim=perp)
    else:
        axs[1].set(xlim=perp[::-1], ylim=para[::-1])

    axs[2].set_title(f"\\textit{{{titles[0]}}}", pad=20)
    if trn:
        axs[2].xaxis.set_ticks_position('bottom')
        axs[2].set(yticks=[], xlabel=labels[1])
    else:
        axs[2].set(xticks=[], ylabel=labels[1])
    mpt.clean_axes(axs[2], **kwds)

    for axh, lab in zip(axs[3:], titles[1:]):
        axh.set_ylabel(labels[2])
        axh.set_xlabel(labels[3])
        axh.set_title(lab)
        axh.xaxis.set_ticks_position('bottom')
        mpt.clean_axes(axh, **kwds)


def label_data(axs: AxList, labels: List[str], leg: bool, **kwds) -> None:
    """Label axes and create colourbar for plasticity-type/readout sequence.

    Parameters
    ----------
    axs : List[Axes] (2,)
        The Axes for the [colourbar/key, the image plot].
    labels : List[str] (P+1,) or (R+1)
        Text for [Axes title, each value on the colourbar/key].
    leg :  bool
        Create a legend?
    trn : bool, optional keyword
        Transpose the layout of the video? By default `False`.
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
    if trn:
        axs[1].set_ylabel(labels[1])
        axs[1].set_xticks([])
    else:
        axs[1].set_xlabel(labels[1])
        axs[1].set_yticks([])
        axs[1].xaxis.set_ticks_position('bottom')
    axs[1].set_aspect('auto')
    mpt.clean_axes(axs[1], **kwds)


def label_sim(axs: AxList, labels: List[str], leg: bool, **kwds) -> None:
    """Label axes and create legend for state sequence.

    Parameters
    ----------
    axs : List[Axes] (2,)
        The axes for the [legend, the true path plot].
    labels : List[str] (3,)
        Text for the [title, time-axis, state-axis, true state line labels].
    leg :  bool
        Create a legend?
    trn : bool, optional keyword
        Transpose the layout of the video? By default `False`.
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
    axs[1].set_ylabel(axlabels[1])
    axs[1].set_aspect('auto')
    mpt.clean_axes(axs[1], **kwds)


def write_info(txt: str, axs: mpl.axes.Axes, **kwds) -> mpl.text.Text:
    """Write info in an axes.

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
    ax_colour = axs.get_facecolor()
    kwds['bbox'] = {'edgecolor': ax_colour, 'facecolor': ax_colour}
    kwds.update(transform=axs.transAxes, ha='center', va='center')
    return axs.text(0.5, 0.5, txt, **kwds)

# =============================================================================
# Private helpers
# =============================================================================


def _model_axes(fig: Figure, gsp: GridSpec, row: int, cols: Sequence[int]
                ) -> AxList:
    """Create axes for a synapse model.

    cols : Sequence[int] (3+P,)
        which column for (colourbar, graph, initial, plasticity matrices).
    """
    mdl = 'fit' if row else 'truth'
    shared = 'x' if isinstance(gsp, TransposeGridSpec) else 'y'
    cax = fig.add_subplot(gsp[row, cols[0]], label=f"{mdl}-cbar")
    gax = fig.add_subplot(gsp[row, cols[1]], label=f"{mdl}-graph")
    iax = fig.add_subplot(gsp[row, cols[2]], label=f"{mdl}-init")
    share = {'share' + shared: iax}
    pax = [fig.add_subplot(gsp[row, i], label=f"{mdl}-plast[{i-1}]", **share)
           for i in cols[3:]]
    return [cax, gax, iax] + pax


def _data_axes(fig: Figure, gsp: GridSpec, rows: ps.Inds, cols: ps.Inds
               ) -> AxList:
    """Create axes for an experiment/simulation.

    rows : Sequence[int] (3,)
        Which row for (plasticity types, readouts, state path).
    cols : Sequence[int] (2,)
        Which column for (legend, heatmap).
    """
    shared = 'y' if isinstance(gsp, TransposeGridSpec) else 'x'
    cpax = fig.add_subplot(gsp[rows[0], cols[0]], label="plast_type legend")
    crax = fig.add_subplot(gsp[rows[1], cols[0]], label="readout legend")
    csax = fig.add_subplot(gsp[rows[2], cols[0]], label="state legend")
    pax = fig.add_subplot(gsp[rows[0], cols[1]], label="plast_type")
    share = {'share' + shared: pax}
    rax = fig.add_subplot(gsp[rows[1], cols[1]], label="readout", **share)
    sax = fig.add_subplot(gsp[rows[2], cols[1]], label="state", **share)
    return [cpax, pax, crax, rax, csax, sax]


# =============================================================================
# Helper class
# =============================================================================


class TransposeGridSpec:
    """A grid spec that transposes subscripts.

    Parameters
    ----------
    gsp : GridSpec
        The Grid spec being transposed.

    When subscripted, `self[i,j] == gsp[j,i]`.
    """
    gsp: GridSpec

    def __init__(self, gsp: GridSpec) -> None:
        self.gsp = gsp

    def __getitem__(self, ind: ps.Inds) -> mpl.gridspec.SubplotSpec:
        return self.gsp[ind[::-1]]

    def __getattr__(self, attr: str) -> Any:
        return getattr(self.gsp, attr)


# =============================================================================
# Hint aliases
# =============================================================================
Figure = mpl.figure.Figure
GridSpec = mpl.gridspec.GridSpec
Disp = Union[mpl.lines.Line2D, ps.Image, mpl.text.Text, gp.GraphPlots]
TxHandle = Union[mpl.axes.Axes, mpl.text.Text]
AxList = List[mpl.axes.Axes]
