# -*- coding: utf-8 -*-
"""Creating video frames for optimal serial models
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib as mpl
import matplotlib.animation as mpla
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import numpy_linalg as la
import sl_py_tools.arg_tricks as ag
# import sl_py_tools.containers as cn
import sl_py_tools.iter_tricks as it
import sl_py_tools.matplotlib_tricks as mpt

# pyright: reportUndefinedVariable=false
import complex_synapse.builders as bld
import complex_synapse.graphs as gr
import complex_synapse.options as op
import complex_synapse.synapse_mem as sm
import complex_synapse.optimise.synapse_opt as so

mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# =============================================================================
# Calculation
# =============================================================================


def serial_eqp(param: la.lnarray) -> la.lnarray:
    """Steady state probabilities for a serial model

    Parameters
    ----------
    param : la.lnarray (2M-2,)
        Transition probabilities, in the order
        mat_01, mat_12, ..., mat_n-2,n-1,
        mat_10, mat_21, ..., mat_n-1,n-2.

    Returns
    -------
    eqp : la.lnarray (M,)
        Steady state probability
    """
    param = la.asarray(param).reshape((2, -1))
    pot, dep = tuple(param)
    nstw = (len(pot) + 1) // 2
    eqp = la.ones(2 * nstw)
    eqp[nstw + 1:] = np.cumprod(pot[nstw:] / dep[nstw:])
    eqp[:nstw - 1] = np.flip(np.cumprod(np.flip(dep[:nstw-1] / pot[:nstw-1])))
    eqp[nstw:] *= pot[nstw - 1] / dep[nstw - 1]
    eqp /= eqp.sum()
    return eqp


def serial_snr(param: la.lnarray, tau: la.lnarray) -> la.lnarray:
    """Signal-to-noise ratio, averaged over exponential distribution in time

    Parameters
    ----------
    params : la.lnarray (2M-2,)
        Transition probabilities, in the order
        mat_01, mat_12, ..., mat_n-2,n-1,
        mat_10, mat_21, ..., mat_n-1,n-2.
    tau : la.lnarray (T,)
        Mean of exponential distribution of time.

    Returns
    -------
    snr : la.lnarray (T,)
        Signal-to-noise ratio
    """
    topology = so.TopologyOptions(serial=True)
    obj = so.SynapseParamModel.from_params(
        param, binary=True, topology=topology).view(sm.SynapseMemoryModel)
    return obj.snr_exp_ave(tau)


def trim_params(param: la.lnarray) -> la.lnarray:
    """Set transition probability between orphaned states to zero

    Parameters
    ----------
    param : la.lnarray (2M-2,)
        Transition probabilities, in the order
        mat_01, mat_12, ..., mat_n-2,n-1,
        mat_10, mat_21, ..., mat_n-1,n-2.

    Returns
    -------
    pot : la.lnarray (M-1,)
        Transition probabilities for potentiation, in the order
        mat_01, mat_12, ..., mat_n-2,n-1,
    dep : la.lnarray (M-1,)
        Transition probabilities for depression, in the order
        mat_10, mat_21, ..., mat_n-1,n-2.
    """
    outward, inward, thresh = 1e-5, 1e-5, 1e-3
    # pylint: disable=unbalanced-tuple-unpacking
    param = la.array(param, copy=True).reshape((2, -1))
    pot, dep = tuple(param)
    nstw = (len(pot) + 1) // 2
    test = pot[nstw:] < thresh
    if test.any():
        ind = test.nonzero()[0][0] + nstw
        pot[ind:] = outward
        dep[ind:] = inward
    test = dep[:nstw-1] < thresh
    if test.any():
        ind = test.nonzero()[0][-1] + 1
        pot[:ind] = inward
        dep[:ind] = outward
    return pot, dep


# =============================================================================
# Fitter video options classes
# =============================================================================


# pylint: disable=too-many-ancestors
class VideoLabels(op.Options):
    """Tiles, axes labels, etc. for EnvelopeFig

    The individual options can be accessed as object instance attributes
    (e.g. `obj.name`) or as dictionary items (e.g. `obj['name']`) for both
    getting and setting.

    Parameters
    ----------
    snr : List[str]
        XY Axes labels for memory curve plots.
        By default: `[r"Time, $\tau$", "SNR"]`.
    env : List[str]
        Legend labels for theoretical envelope, numerical envelope and example.
        By default: `["Theory", "Numeric", "One model"]`.
    potdep : List[str]
        Y axis labels for transition probability bar charts.
        By default: `["Pot. prob.", "Dep. prob."]`.
    eqp : str
        Colourbar label for steady state distribution.
        By default: `"Equilibrium probability"`

    All parameters are optional keywords. Any dictionary passed as positional
    parameters will be popped for the relevant items. Keyword parameters must
    be valid keys, otherwise a `KeyError` is raised.
    """
    # Text for labels
    snr: List[str]
    env: List[str]
    potdep: List[str]
    eqp: str

    def __init__(self, *args, **kwds) -> None:
        self.snr = [r"Time, $\tau$", "SNR"]
        self.env = ["Theory", "Numeric", "One model"]
        self.potdep = ["Pot. prob.", "Dep. prob."]
        self.eqp = "Equilibrium probability"
        super().__init__(*args, **kwds)


# pylint: disable=too-many-ancestors
class VideoLayout(op.Options):
    """Tiles, axes labels, etc. for FitterVideo

    The individual options can be accessed as object instance attributes
    (e.g. `obj.name`) or as dictionary items (e.g. `obj['name']`) for both
    getting and setting.

    Parameters
    ----------
    height_ratios: Tuple[int, ...]
        Grid size ratios for various axes.
    width_ratios: Tuple[int, ...]
        Grid size ratios for various axes.
    ax_pos: Dict[str, Inds]
        which rows/cols for each Axes
    scale : float
        Controls figure size: inches per grid-ratio unit.

    All parameters are optional keywords. Any dictionary passed as positional
    parameters will be popped for the relevant items. Keyword parameters must
    be valid keys, otherwise a `KeyError` is raised.

    Notes
    -----
    Attributes can also be referenced and modified by subscripting with the
    attribute name. If the name is not found, it will search `sizes`.
    """
    map_attributes: op.Attrs = ('ax_pos',)
    # height/width ratios for rows/columns
    height_ratios: Tuple[int, ...]
    width_ratios: Tuple[int, ...]
    # row/column assignments for models/data
    ax_pos: Dict[str, Inds]
    scale: float

    def __init__(self, *args, **kwds) -> None:
        self.height_ratios = (2, 2, 1, 2, 2)
        self.width_ratios = (12, 9, 1)
        self.ax_pos = {'snr': np.s_[1:, 0], 'eqp': np.s_[2, 1],
                       'pot': np.s_[:2, 1], 'dep': np.s_[3:, 1],
                       'cbr': np.s_[:, 2], 'grf': np.s_[0, 0]}
        self.scale = 0.6
        super().__init__(*args, **kwds)

    def gspec_opts(self) -> Tuple[Tuple[float, float], Tuple[int, int],
                                  Dict[str, Tuple[int, ...]]]:
        """Options for creating a gridspec"""
        kwargs = {'height_ratios': self.height_ratios,
                  'width_ratios': self.width_ratios}
        gargs = tuple(map(len, kwargs.values()))
        fsiz = tuple(self.scale * sum(g) for g in reversed(kwargs.values()))
        return fsiz, gargs, kwargs


# pylint: disable=too-many-ancestors
class VideoOptions(op.MasterOptions, fallback='im_opt'):
    """Visual options for FitterVideo

    The individual options can be accessed as object instance attributes
    (e.g. `obj.name`) or as dictionary items (e.g. `obj['name']`) for both
    getting and setting.

    Parameters
    ----------
    txt : VideoLabels
        Title, axis, legend and colourbar labels.
    layout : VideoLayout
        Options for arrangement and sizes of video frames.
    graph : bool
        Transpose the layout of the video? Stored in `txt_opt` and `im_opt`
        under the key `'trn'`. By default: `False`
    ax_opt : Dict[str, Any]
         Options for `Axes`. See `sl_py_tools.matplotlib_tricks.clean_axes`.
         By default: `{'box': False, 'tight': False}`.
    im_opt : ImageOptions
        Options for heatmaps. By default: `ImageOptions(edgecolors='black')`.
    an_opt : AnimationOptions
        Options for animating the video.

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
    map_attributes: op.Attrs = ('txt', 'layout', 'graph', 'im_opt', 'ax_opt')
    # Text for labels
    txt: VideoLabels
    # Layout
    layout: VideoLayout
    # Graph
    graph: gr.GraphOptions
    # keyword options
    ax_opt: Dict[str, Any]
    im_opt: op.ImageOptions
    an_opt: op.AnimationOptions

    def __init__(self, *args, **kwds) -> None:
        self.txt = VideoLabels()
        self.layout = VideoLayout()
        self.graph = gr.GraphOptions()
        self.ax_opt = {'box': False, 'tight': False}
        self.im_opt = op.ImageOptions(edgecolors='black')
        self.an_opt = op.AnimationOptions()

        self.ax_opt.update(mpt.clean_axes_keys(self.ax_opt))
        super().__init__(*args, **kwds)
# pylint: enable=too-many-ancestors


# =============================================================================
# Envelope video figure
# =============================================================================


def _make_fig(opt: VideoLayout) -> Tuple[Figure, AxList]:
    """Create a figure with axes for an experiment/simulation fit.

    Returns
    -------
    fig : Figure
        Figure object
    axs : List[Axes]
        Axes for [SNR, Equilibrium, Colourbar, Potentiation, Depression, Graph]
    """
    fsiz, gsiz, kwargs = opt.gspec_opts()

    fig = plt.figure(figsize=fsiz, frameon=True, constrained_layout=True)
    gsp = fig.add_gridspec(*gsiz, **kwargs)
    fig.set_constrained_layout_pads(hpad=0, hspace=0)

    ax_snr = fig.add_subplot(gsp[opt.ax_pos['snr']])
    ax_eqp = fig.add_subplot(gsp[opt.ax_pos['eqp']])
    ax_pot = fig.add_subplot(gsp[opt.ax_pos['pot']], sharex=ax_eqp)
    ax_dep = fig.add_subplot(gsp[opt.ax_pos['dep']], sharex=ax_eqp)
    ax_cbr = fig.add_subplot(gsp[opt.ax_pos['cbr']])
    ax_grf = fig.add_subplot(gsp[opt.ax_pos['grf']])

    return fig, [ax_snr, ax_eqp, ax_cbr, ax_pot, ax_dep, ax_grf]


def _format_axes_snr(ax_snr: mpl.axes.Axes, opt: VideoOptions) -> None:
    """Format axes for SNR plots

    Parameters
    ----------
    ax_snr : Axes
        Axes for SNR
    """
    mpt_opts = opt.ax_opt.copy()
    mpt_opts.update(legendbox=False, box=True)

    ax_snr.set_xlabel(opt.txt.snr[0])
    ax_snr.set_ylabel(opt.txt.snr[1])
    ax_snr.legend(loc='lower left', edgecolor=ax_snr.get_facecolor())
    mpt.clean_axes(ax_snr, **mpt_opts)


def _format_axes_eqp(axs: AxList, imh: mpl.collections.QuadMesh, nst: int,
                     opt: VideoOptions) -> None:
    """Format axes for equilibrium distribution plots

    Parameters
    ----------
    axs : List[Axes]
        Axes for [Equilibrium, Colourbar]
    imh : QuadMesh
        Heatmap of equilibrium distribution.
    nst : int
        Number of states
    """
    mpt_opts = opt.ax_opt.copy()
    mpt_opts.update(box=False)

    ax_eqp, ax_cbr = axs
    fig = ax_cbr.get_figure()
    fig.colorbar(imh, cax=ax_cbr)
    ax_eqp.set_xticks(np.arange(nst))
    ax_eqp.set_xticklabels([''] * nst)
    ax_eqp.set_yticks([])
    ax_eqp.xaxis.set_ticks_position('both')
    mpt.clean_axes(ax_eqp, **mpt_opts)

    ax_cbr.set_ylabel(opt.txt.eqp)
    mpt.clean_axes(ax_cbr, **mpt_opts)


def _format_axes_bar(ax_bar: mpl.axes.Axes, nst: int, dep: bool,
                     opt: VideoOptions) -> None:
    """Format axes for potentiation probability plots

    Parameters
    ----------
    axs : Axes
        Axes for Potentiation
    nst : int
        Number of states
    dep : bool
        Is this the bar chart for depression transitions?
    """
    mpt_opts = opt.ax_opt.copy()
    mpt_opts.update(box=False)

    ax_bar.set_ylim([0, 1])
    if dep:
        ax_bar.invert_yaxis()
        ax_bar.xaxis.set_ticks_position('top')
        ax_bar.set_ylabel(opt.txt.potdep[1])
        ax_bar.spines['bottom'].set_visible(False)
    else:
        ax_bar.set_ylabel(opt.txt.potdep[0])
        ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)
    ax_bar.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax_bar.set_yticklabels(['0', '', '', '', '1'])
    ax_bar.set_xticks(np.arange(nst))
    ax_bar.set_xticklabels([''] * nst)
    mpt.clean_axes(ax_bar, **mpt_opts)


# =============================================================================
# Plotting
# =============================================================================


class GraphPlots:
    """Class for plotting model as a graph.

    Parameters
    ----------
    graph : nx.DiGraph
        Graph object describing model. Nodes have attributes `kind` and
        `value`. Edges have attributes `kind`, `value` and `pind` (if the
        model was a `SynapseParamModel`).
    opt : GraphOptions|None, optional
        Options for plotting the graph, by default `None -> GraphOptions()`.
    Other keywords passed to `opt` or `complex_synapse.graph.draw_graph`.
    """
    nodes: gr.NodePlots
    edges: gr.DirectedEdgeCollection
    pinds: np.ndarray
    opt: gr.GraphOptions

    def __init__(self, graph: nx.DiGraph,
                 opt: Optional[gr.GraphOptions] = None, **kwds) -> None:
        """Class for plotting model as a graph.

        Parameters
        ----------
        graph : nx.DiGraph
            Graph object describing model. Nodes have attributes `kind` and
            `value`.  Edges have attributes `kind`, `value` and `pind` (if the
            model was a `SynapseParamModel`).
        opt : GraphOptions|None, optional
            Options for drawing the graph, by default `None -> GraphOptions()`.
        Other keywords passed to `opt` or `complex_synapse.graph.draw_graph`.
        """
        self.opt = ag.default_eval(opt, gr.GraphOptions)
        self.opt.pop_my_args(kwds)
        self.nodes, self.edges = gr.draw_graph(graph, opts=self.opt, **kwds)
        if gr.has_edge_attr(graph, 'pind'):
            self.pinds = gr.edge_attr_vec(graph, 'pind')
        else:
            self.pinds = np.arange(len(self.edges))

    def update(self, params: la.lnarray, peq: Optional[la.lnarray] = None
               ) -> None:
        """Update plots.

        Parameters
        ----------
        params : la.lnarray (E,)
            Transition probabilities, for edge line widdths.
        peq : None|la.lnarray (N,), optional
            Equilibrium distribution,for nodes sizes (area), by default `None`
            -> calculate from `params`.
        """
        params = la.asarray(params).ravel()
        peq = ag.default_eval(peq, lambda: serial_eqp(params))
        self.nodes.set_sizes(peq * self.opt.size)
        self.edges.set_widths(params[self.pinds] * self.opt.width)
        self.edges.set_node_sizes(peq * self.opt.size)

    def update_from(self, graph: nx.DiGraph) -> None:
        """Update plots using a graph object.

        Parameters
        ----------
        graph : nx.DiGraph
            Graph object describing model. Nodes have attributes `kind` and
            `value`.  Edges have attributes `kind`, `value`.
        """
        params = gr.edge_attr_vec(graph, 'value')
        peq = gr.node_attr_vec(graph, 'value')
        self.update(params, peq)

    @classmethod
    def from_data(cls, axs: mpl.axes.Axes, model: la.lnarray,
                  peq: Optional[la.lnarray] = None,
                  weight: Optional[la.lnarray] = None, **kwds) -> GraphPlots:
        """Plot data from model parameters to make an instance

        Parameters
        ----------
        axs : mpl.axes.Axes
            Axes on which we plot
        model : la.lnarray (2M-2,)
            Parameters of serial model
        peq : None|la.lnarray (M,), optional
            Equilbrium distribution, by default `None` -> `serial_eqp(model)`.
        weight : None|la.lnarray (M,), optional
            Maps state to synaptic weight, by default `None` -> `binary_weights`.

        Returns
        -------
        obj : GraphPlots
            A `GraphPlots` instance containing plotted objects.
        """
        opts: gr.GraphOptions = kwds.pop('opts', gr.GraphOptions())
        opts.pop_my_args(kwds)
        model = la.asarray(model).reshape((2, -1))
        peq = ag.default_eval(peq, lambda: serial_eqp(model))
        weight = ag.default_eval(weight, lambda: bld.binary_weights(len(peq)))
        graph = gr.param_to_graph(model, peq, weight, opts.topology)
        return cls(graph, axs=axs, opt=opts, **kwds)


class ModelPlots:
    """The plots associated with a model

    Parameters
    ----------
    lns : Sequence[mpl.lines.Line2D]
        The plots for [model's memory curve, time at which it's optimal].
    imh : mpl.image.AxesImage
        Heatmap of equilibrium distribution.
    brs : List[mpl.container.BarContainer]
        Bar charts of transition rates under [potentiation, depression].
    graph : GraphPlots
        Plot of graph describing the model.
    """
    lnh: mpl.lines.Line2D
    vln: mpl.lines.Line2D
    # imh: mpl.image.AxesImage
    imh: mpl.collections.QuadMesh
    brp: mpl.container.BarContainer
    brd: mpl.container.BarContainer
    graph: GraphPlots

    def __init__(self, lns: Sequence[mpl.lines.Line2D],
                 imh: mpl.image.AxesImage,
                 brs: List[mpl.container.BarContainer],
                 graph: GraphPlots) -> None:
        """The plots associated with a model

        Parameters
        ----------
        lns : Sequence[mpl.lines.Line2D]
            The plots for [model's memory curve, time at which it's optimal].
        imh : mpl.image.AxesImage
            Heatmap of equilibrium distribution.
        brs : List[mpl.container.BarContainer]
            Bar charts of transition rates under [potentiation, depression].
        graph : GraphPlots
            Plot of graph describing the model.
        """
        self.lnh = lns[0]
        self.vln = lns[1]
        self.imh = imh
        self.brp = brs[0]
        self.brd = brs[1]
        self.graph = graph

    def update_snr(self, snr: la.lnarray, tau: float):
        """Update model's SNR plot

        Parameters
        ----------
        snr : la.lnarray
            SNR of model
        tau : float
            Time at which this model is optimal
        """
        self.lnh.set_ydata(snr)
        self.vln.set_xdata([tau, tau])

    def update_eqp(self, eqp: la.lnarray):
        """Update model's wquilibrium distribution heatmap

        Parameters
        ----------
        eqp : la.lnarray (M,)
            Equilibrium probability distribution
        """
        # self.imh.set_data(eqp.r)
        self.imh.set_array(eqp)

    def update_potdep(self, pot: la.lnarray, dep: la.lnarray):
        """Update model's potentiation bar chart

        Parameters
        ----------
        pot : la.lnarray (M-1,)
            Potentiation transition probability, M_i,i+1
        dep : la.lnarray
            Depression transition probability, M_i+1,i
        """
        for rect, height in zip(self.brp, pot):
            rect.set_height(height)

        for rect, height in zip(self.brd, dep):
            rect.set_height(height)

    def update_plots(self, time: la.lnarray, model: la.lnarray, tau: float):
        """Update plots with new model

        Parameters
        ----------
        time : la.lnarray (T,)
            Array of times for memory curve
        model : la.lnarray (2M-2,)
            Parameters of serial model.
        tau : float
            Time at which `model` is optimal.
        """
        snr = serial_snr(model, time)
        eqp = serial_eqp(model)
        pot, dep = trim_params(model)

        self.update_snr(snr, tau)
        self.update_eqp(eqp)
        self.update_potdep(pot, dep)
        self.graph.update(np.r_[pot, dep], eqp)

    @classmethod
    def from_data(cls, axs: AxList, time: la.lnarray, model: la.lnarray,
                  tau: float, **kwds) -> ModelPlots:
        """Plot data to make an instance

        Parameters
        ----------
        axs : AxList
            Array of Axes on which we plot.
        time : la.lnarray (T,)
            Array of times for memory curve
        model : la.lnarray (2M-2,)
            Parameters of serial model.
        tau : float
            Time at which `model` is optimal.

        Returns
        -------
        obj : ModelPlots
            Object holding the plot objects associated with `model`.
        """
        opt: VideoOptions = kwds.pop('opt', VideoOptions())
        opt.pop_my_args(kwds)

        snr = serial_snr(model, time)
        eqp = serial_eqp(model)
        pot_dep = trim_params(model)
        njmp = len(pot_dep[0])


        lns = axs[0].loglog(time, snr, label=opt.txt.env[2])
        lns.append(axs[0].axvline(tau))

        # imh = axs[1].imshow(eqp.r, norm=mpl.colors.Normalize(0, 1))
        imh = axs[1].pcolormesh(np.arange(-0.5, njmp+1), [0, 1], eqp.r,
                                **opt.im_opt)

        graph = GraphPlots.from_data(axs[5], pot_dep, eqp, opts=opt.graph,
                                     **kwds)

        brs = [axs[3].bar(np.arange(0.5, njmp), pot_dep[0],
                          color=graph.opt.edge_cmap(1.0)),
               axs[4].bar(np.arange(0.5, njmp), pot_dep[1],
                          color=graph.opt.edge_cmap(0.0))]

        return cls(lns, imh, brs, graph)


class EnvelopeFig:
    """Data and figure objects for an envelope plot

    Parameters
    ----------
    rate : la.lnarray (T,)
        Rate parameter of Laplace transform (stored as `time = 1/rate`)
    env_th : la.lnarray (T,)
        Theoretical Laplacian envelope (stored as exponential running average,
        `env_th -> env_th * rate`)
    env_num : la.lnarray (T,)
        Numerical Laplacian envelope (stored as exponential running average,
        `env_num -> env_num * rate`)
    models : la.lnarray (T,2N-2)
        Parameters of serial models that form envelope, in the order
        mat_01, mat_12, ..., mat_n-2,n-1,
        mat_10, mat_21, ..., mat_n-1,n-2.
    """
    time: la.lnarray
    env_th: la.lnarray
    env_num: la.lnarray
    models: la.lnarray
    fig: Figure
    model_plots: ModelPlots
    opt: VideoOptions

    def __init__(self, rate: la.lnarray, env_th: la.lnarray,
                 env_num: la.lnarray, models: la.lnarray, **kwds) -> None:
        """Data and figure objects for an envelope plot

        Parameters
        ----------
        rate : la.lnarray (T,)
            Rate parameter of Laplace transform (stored as `time = 1/rate`)
        env_th : la.lnarray (T,)
            Theoretical Laplacian envelope (stored as exponential running
            average, `env_th -> env_th * rate`)
        env_num : la.lnarray (T,)
            Numerical Laplacian envelope (stored as exponential running
            average, `env_num -> env_num * rate`)
        models : la.lnarray (T,2M-2)
            Parameters of serial models that form envelope, in the order
            mat_01, mat_12, ..., mat_n-2,n-1,
            mat_10, mat_21, ..., mat_n-1,n-2.
        """
        self.opt = kwds.pop('opt', VideoOptions())
        self.opt.pop_my_args(kwds)
        self.time = 1 / rate
        self.env_th = env_th * rate
        self.env_num = env_num * rate
        self.models = models
        if self.time[0] > self.time[-1]:
            self.flip()

        self.fig, axs = _make_fig(self.opt.layout)

        axs[0].loglog(self.time, self.env_th, label=self.opt.txt.env[0])
        axs[0].loglog(self.time, self.env_num, label=self.opt.txt.env[1])

        self.model_plots = ModelPlots.from_data(axs, self.time, self.models[0],
                                                self.time[0], opt=self.opt)

        _format_axes_snr(axs[0], self.opt)
        _format_axes_eqp(axs[1:3], self.model_plots.imh, self.nstate, self.opt)
        _format_axes_bar(axs[3], self.nstate, False, self.opt)
        _format_axes_bar(axs[4], self.nstate, True, self.opt)
        axs[5].set_frame_on(False)

        # fig.set_constrained_layout(True)
        self.fig.set_constrained_layout_pads(hpad=0, hspace=0)
        plt.draw()

    @property
    def nstate(self) -> int:
        """Number of states"""
        return self.models.shape[-1] // 2 + 1

    @property
    def ntime(self) -> int:
        """Number of time points"""
        return self.models.shape[0]

    def update(self, ind: int):
        """Change which model we're plotting.

        Parameters
        ----------
        ind : int
            Which row of `models` should we plot?
        """
        self.model_plots.update_plots(self.time, self.models[ind],
                                      self.time[ind])

    def flip(self) -> None:
        """Flip arrays along time axis
        """
        self.time = np.flip(self.time, axis=0)
        self.env_th = np.flip(self.env_th, axis=0)
        self.env_num = np.flip(self.env_num, axis=0)
        self.models = np.flip(self.models, axis=0)


# =============================================================================
# Animation
# =============================================================================


def animate(env_fig: EnvelopeFig, **kwargs) -> mpla.FuncAnimation:
    """Animate a fitter video

    Parameters
    ----------
    env_fig : EnvelopeFig
        The envelope figure we will animate. Default options are in
        `env_fig.opt.an_opt`.
    Other keywords
        Passed to `FuncAnimation`.

    Returns
    -------
    ani : FuncAnimation
        The animation.

    To view the animation call `plt.show()`. To save a video, call `ani.save()`
    (see https://matplotlib.org/api/_as_gen/matplotlib.animation.Animation.html
    #matplotlib.animation.Animation.save)
    """
    kwargs.setdefault('init_func', None)
    kwargs.setdefault('frames', it.erange(env_fig.ntime))
    opt = env_fig.opt.an_opt.copy()
    opt.update(kwargs)
    return mpla.FuncAnimation(env_fig.fig, env_fig.update, **opt)


# =============================================================================
# Hint aliases
# =============================================================================
Figure = mpl.figure.Figure
AxList = List[mpl.axes.Axes]
Ratios = Tuple[List[int], List[int]]
Ind = Union[int, slice]
Inds = Tuple[Ind, ...]
