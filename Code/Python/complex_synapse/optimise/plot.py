"""Creating figures for Envelopes
"""
from __future__ import annotations
import typing as ty
from typing import Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import numpy_linalg as la
import sl_py_tools.arg_tricks as ag
import sl_py_tools.matplotlib_tricks as mplt
import sl_py_tools.dict_tricks as dt
import sl_py_tools.numpy_tricks.markov as ma
import sl_py_tools.graph_plots as gp
import sl_py_tools.graph_tricks as gt

import complex_synapse.builders as bld
import complex_synapse.optimise.optimise as cso
import complex_synapse.optimise.shorten as sh
import complex_synapse.synapse_mem as sm

np.set_printoptions(precision=3, suppress=False, linewidth=90)
mplt.rc_fonts()
mplt.rc_colours()


# =============================================================================
# Graph plotting
# =============================================================================


class GraphPlots(gp.GraphPlots):
    """Class for plotting model as a graph.

    Parameters
    ----------
    graph : nx.DiGraph
        Graph object describing model. Nodes have attributes `kind` and
        `value`. Edges have attributes `kind`, `value` and `pind` (if the
        model was a `SynapseParamModel`).
    opts : GraphOptions|None, optional
        Options for plotting the graph, by default `None -> GraphOptions()`.
    Other keywords passed to `opt` or `complex_synapse.graph.draw_graph`.
    """

    def update(self, edge_vals: la.lnarray,
               node_vals: ty.Optional[la.lnarray] = None) -> None:
        """Update plots.

        Parameters
        ----------
        params : la.lnarray (E,)
            Transition probabilities, for edge line widdths.
        peq : None|la.lnarray (N,), optional
            Equilibrium distribution,for nodes sizes (area), by default `None`
            -> calculate from `params`.
        """
        edge_vals = la.asarray(edge_vals).ravel()
        node_vals = ag.default_eval(node_vals, lambda: serial_eqp(edge_vals))
        super().update(edge_vals, node_vals)

    @classmethod
    def from_data(cls, axs: mpl.axes.Axes, model: la.lnarray,
                  peq: ty.Optional[la.lnarray] = None,
                  weight: ty.Optional[la.lnarray] = None,
                  **kwds) -> GraphPlots:
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
            Synaptic weight of states, by default `None` -> `binary_weights`.

        Returns
        -------
        obj : GraphPlots
            A `GraphPlots` instance containing plotted objects.
        """
        opts: gp.GraphOptions = kwds.pop('opts', gp.GraphOptions())
        opts.pop_my_args(kwds)
        model = la.asarray(model).reshape((2, -1))
        peq = ag.default(peq, serial_eqp(model))
        weight = ag.default(weight, bld.binary_weights(len(peq)))
        graph = gt.param_to_graph(model, peq, weight, None, opts.topology)
        return cls(graph, axs=axs, opts=opts, **kwds)


def serial_eqp(param: la.lnarray) -> la.lnarray:
    """Steady state probabilities for a serial model

    Parameters
    ----------
    param : la.lnarray (2,M-1)
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


# =============================================================================
# ## All data
# =============================================================================


def all_data(nst: int, sss: Array, **opts) -> ty.Tuple[Arrays, Arrays]:
    """All of the data generation

    Parameters
    ----------
    nst : int
        Number of states
    sss : Array (T,)
        Rate parameter of Laplace transform of SNR curve

    Returns
    -------
    envs : Tuple[Array] (3,)(T,)
        Envelopes: maximum Laplace(SNR) at varios times.
    models : Tuple[Array] (3,)(T,Q)
        Parameters of the models that achieve the envelopes.
        `Q = 2M(M-1) or 2(M-1)`.
    Tuple members:-
        * Arbritrary topology, normal problem.
        * Serial topology, normal problem.
        * Arbritrary topology, shifted problem.
    """
    env_gen, mods_gen = cso.optim_laplace_range(sss, nst, **opts)
    with dt.updated(opts, serial=True) as opt:
        env_srl, mods_srl = cso.optim_laplace_range(sss, nst, **opt)
    with dt.updated(opts, cond_lim=True, maker=cso.shifted_problem) as opt:
        env_shf, mods_shf = cso.optim_laplace_range(sss, nst, **opt)
    return (env_gen, env_srl, env_shf), (mods_gen, mods_srl, mods_shf)


# =============================================================================
# ## Generic
# =============================================================================


def mem_plot(sss: Array, nsyn: float = 1, **mem: Array
             ) -> ty.Tuple[Figure, Axes, Line]:
    """Plot memory curves

    Parameters
    ----------
    sss : Array
        1/time
    nsyn : float, optional
        Number of synapses, by default 1.
    mem : Dict[str, Array]
        Laplace transforms of memory curves to plot, keyed by plot label

    Returns
    -------
    fig : mpl.figure.Figure
        Figure containing `axs`.
    axs : mpl.axes.Axes
        Axes containing plots.
    """
    pref = np.sqrt(nsyn)
    fig, axs = plt.subplots()
    lns = []
    for label, mem in mem.items():
        lns.extend(axs.loglog(1/sss, mem * sss * pref, label=label))
    axs.set_xlabel(r"Time, $r\tau$")
    axs.set_ylabel("SNR")
    return fig, axs, lns


# =============================================================================
# Examples
# =============================================================================


def ser_casc_plot(nst: int, jmp: float, sss: la.lnarray, nsyn: float = 1
                  ) -> mpl.figure.Figure:
    """Example plot

    Parameters
    ----------
    nst : int
        Number of states
    jmp : float
        Model parameter
    time : la.lnarray
        Vector of times to evaluate SNR
    nsyn : float, optional
        Number of synapses, by default 1.

    Returns
    -------
    fig : Figure
        The figure object containing the plot.
    """
    nst = ag.default(nst, 10)
    sss = ag.default(sss, la.geomspace(1e-4, 10, 50))
    jmp = ag.default(jmp, 0.7)

    serial = sm.SynapseMemoryModel.build(bld.build_serial, nst, jmp=jmp)
    cascade = sm.SynapseMemoryModel.build(bld.build_cascade, nst, jmp=jmp)

    serial_snr = serial.snr_laplace(sss)
    cascade_snr = cascade.snr_laplace(sss)

    fig, axs, _ = mem_plot(sss, nsyn, serial=serial_snr, cascade=cascade_snr)
    axs.set_xlim(1e-1, 1e4)
    axs.legend(loc="lower left")
    mplt.clean_axes(axs)
    return fig


# =============================================================================
# ## Theory
# =============================================================================


def theory_plot(nst: int, sss: Array, nsyn: float = 1) -> Figure:
    """Theoretical envelope plot

    Parameters
    ----------
    nst : int
        Number of states
    sss : Array (T,)
        Rate parameter of Laplace transform of SNR curve

    Returns
    -------
    fig : Figure
        The figure object containing the plot.
    """
    env_th = cso.proven_envelope_laplace(sss, nst)
    fig, axs, _ = mem_plot(sss, nsyn, theory=env_th)
    t_pts = np.array([1e-1, 0.9e1, 1.1e1, 1e4])
    axs.loglog(t_pts[:2], np.ones(2), 'k:')
    axs.loglog(t_pts[-2:], (nst-1) / t_pts[-2:], 'k:')
    axs.set_xlim(1e-1, 1e3)
    axs.set_ylim(1e-2 * np.sqrt(nsyn), 3 * np.sqrt(nsyn))
    ind = 30
    axs.text(t_pts[1]/2, 1.1, r"$\sqrt{N}$", fontsize=18,
             ha='right', va='bottom')
    axs.text(t_pts[-2]*4, (nst-1) / (t_pts[-2] * 4),
             r"$\frac{\sqrt{N}(M-1)}{r\tau}$", fontsize=20,
             ha='left', va='bottom')
    axs.text(1/sss[ind], env_th[ind] * sss[ind],
             r"$\frac{\sqrt{N}(M-1)}{r\tau + (M-1)}$", fontsize=24,
             ha='right', va='top')
    mplt.clean_axes(axs)
    return fig


def equilibrium_plot(nst: int, sss: Array, nsyn: float = 1) -> Figure:
    """Theoretical envelope plot

    Parameters
    ----------
    nst : int
        Number of states
    sss : Array (T,)
        Rate parameter of Laplace transform of SNR curve
    nsyn : float, optional
        Number of synapses, by default 1.

    Returns
    -------
    fig : Figure
        The figure object containing the plot.
    """
    env_th = cso.proven_envelope_laplace(sss, nst)
    env_eq = cso.equlibrium_envelope_laplace(sss, nst)
    fig, axs, lns = mem_plot(sss, nsyn, **{"general": env_th,
                                           "det. bal.": env_eq})
    lns[0].set_ls(":")
    lns[1].set_c(lns[0].get_c())
    # t_pts = np.array([1e-1, 0.9e1, 1.1e1, 1e4])
    # axs.loglog(t_pts[:2], np.ones(2), 'k:')
    # axs.loglog(t_pts[-2:], (nst-1) / t_pts[-2:], 'k:')
    axs.set_xlim(1e-1, 1e3)
    ylim = axs.set_ylim(1e-2 * np.sqrt(nsyn), 3 * np.sqrt(nsyn))
    s_two, s_sticky = 0.5, 2 / (nst-1)**2
    _annotate_axis(axs, nst, ylim[0], [s_two, s_sticky],
                ["2", r"\frac{(M-1)^2}{2}"])

    ind = 40
    axs.text(1/sss[ind], env_eq[ind] * sss[ind],
             r"$\frac{2\sqrt{N}}{2+r\tau}$", fontsize=20,
             ha='right', va='top')
    ind = 30
    axs.text(1/sss[ind], env_eq[ind] * sss[ind],
             r"$\sqrt{\frac{N}{2r\tau}}$", fontsize=20,
             ha='right', va='top')
    ind = 20
    axs.text(1/sss[ind], env_eq[ind] * sss[ind],
             r"$\frac{2\sqrt{N}(M-1)}{(M-1)^2 + 2r\tau}$", fontsize=20,
             ha='left', va='bottom')
    axs.legend(loc="upper right")
    mplt.clean_axes(axs)
    return fig


# =============================================================================
# ## Optimisation
# =============================================================================


def optim_plot(nst: int, sss: Array, env_gen: Array, env_srl: Array,
               nsyn: float = 1) -> Figure:
    """Optimisation plot.

    Comparison of theoretical envelope and numerical optimisation over models
    with aribtrary or serial topologies.

    Parameters
    ----------
    nst : int
        Number of states
    sss : Array (T,)
        Rate parameter of Laplace transform of SNR curve
    env_gen : Array (T,)
        Envelope for arbitrary topologies.
    env_srl : Array (T,)
        Envelope for serial topologies only.
    nsyn : float, optional
        Number of synapses, by default 1.

    Returns
    -------
    fig : Figure
        The figure object containing the plot.
    """
    env_th = cso.proven_envelope_laplace(sss, nst)
    env_eq = cso.equlibrium_envelope_laplace(sss, nst)
    fig, axs, lns = mem_plot(sss, nsyn, **{"theory: general": env_th,
                                           "conjecture: det. bal.": env_eq,
                                           'numeric: general': env_gen,
                                           'numeric: serial': env_srl})
    lns[0].set_ls(":")
    lns[2].set_c(lns[1].get_c())
    lns[1].set_c(lns[0].get_c())
    axs.set_xlim(1/sss[-1], 1/sss[0])
    axs.legend(loc="lower left")
    mplt.clean_axes(axs)
    return fig


# =============================================================================
# ### Shifted optimisation
# =============================================================================


def shift_plot(sss: Array, env_gen: Array, env_shf: Array, nsyn: float = 1
               ) -> Figure:
    """Shifted optimisation plot.

    Comparison of numerical optimisation of the normal and shifted problems.

    Parameters
    ----------
    sss : Array (T,)
        Rate parameter of Laplace transform of SNR curve
    env_gen : Array (T,)
        Envelope from normal problem.
    env_shf : Array (T,)
        Envelope from shifted problem.
    nsyn : float, optional
        Number of synapses, by default 1.

    Returns
    -------
    fig : Figure
        The figure object containing the plot.
    """
    fig, axs, _ = mem_plot(sss, nsyn, normal=env_gen, shifted=env_shf)
    axs.set_xlim(1/sss[-1], 1/sss[0])
    axs.legend(loc="lower left")
    mplt.clean_axes(axs)
    return fig


# =============================================================================
# ## Heuristic
# =============================================================================


def heuristic_plot(nst: int, sss: Array, env_srl: Array, models_srl: Array,
                   nsyn: float = 1, **kwds) -> Figure:
    """Heuristic envelope plot

    Parameters
    ----------
    nst : int
        Number of states
    sss : Array (T,)
        Rate parameter of Laplace transform of SNR curve
    env_srl : Array (T,)
        Numeric envelope.
    models_srl : Arrays (3,)(T,2(M-1))
        Parameters of the serial models that achieve the envelope.
    nsyn : float, optional
        Number of synapses, by default 1.

    Returns
    -------
    fig : Figure
        The figure object containing the plot.
    """
    kwds.setdefault('model_inds', (43, 25, 9))
    kwds.setdefault('model_above', (0, 0, 0))
    kwds.setdefault('model_siz', (2.3, 1.7))
    kwds.setdefault('note_inds', (43, 30, 9))
    kwds.setdefault('note_above', (1, 1, 1))
    kwds.setdefault('legend_loc', "upper right")
    kwds.setdefault('legendfontscale', 0.8)
    kwds.setdefault('legendbox', False)
    ax_opts = mplt.AxesOptions(kwds)

    env_th = cso.proven_envelope_laplace(sss, nst)
    env_eq = cso.equlibrium_envelope_laplace(sss, nst)
    env_heu = cso.heuristic_envelope_laplace(sss, nst)
    envs = {"theory": env_th, "det. bal.": env_eq}

    fig, axs, lns = mem_plot(sss, nsyn, **envs, numeric=env_srl,
                             heuristic=env_heu)
    lns[0].set_ls(":")
    lns[2].set_c(lns[1].get_c())
    lns[1].set_c(lns[0].get_c())
    axs.set_xlim(1/sss[-1], 1e4)
    ylim = axs.set_ylim(2e-3 * np.sqrt(nsyn), 1.1 * np.sqrt(nsyn))

    _annotate_curves(axs, 1/sss, env_heu * sss, models_srl, **kwds)
    _annotate_axis(axs, nst, ylim[0], [sh.s_star(2), sh.s_star(nst)],
                   ["0.73", "0.22M^2"])

    axs.legend(loc=kwds['legend_loc'], edgecolor=axs.get_facecolor())
    mplt.clean_axes(axs, **ax_opts)
    return fig


def _annotate_curves(axs: mpl.axes.Axes, time: Array, env: Array,
                     models: Array, **kwds) -> None:
    """Add text and graphs along envelope
    """
    args = ({'ha': 'right', 'va': 'top'}, {'ha': 'left', 'va': 'bottom'})
    txt_opt = {'size': 24}
    notes = (r"$\sqrt{N}$",
             r"$\frac{\sqrt{N}(0.54)}{\sqrt{r\tau}}$",
             r"$\frac{\sqrt{N}(M-1)}{r\tau}$")

    for ind, above in zip(kwds['model_inds'], kwds['model_above']):
        add_graph(axs, [time[ind-1], env[ind], *kwds['model_siz']],
                  models[ind], **args[above])

    for ind, above, note in zip(kwds['note_inds'], kwds['note_above'], notes):
        axs.text(time[ind], env[ind-1], note, **txt_opt, **args[above])


def _annotate_axis(axs: mpl.axes.Axes, nst: int, y_0: float,
                   rates: Sequence[float], labels: Sequence[str]) -> None:
    """Add text and graphs along time axis
    """
    t_two = 1 / rates[0]
    t_sticky = 1 / rates[1]
    txt_opt = {'ha': 'center', 'va': 'bottom', 'size': 16}

    axs.axvline(x=t_two, color='k', linestyle=':')
    axs.axvline(x=t_sticky, color='k', linestyle=':')

    axs.text(t_two, y_0, r"$r\tau = {}$".format(labels[0]), **txt_opt)
    axs.text(t_sticky, y_0 * 1.1, r"$r\tau = {}$".format(labels[1]), **txt_opt)


def add_graph(axs: mpl.axes.Axes, bounds: ty.Sequence[float], model: Array,
              horizontalalignment: str = 'left',
              verticalalignment: str = 'bottom',
              **kwds) -> ty.Tuple[mpl.axes.Axes, GraphPlots]:
    """Draw the graph of a model on some axes

    Parameters
    ----------
    axs : mpl.axes.Axes
        The axes upon which we draw.
    bounds : Sequence[float]
        `[left, bottom, width, height]`, where `width` and `height` are
        multiplicative when the corresponding axis is logarithmic.
    model : Array
        The parameters of the serial model we draw.

    Returns
    -------
    axin : mpl.axes.Axes
        The inset axes on which the graph is drawn.
    graph : GraphPlots
        An object containing all the plot objects for the graph.
    """
    hal = kwds.pop('ha', horizontalalignment)
    val = kwds.pop('va', verticalalignment)
    model = trim_params(model)
    nst = ma.params.num_state(model.shape[-1], serial=True)
    bounds, scale = scale_bounds(bounds, nst, (3, 2), axs)
    bounds = logify(bounds, axs, hal, val)
    kwds.setdefault('edges.mut_scale', 0.1 * scale)
    kwds.setdefault('edges.mult', 2 * scale)
    kwds.setdefault('nodes.mult', 100 * scale**2)
    opts = gp.GraphOptions(kwds)

    axin = axs.inset_axes(bounds, transform=axs.transData, frame_on=False)
    graph = GraphPlots.from_data(axin, model, opts=opts)
    axin.set_xlim(-0.5, nst + 0.5)
    axin.set_ylim(-0.75, 0.75)
    axin.set_clip_on(False)
    # con = ConnectionPatch(xyA=(1.0, 0.5), coordsA=axin.transAxes,
    #                       xyB=xy, coordsB=axs.transData)
    return axin, graph


def logify(bounds: ty.Sequence[float], axs: mpl.axes.Axes,
           horizontalalignment: str = 'left',
           verticalalignment: str = 'bottom', **kwds) -> ty.List[float]:
    """Adjust bounds for logarithmic axes

    Parameters
    ----------
    bounds : Sequence[float]
        `[left, bottom, width, height]`, where `width` and `height` are
        multiplicative when the corresponding axis is logarithmic.
    axs : mpl.axes.Axes
        The axes to which the bounds refer. Used to check for logarithmic axes.

    Returns
    -------
    bounds : List[float]
        New bounds, with positive additive `width` and `height`
    """
    hal = kwds.pop('ha', horizontalalignment)
    val = kwds.pop('va', verticalalignment)
    x_0, x_1 = _logify(bounds[::2], axs.get_xscale(), hal)
    y_0, y_1 = _logify(bounds[1::2], axs.get_yscale(), val)
    return [x_0, y_0, x_1 - x_0, y_1 - y_0]


def _logify(bounds: ty.Sequence[float], scale: str, align: str) -> Array:
    """Adjust pair of bounds for log scales"""
    return _BWD[scale](np.sort(np.cumsum(_align(_FWD[scale](bounds), align))))


def _align(bounds: Array, align: str) -> Array:
    """Set alignment, modifies in place and returns"""
    bounds[1] = abs(bounds[1])
    bounds[0] -= _ALGN[align] * bounds[1]
    return bounds


def scale_bounds(bounds: ty.Sequence[float], nst: int, std: ty.Sequence[float],
                 axs: mpl.axes.Axes) -> ty.Tuple[ty.List[float], float]:
    """Scale bounds for the number of states we have

    Parameters
    ----------
    bounds : Sequence[float]
        `[left, bottom, width, height]`
    nst : int
        Number of states minus one
    std : ty.Sequence[float]
        The `(width, height)`, for which node/edge options are designed.
    axs : mpl.axes.Axes
        The axes to which the bounds refer. Used to check for logarithmic axes.

    Returns
    -------
    bounds : List[float]
        Modified input, with width adjusted for nst`
    scale : float
        Scale factor to be applied to node/edge options.
    """
    bounds = list(bounds)
    funx, funy = _FWD[axs.get_xscale()], _FWD[axs.get_yscale()]
    scalex = funx(bounds[2]) / funx(std[0])
    scaley = funy(bounds[3]) / funy(std[1])
    ifun = _BWD[axs.get_xscale()]
    # half of width is end-node buffer, rest for transition, need nst-1 more
    bounds[2] = ifun(((nst - 1) // 2 + 1) * funx(bounds[2]))
    return bounds, np.sqrt(scalex * scaley)


def trim_params(param: la.lnarray) -> la.lnarray:
    """Remove transitions between orphaned states

    Parameters
    ----------
    param : la.lnarray (2M-2,)
        Transition probabilities, in the order
        mat_01, mat_12, ..., mat_n-2,n-1,
        mat_10, mat_21, ..., mat_n-1,n-2.

    Returns
    -------
    param : la.lnarray (2M'-2,)
        Transition probabilities, in the order
        mat_01, mat_12, ..., mat_n-2,n-1,
        mat_10, mat_21, ..., mat_n-1,n-2.
    """
    thresh = 1e-3
    # pylint: disable=unbalanced-tuple-unpacking
    pot, dep = la.array(param, copy=True).reshape((2, -1))
    nstw = (len(pot) + 1) // 2
    test = pot[nstw:] < thresh
    if test.any():
        ind = test.nonzero()[0][0] + nstw
        pot, dep = pot[:ind], dep[:ind]
    test = dep[:nstw-1] < thresh
    if test.any():
        ind = test.nonzero()[0][-1] + 1
        pot, dep = pot[ind:], dep[ind:]
    return np.concatenate((pot, dep))


# =============================================================================
# ## Diagnosis and treatment
# =============================================================================


def reoptim(env_shf: Array, env_gen: Array, models_shf: Array, sss: Array,
            **options) -> Arrays:
    """Reattempt failed optimisations

    Parameters
    ----------
    env_shf : Array (T,)
        Envelope from shifted problem.
    env_gen : Array (T,)
        Envelope from normal problem.
    models_shf : Arrays (3,)(T,2M(M-1))
        Parameters of the shifted models that achieve the envelope.
    sss : Array (T,)
        Rate parameter of Laplace transform of SNR curve

    Returns
    -------
    env_shf : Array (T,)
        Envelope from shifted problem.
    models_shf : Arrays (3,)(T,2M(M-1))
        Parameters of the shifted models that achieve the envelope.
    """
    bad_inds, = (env_shf / env_gen < 0.6).nonzero()
    options['repeats'] = 20
    env_shf, models_shf = cso.reoptim_laplace_range(
        bad_inds, sss, env_shf, models_shf,
        maker=cso.shifted_problem, cond_thresh=1e3, cond_lim=False, **options
    )
    return env_shf, models_shf


def cond_plot(sss: Array, models_shf: Array, env_shf: Array, env_gen: Array
              ) -> mpl.figure.Figure:
    """Plot discrepancy vs condition number

    Parameters
    ----------
    sss : Array (T,)
        Rate parameter of Laplace transform of SNR curve
    models_shf : Arrays (3,)(T,2M(M-1))
        Parameters of the shifted models that achieve the envelope.
    env_shf : Array (T,)
        Envelope from shifted problem.
    env_gen : Array (T,)
        Envelope from normal problem.

    Returns
    -------
    fig : Figure
        The figure object containing the plot.
    """
    conds = cso.check_cond_range(sss, models_shf)

    fig, axs = plt.subplots()
    # good = senvelope_gen < envelope_gen
    # bad = np.logical_not(good)
    # ax.semilogx(conds[good], senvelope_gen[good] / envelope_gen[good], '+',
    #             conds[bad], senvelope_gen[bad] / envelope_gen[bad], '+')
    axs.semilogx(conds, env_shf / env_gen, '+')
    axs.axhline(y=1, color='k', linestyle=':')
    axs.set_xlabel('Condition number of $Z$')
    axs.set_ylabel('Shifted / Normal')
    # mplt.centre_spines(ax, 1, 1, arrow=False, centre_tick='x')
    mplt.clean_axes(axs)
    return fig


# =============================================================================
# ## Saving/loading
# =============================================================================


def save_data(fname: str, sss: Array, options: OptDict,
              envs: Arrays, models: Arrays) -> None:
    """Save generated data

    Parameters
    ----------
    fname : str
        Name of `.npz` file to save as, without extension.
    sss : Array (T,)
        Rate parameter of Laplace transform of SNR curve
    envs : Tuple[Array] (3,)(T,)
        Envelopes: maximum Laplace(SNR) at varios times.
    models : Tuple[Array] (3,)(T,Q)
        Parameters of the models that achieve the envelopes.
        `Q = 2M(M-1) or 2(M-1)`.
    Tuple members:-
        * Arbritrary topology, normal problem.
        * Serial topology, normal problem.
        * Arbritrary topology, shifted problem.
    """
    env_gen, env_srl, env_shf = envs
    models_gen, models_srl, models_shf = models
    nst = ma.params.num_state(models_srl.shape[-1], serial=True)
    np.savez_compressed(fname, nstate=nst, s=sss, options=options,
                        envelope_gen=env_gen, envelope_srl=env_srl,
                        envelope_shf=env_shf, models_gen=models_gen,
                        models_srl=models_srl, models_shf=models_shf)

def load_data(fname: str = 'optim') -> ty.Tuple[int, Array, OptDict,
                                                Arrays, Arrays]:
    """Load previously generated data

    Returns
    -------
    nst : int
        Number of states
    sss : Array (T,)
        Rate parameter of Laplace transform of SNR curve
    envs : Tuple[Array] (3,)(T,)
        Envelopes: maximum Laplace(SNR) at varios times.
    models : Tuple[Array] (3,)(T,Q)
        Parameters of the models that achieve the envelopes.
        `Q = 2M(M-1) or 2(M-1)`.
    Tuple members:-
        * Arbritrary topology, normal problem.
        * Serial topology, normal problem.
        * Arbritrary topology, shifted problem.
    """
    with la.load(fname + '.npz', allow_pickle=True) as saved:
        nst = saved['nstate'][()]
        sss = saved['s']
        options = saved['options'][()]
        envs = tuple(saved['envelope_' + x] for x in ['gen', 'srl', 'shf'])
        mods = tuple(saved['models_' + x] for x in ['gen', 'srl', 'shf'])
    return nst, sss, options, envs, mods


# =============================================================================
def main(save_npz: bool, save_pdf: bool, nsyn: float = 1):
    """Execute plot creation & saving
    """
    if save_npz:
        nst = 10
        sss = la.geomspace(1e-4, 10, 50)
        opts = {'repeats': 10, 'method': 'SLSQP'}
        # (gen, srl, shf), (gen, srl, shf)
        envs, mods = all_data(nst, sss, **opts)
        save_data('optim', sss, opts, envs, mods)
    else:
        nst, sss, opts, envs, mods = load_data('optim')
    jmp = 0.7

    fig_sc = ser_casc_plot(nst, jmp, sss, nsyn)
    fig_th = theory_plot(nst, sss, nsyn)
    fig_eq = equilibrium_plot(nst, sss, nsyn)
    fig_num = optim_plot(nst, sss, envs[0], envs[1], nsyn)
    fig_shift = shift_plot(sss, envs[0], envs[2], nsyn)
    fig_cond = cond_plot(sss, mods[2], envs[2], envs[0], nsyn)
    fig_heuristic = heuristic_plot(nst, sss, envs[1], mods[1], nsyn)

    if save_pdf:
        fig_sc.savefig("../../Notes/Figures/serial_vs_cascade.pdf")
        fig_th.savefig("../../Notes/Figures/LenvProven.pdf")
        fig_eq.savefig("../../Notes/Figures/LenvConj.pdf")
        fig_num.savefig("../../Notes/Figures/LenvNum.pdf")
        fig_shift.savefig("../../Notes/Figures/LenvShift.pdf")
        fig_cond.savefig("../../Notes/Figures/shift_cond.pdf")
        fig_heuristic.savefig("../../Notes/Figures/LenvHeuristic.pdf")


# =============================================================================
# Constants and Aliases
# =============================================================================
_FWD = {'linear': np.positive, 'log': np.log, 'symlog': np.log,
        'logit': lambda x: np.log(x/(1 - x))}
_BWD = {'linear': np.positive, 'log': np.exp, 'symlog': np.exp,
        'logit': lambda x: 1/(1 + np.exp(-x))}
_ALGN = {'left': 0, 'center': 0.5, 'centre': 0.5, 'right': 1,
         'bottom': 0, 'middle': 0.5, 'top': 1}
Array = la.lnarray
Arrays = ty.Tuple[Array, ...]
OptDict = ty.Dict[str, ty.Any]
Figure, Axes, Line = mpl.figure.Figure, mpl.axes.Axes, mpl.lines.Line2D

# =============================================================================
if __name__ == "__main__":
    main(False, False, 1e4)
