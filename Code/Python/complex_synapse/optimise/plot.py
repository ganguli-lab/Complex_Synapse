"""Creating figures for Envelopes
"""
from __future__ import annotations
import typing as ty

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
    """
    env_th = cso.proven_envelope_laplace(sss, nst)
    env_gen, mods_gen = cso.optim_laplace_range(sss, nst, **opts)
    with dt.updated(opts, serial=True) as opt:
        env_srl, mods_srl = cso.optim_laplace_range(sss, nst, **opt)
    with dt.updated(opts, cond_lim=True, maker=cso.shifted_problem) as opt:
        env_shf, mods_shf = cso.optim_laplace_range(sss, nst, **opt)
    return (env_gen, env_srl, env_shf, env_th), (mods_gen, mods_srl, mods_shf)


# =============================================================================
# ## Generic
# =============================================================================


def mem_plot(sss: Array, mems: ty.Dict[str, Array]) -> ty.Tuple[Figure, Axes]:
    """Plot memory curves

    Parameters
    ----------
    sss : Array
        1/time
    mems : Dict[str, Array]
        Laplace transforms of memory curves to plot, keyed by plot label

    Returns
    -------
    fig : mpl.figure.Figure
        Figure containing `axs`.
    axs : mpl.axes.Axes
        Axes containing plots.
    """
    fig, axs = plt.subplots()
    for label, mem in mems.items():
        axs.loglog(1/sss, mem * sss, label=label)
    axs.set_xlabel(r"Time, $r\tau$")
    axs.set_ylabel("SNR")
    return fig, axs



# =============================================================================
# ## Theory
# =============================================================================


def theory_plot(nst: int, sss: Array, env_th: Array) -> Figure:
    """Theoretical envelope plot
    """
    fig, axs = mem_plot(sss, {'theory': env_th})
    t_pts = np.array([1e-1, 0.9e1, 1.1e1, 1e4])
    axs.loglog(t_pts[:2], np.ones(2), 'k:')
    axs.loglog(t_pts[-2:], (nst-1) / t_pts[-2:], 'k:')
    axs.set_xlim(1e-1, 1e3)
    axs.set_ylim(1e-2, 3)
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


# =============================================================================
# ## Optimisation
# =============================================================================


def optim_plot(sss: Array, envelope_theory: Array, envelope_gen: Array,
               envelope_srl: Array) -> Figure:
    """Optimisation plot
    """
    fig, axs = mem_plot(sss, {'theory': envelope_theory,
                              'numeric: general': envelope_gen,
                              'numeric: serial': envelope_srl})
    axs.set_xlim(1/sss[-1], 1/sss[0])
    axs.legend(loc="lower left")
    mplt.clean_axes(axs)
    return fig


# =============================================================================
# ### Shifted optimisation
# =============================================================================


def shift_plot(sss: Array, envelope_gen: Array, senvelope_gen: Array) -> Figure:
    """Shifted optimisation plot
    """
    fig, axs = mem_plot(sss, {'normal': envelope_gen,
                              'shifted': senvelope_gen})
    axs.set_xlim(1/sss[-1], 1/sss[0])
    axs.legend(loc="lower left")
    mplt.clean_axes(axs)
    return fig


# =============================================================================
# ## Heuristic
# =============================================================================


def heuristic_plot(nst: int, sss: Array, envelope_theory: Array,
                   envelope_gen: Array, models_srl: Array, **kwds) -> Figure:
    """Heuristic envelope plot
    """
    kwds.setdefault('model_inds', (43, 25, 9))
    kwds.setdefault('model_above', (0, 0, 0))
    kwds.setdefault('model_siz', (2.3, 1.7))
    kwds.setdefault('note_inds', (43, 25, 9))
    kwds.setdefault('note_above', (1, 1, 1))
    kwds.setdefault('legend_loc', "upper right")
    kwds.setdefault('legendfontscale', 0.8)
    kwds.setdefault('legendbox', False)
    ax_opts = mplt.clean_axes_keys(kwds)

    env_heuristic = cso.heuristic_envelope_laplace(sss, nst)

    fig, axs = mem_plot(sss, {'theory': envelope_theory,
                              'numeric': envelope_gen,
                              'heuristic': env_heuristic})
    axs.set_xlim(1/sss[-1], 1e4)
    ylim = axs.set_ylim(2e-3, 1.1)

    _annotate_curves(axs, 1/sss, env_heuristic * sss, models_srl, **kwds)
    _annotate_axis(axs, nst, ylim[0])

    axs.legend(loc=kwds['legend_loc'], edgecolor=axs.get_facecolor())
    mplt.clean_axes(axs, **ax_opts)
    return fig


def _annotate_curves(axs, time, env, models, **kwds) -> None:
    """Add text and graphs along envelope
    """
    args = (('right', 'top'), ('left', 'middle'))
    for ind, above in zip(kwds['model_inds'], kwds['model_above']):
        add_graph(axs, [time[ind-1], env[ind], *kwds['model_siz']],
                  models[ind], *args[above])

    txt_opt = {'size': 24}
    txt_opt = ({'ha': 'right', 'va': 'top', **txt_opt},
               {'ha': 'left', 'va': 'bottom', **txt_opt})
    notes = (r"$\sqrt{N}$", r"$\frac{\sqrt{N}(0.54)}{\sqrt{r\tau}}$",
             r"$\frac{\sqrt{N}(M-1)}{r\tau}$")
    for ind, above, note in zip(kwds['note_inds'], kwds['note_above'], notes):
        axs.text(time[ind], env[ind-1], note, **txt_opt[above])


def _annotate_axis(axs, nst, y_0) -> None:
    """Add text and graphs along time axis
    """
    t_two = 1 / sh.s_star(2)
    t_sticky = 1 / sh.s_star(nst)

    axs.axvline(x=t_two, color='k', linestyle=':')
    axs.axvline(x=t_sticky, color='k', linestyle=':')

    txt_opt = {'ha': 'center', 'va': 'bottom', 'size': 16}
    axs.text(t_two, y_0, r"$r\tau = 0.73$", **txt_opt)
    axs.text(t_sticky, y_0 * 1.1, r"$r\tau = 0.22M^2$", **txt_opt)


def add_graph(axs: mpl.axes.Axes, bounds: ty.Sequence[float], model: Array,
              hal: str = 'left', val: str = 'bottom',
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
    model = trim_params(model)
    nst = ma.params.num_state(model.shape[-1], serial=True)
    bounds, scale = scale_bounds(bounds, nst, (3, 2), axs)
    bounds = logify(bounds, axs, hal, val)
    kwds.setdefault('edge_style.mut_scale', 0.1 * scale)
    kwds.setdefault('edge_style.mult', 2 * scale)
    kwds.setdefault('node_style.mult', 150 * scale**2)
    opts = gp.GraphOptions(kwds)

    axin = axs.inset_axes(bounds, transform=axs.transData, frame_on=False)
    graph = GraphPlots.from_data(axin, model, opts=opts)
    axin.set_xlim(-0.5, nst + 0.5)
    axin.set_ylim(-0.75, 0.75)
    axin.set_clip_on(False)
    return axin, graph


def logify(bounds: ty.Sequence[float], axs: mpl.axes.Axes,
           hal: str = 'left', val: str = 'bottom') -> ty.List[float]:
    """Adjust bounds for logarithmic axes

    Parameters
    ----------
    bounds : Sequence[float]
        `[left, bottom, width, height]`, where `width` and `height` are
        multiplicative when the corresponding axis is logarithmic.
    axs : mpl.axes.Axes
        The axes to which the bounds refer.

    Returns
    -------
    bounds : List[float]
        New bounds, with positive additive `width` and `height`
    """
    x_0, x_1 = _logify(bounds[::2], axs.get_xscale(), hal)
    y_0, y_1 = _logify(bounds[1::2], axs.get_yscale(), val)
    return [x_0, y_0, x_1 - x_0, y_1 - y_0]


_FWD = {'linear': np.positive, 'log': np.log, 'symlog': np.log,
        'logit': lambda x: np.log(x/(1 - x))}
_BWD = {'linear': np.positive, 'log': np.exp, 'symlog': np.exp,
        'logit': lambda x: 1/(1 + np.exp(-x))}
_ALGN = {'left': 0, 'center': 0.5, 'centre': 0.5, 'right': 1,
         'bottom': 0, 'middle': 0.5, 'top': 1}

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


def reoptim(senvelope_gen: Array, envelope_gen: Array, smodels_gen: Array,
            sss: Array, **options) -> Arrays:
    """Reattempt failed optimisations
    """
    bad_inds, = (senvelope_gen / envelope_gen < 0.6).nonzero()
    reoptions = options.copy()
    reoptions['repeats'] = 20
    senvelope_gen, smodels_gen = cso.reoptim_laplace_range(
        bad_inds, sss, senvelope_gen, smodels_gen,
        maker=cso.shifted_problem, cond_thresh=1e3, cond_lim=False, **reoptions
    )
    return senvelope_gen, smodels_gen


def cond_plot(sss: Array, smodels_gen: Array, senvelope_gen: Array,
              envelope_gen: Array) -> mpl.figure.Figure:
    """Plot discrepancy vs condition number
    """
    conds = cso.check_cond_range(sss, smodels_gen)

    fig, axs = plt.subplots()
    # good = senvelope_gen < envelope_gen
    # bad = np.logical_not(good)
    # ax.semilogx(conds[good], senvelope_gen[good] / envelope_gen[good], '+',
    #             conds[bad], senvelope_gen[bad] / envelope_gen[bad], '+')
    axs.semilogx(conds, senvelope_gen / envelope_gen, '+')
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
              envs: Arrays, models: Arrays):
    """Save generated data
    """
    env_gen, env_srl, env_shf, env_theory = envs
    models_gen, models_srl, models_shf = models
    nst = ma.params.num_state(models_srl.shape[-1], serial=True)
    np.savez_compressed(fname, nstate=nst, s=sss, options=options,
                        envelope_gen=env_gen, envelope_srl=env_srl,
                        envelope_shf=env_shf, envelope_theory=env_theory,
                        models_gen=models_gen, models_srl=models_srl,
                        models_shf=models_shf)

def load_data(fname: str = 'optim') -> ty.Tuple[int, Array, OptDict,
                                                Arrays, Arrays]:
    """Load previously generated data
    """
    with la.load(fname + '.npz', allow_pickle=True) as saved:
        nst = saved['nstate'][()]
        sss = saved['s']
        options = saved['options'][()]
        envs = tuple(saved[x] for x in ['envelope_gen', 'envelope_srl',
                                        'envelope_shf', 'envelope_theory'])
        mods = tuple(saved[x] for x in ['models_gen', 'models_srl',
                                        'models_shf'])
    return nst, sss, options, envs, mods


# =============================================================================
def main(save_npz: bool, save_pdf: bool):
    """Execute polt creation & saving
    """
    nst = 10
    sss = la.geomspace(1e-4, 10, 50)
    opts = {'repeats': 10, 'method': 'SLSQP'}
    # (gen, srl, shf, th), (gen, srl, shf)
    envs, mods = all_data(nst, sss, **opts)

    if save_npz:
        save_data('optim', sss, opts, envs, mods)

    fig_th = theory_plot(nst, sss, envs[-1])
    fig_num = optim_plot(sss, envs[-1], envs[0], envs[1])
    fig_shift = shift_plot(sss, envs[0], envs[2])
    fig_cond = cond_plot(sss, mods[2], envs[2], envs[0])
    fig_heuristic = heuristic_plot(nst, sss, envs[-1], envs[0], mods[1])

    if save_pdf:
        fig_th.savefig("../../Notes/Figures/LenvProven.pdf")
        fig_num.savefig("../../Notes/Figures/LenvNum.pdf")
        fig_shift.savefig("../../Notes/Figures/LenvShift.pdf")
        fig_cond.savefig("../../Notes/Figures/shift_cond.pdf")
        fig_heuristic.savefig("../../Notes/Figures/LenvHeuristic.pdf")


# =============================================================================
# Aliases
# =============================================================================
Array = la.lnarray
Arrays = ty.Tuple[Array, ...]
OptDict = ty.Dict[str, ty.Any]
Figure, Axes = mpl.figure.Figure, mpl.axes.Axes

# =============================================================================
if __name__ == "__main__":
    main(False, False)
