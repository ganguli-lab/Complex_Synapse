# -*- coding: utf-8 -*-
"""Creating video frames for optimal serial models
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import matplotlib as mpl
import matplotlib.animation as mpla
import matplotlib.pyplot as plt
import numpy as np

import numpy_linalg as la
import sl_py_tools.arg_tricks as ag
# import sl_py_tools.containers as cn
import sl_py_tools.iter_tricks as it
import sl_py_tools.matplotlib_tricks as mpt

from . import synapse_opt as so
from .. import synapse_mem as sm

mpt.rc_colours()
mpt.rc_fonts('sans-serif')
Figure = mpl.figure.Figure
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
    pot, dep = np.split(param, 2)  # pylint: disable=unbalanced-tuple-unpacking
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
    # pylint: disable=unbalanced-tuple-unpacking
    pot, dep = np.split(param.copy(), 2)
    nstw = (len(pot) + 1) // 2
    test = pot[nstw:] < 1e-3
    if test.any():
        ind = test.nonzero()[0][0] + nstw
        pot[ind:] = 1e-5
        dep[ind:] = 1e-5
    test = dep[:nstw-1] < 1e-3
    if test.any():
        ind = test.nonzero()[0][-1] + 1
        pot[:ind] = 1e-5
        dep[:ind] = 1e-5
    return pot, dep


# =============================================================================
# Fitter video figure
# =============================================================================


def _make_fig(ratios: Optional[Ratios] = None, scale: Optional[float] = None
              ) -> Tuple[Figure, AxList]:
    """Create a figure with axes for an experiment/simulation fit.

    Returns
    -------
    fig : Figure
        Figure object
    axs : List[Axes]
        Axes for [SNR, Equilibrium, Colourbar, Potentiation, Depression]
    """
    ratios = ag.default(ratios, ([4, 1, 4], [12, 9, 1]))
    scale = ag.default(scale, 0.6)
    keys = ['height_ratios', 'width_ratios']

    fsiz = tuple(scale * sum(g) for g in reversed(ratios))
    gsiz, kwargs = tuple(map(len, ratios)), dict(zip(keys, ratios))

    fig = plt.figure(figsize=fsiz, frameon=True, constrained_layout=True)
    gsp = fig.add_gridspec(*gsiz, **kwargs)
    fig.set_constrained_layout_pads(hpad=0, hspace=0)
    ax_snr = fig.add_subplot(gsp[:, 0])
    ax_eqp = fig.add_subplot(gsp[1, 1])
    ax_pot = fig.add_subplot(gsp[0, 1], sharex=ax_eqp)
    ax_dep = fig.add_subplot(gsp[2, 1], sharex=ax_eqp)
    ax_cbr = fig.add_subplot(gsp[:, 2])

    return fig, [ax_snr, ax_eqp, ax_cbr, ax_pot, ax_dep]


def _format_axes_snr(ax_snr: mpl.axes.Axes):
    """Format axes for SNR plots

    Parameters
    ----------
    ax_snr : Axes
        Axes for SNR
    """
    ax_snr.set_xlabel(r"Time, $\tau$")
    ax_snr.set_ylabel(r"SNR")
    ax_snr.legend(loc='lower left', edgecolor=ax_snr.get_facecolor())
    mpt.clean_axes(ax_snr, legendbox=False, tight=False)


def _format_axes_eqp(axs: AxList, imh: mpl.image.AxesImage, nst: int):
    """Format axes for equilibrium distribution plots

    Parameters
    ----------
    axs : List[Axes]
        Axes for [Equilibrium, Colourbar]
    nst : int
        Number of states
    """
    ax_eqp, ax_cbr = axs
    fig = ax_cbr.get_figure()
    fig.colorbar(imh, cax=ax_cbr)
    ax_eqp.set_xticks(np.arange(nst))
    ax_eqp.set_xticklabels([''] * nst)
    ax_eqp.set_yticks([])
    ax_eqp.xaxis.set_ticks_position('both')
    mpt.clean_axes(ax_eqp, box=False, tight=False)

    ax_cbr.set_ylabel("Equilibrium probability")
    mpt.clean_axes(ax_cbr, box=False, tight=False)


def _format_axes_bar(ax_bar: mpl.axes.Axes, nst: int, dep: bool):
    """Format axes for potentiation probability plots

    Parameters
    ----------
    axs : Axes
        Axes for Potentiation
    nst : int
        Number of states
    """
    ax_bar.set_ylim([0, 1])
    if dep:
        ax_bar.invert_yaxis()
        ax_bar.xaxis.set_ticks_position('top')
    ax_bar.set_ylabel("Trans. prob.")
    ax_bar.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax_bar.set_yticklabels(['0', '', '', '', '1'])
    ax_bar.set_xticks(np.arange(nst))
    ax_bar.set_xticklabels([''])
    mpt.clean_axes(ax_bar, tight=False)
    if dep:
        ax_bar.spines['top'].set_visible(True)
        ax_bar.spines['bottom'].set_visible(False)


# =============================================================================
# Plot
# =============================================================================


def make_plots(rate: la.lnarray, envs: Sequence[la.lnarray], jmps: la.lnarray,
               tau: float) -> Tuple[Figure, ModelPlots]:
    """Create the original plots
    """
    fig, axs = _make_fig()

    axs[0].loglog(1/rate, rate * envs[1], label="Theory")
    axs[0].loglog(1/rate, rate * envs[0], label="Numeric")

    obj = make_model_plots(axs, 1/rate, jmps, tau)
    nst = len(obj.brp) + 1

    _format_axes_snr(axs[0])
    _format_axes_eqp(axs[1:3], obj.imh, nst)
    _format_axes_bar(axs[3], nst, False)
    _format_axes_bar(axs[4], nst, True)

    # fig.set_constrained_layout(True)
    fig.set_constrained_layout_pads(hpad=0, hspace=0)
    plt.draw()
    # fig.set_constrained_layout(False)
    return fig, obj


def make_model_plots(axs: AxList, time: la.lnarray,
                     model: la.lnarray, tau: float) -> ModelPlots:
    """Create the original plots for a model
    """
    snr = serial_snr(model, time)
    eqp = serial_eqp(model)
    pot, dep = trim_params(model)
    njmp = len(pot)

    lnh = axs[0].loglog(time, snr, label="One model")
    vln = axs[0].axvline(tau)

    # imh = axs[1].imshow(eqp.r, norm=mpl.colors.Normalize(0, 1))
    imh = axs[1].pcolormesh(np.arange(-0.5, njmp+1), [0, 1], eqp.r,
                            norm=mpl.colors.Normalize(0, 1),
                            edgecolors='black')

    brp = axs[3].bar(np.arange(0.5, njmp), pot)
    brd = axs[4].bar(np.arange(0.5, njmp), dep)

    return ModelPlots([lnh[0], vln], imh, brp, brd)


class ModelPlots:
    """The plots associated with a model
    """
    lnh: mpl.lines.Line2D
    vln: mpl.lines.Line2D
    # imh: mpl.image.AxesImage
    imh: mpl.collections.QuadMesh
    brp: mpl.container.BarContainer
    brd: mpl.container.BarContainer

    def __init__(self, lns: Sequence[mpl.lines.Line2D],
                 imh: mpl.image.AxesImage,
                 brp: mpl.container.BarContainer,
                 brd: mpl.container.BarContainer) -> None:
        self.lnh = lns[0]
        self.vln = lns[1]
        self.imh = imh
        self.brp = brp
        self.brd = brd

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
        """
        snr = serial_snr(model, time)
        eqp = serial_eqp(model)
        pot, dep = trim_params(model)

        self.update_snr(snr, tau)
        self.update_eqp(eqp)
        self.update_potdep(pot, dep)

        plt.draw()

    @classmethod
    def from_data(cls, axs: AxList, time: la.lnarray, model: la.lnarray,
                  tau: float) -> ModelPlots:
        """Plot data to make an instance"""
        snr = serial_snr(model, time)
        eqp = serial_eqp(model)
        pot, dep = trim_params(model)
        njmp = len(pot)

        lnh = axs[0].loglog(time, snr, label="One model")
        vln = axs[0].axvline(tau)

        # imh = axs[1].imshow(eqp.r, norm=mpl.colors.Normalize(0, 1))
        imh = axs[1].pcolormesh(np.arange(-0.5, njmp+1), [0, 1], eqp.r,
                                norm=mpl.colors.Normalize(0, 1),
                                edgecolors='black')

        brp = axs[3].bar(np.arange(0.5, njmp), pot)
        brd = axs[4].bar(np.arange(0.5, njmp), dep)

        return cls([lnh[0], vln], imh, brp, brd)


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

    def __init__(self, rate: la.lnarray, env_th: la.lnarray,
                 env_num: la.lnarray, models: la.lnarray) -> None:
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
        models : la.lnarray (T,2N-2)
            Parameters of serial models that form envelope, in the order
            mat_01, mat_12, ..., mat_n-2,n-1,
            mat_10, mat_21, ..., mat_n-1,n-2.
        """
        self.time = 1 / rate
        self.env_th = env_th * rate
        self.env_num = env_num * rate
        self.models = models
        if self.time[0] > self.time[-1]:
            self.flip()

        self.fig, axs = _make_fig()

        axs[0].loglog(self.time, self.env_th, label="Theory")
        axs[0].loglog(self.time, self.env_num, label="Numeric")

        self.model_plots = ModelPlots.from_data(axs, self.time, self.models[0],
                                                self.time[0])

        _format_axes_snr(axs[0])
        _format_axes_eqp(axs[1:3], self.model_plots.imh, self.nstate)
        _format_axes_bar(axs[3], self.nstate, False)
        _format_axes_bar(axs[4], self.nstate, True)

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
        """Change which model we're plotting
        """
        self.model_plots.update_plots(self.time, self.models[ind],
                                      self.time[ind])

    def savefig(self, fname: str):
        """Save the current figure
        """
        self.fig.savefig(fname)

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
    fitter : SynapseFitter
        The synapse fitter we will animate. An instance of `FitterReplay` or
        one of its subclasses is recommended.
    Other keywords
        Passed to `FuncAnimation`.

    Returns
    -------
    ani : FuncAnimation
        The animation.
    """
    kwargs.setdefault('init_func', None)
    kwargs.setdefault('frames', it.erange(env_fig.ntime))
    kwargs.setdefault('interval', 500)
    kwargs.setdefault('repeat_delay', 1000)
    kwargs.setdefault('blit', False)
    return mpla.FuncAnimation(env_fig.fig, env_fig.update, **kwargs)


# =============================================================================
# Hint aliases
# =============================================================================
AxList = List[mpl.axes.Axes]
Ratios = Tuple[List[int], List[int]]
