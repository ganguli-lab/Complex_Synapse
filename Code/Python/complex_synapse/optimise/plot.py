"""Creating figures for Envelopes
"""
import typing as ty
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy_linalg as la
import sl_py_tools.matplotlib_tricks as mplt
import complex_synapse.optimise as cso
import complex_synapse.optimise.shorten as sh

np.set_printoptions(precision=3, suppress=False, linewidth=90)
mplt.rc_fonts()
mplt.rc_colours()


# =============================================================================
# ## Theory
# =============================================================================


def theory_plot(sss, envelope_theory) -> mpl.figure.Figure:
    """Theoretical envelope plot
    """
    fig, axs = plt.subplots()
    t_pts = np.array([1e-1, 0.9e1, 1.1e1, 1e4])
    axs.loglog(1/sss, envelope_theory * sss, label='theory')
    axs.loglog(t_pts[:2], np.ones(2), 'k:')
    axs.loglog(t_pts[-2:], (M-1) / t_pts[-2:], 'k:')
    axs.set_xlabel(r"Time, $r\tau$")
    axs.set_ylabel("SNR")
    axs.set_xlim(1e-1, 1e3)
    axs.set_ylim(1e-2, 3)
    ind = 30
    axs.text(t_pts[2]/2, 1.1, r"$\sqrt{N}$", fontsize=18,
             horizontalalignment='right', verticalalignment='bottom')
    axs.text(t_pts[-2]*4, (M-1) / (t_pts[-2] * 4),
             r"$\frac{\sqrt{N}(M-1)}{r\tau}$", fontsize=20,
             horizontalalignment='left', verticalalignment='bottom')
    axs.text(1/sss[ind], envelope_theory[ind] * sss[ind],
             r"$\frac{\sqrt{N}(M-1)}{r\tau + (M-1)}$",
             horizontalalignment='right', verticalalignment='top', fontsize=24)
    mplt.clean_axes(axs)
    return fig


# =============================================================================
# ## Optimisation
# =============================================================================


def optim_data(nst, sss, **options):
    """Data for optimisation plot
    """
    envelope_gen, models_gen = cso.optim_laplace_range(
        sss, nst, serial=False, **options)
    envelope_srl, models_srl = cso.optim_laplace_range(
        sss, nst, serial=True, **options)
    return envelope_gen, envelope_srl, models_gen, models_srl


def optim_plot(sss, envelope_theory, envelope_gen, envelope_srl):
    """Optimisation plot
    """
    fig, axs = plt.subplots()
    axs.loglog(1/sss, envelope_theory * sss, label='theory')
    axs.loglog(1/sss, envelope_gen * sss, label='numeric: general')
    axs.loglog(1/sss, envelope_srl * sss, label='numeric: serial')
    axs.set_xlabel(r"Time, $r\tau$")
    axs.set_ylabel("SNR")
    axs.set_xlim(1/sss[-1], 1/sss[0])
    axs.legend(loc="lower left")
    mplt.clean_axes(axs)
    return fig


# =============================================================================
# ### Shifted optimisation
# =============================================================================


def shift_data(nst, sss, **options):
    """Data for shifted optimisation plot
    """
    senvelope_gen, smodels_gen = cso.optim_laplace_range(
        sss, nst, maker=cso.shifted_problem, cond_lim=True,
        **options)
    return senvelope_gen, smodels_gen


def shift_plot(sss, envelope_gen, senvelope_gen):
    """Shifted optimisation plot
    """
    fig, axs = plt.subplots()
    axs.loglog(1/sss, envelope_gen * sss, label='normal')
    axs.loglog(1/sss, senvelope_gen * sss, label='shifted')
    axs.set_xlabel(r"Time, $r\tau$")
    axs.set_ylabel("SNR")
    axs.set_xlim(1/sss[-1], 1/sss[0])
    axs.legend(loc="lower left")
    mplt.clean_axes(axs)
    return fig


# =============================================================================
# ## Heuristic
# =============================================================================


def heuristic_plot(nst, envelope_theory, envelope_gen):
    """Heuristic envelope plot
    """
    s_two = sh.s_star(2)
    s_sticky = sh.s_star(nst)
    env_heuristic = cso.heuristic_envelope_laplace(s, nst)

    fig, axs = plt.subplots()
    axs.loglog(1/s, envelope_theory * s, label='theory')
    axs.loglog(1/s, envelope_gen * s, label='numeric')
    axs.loglog(1/s, env_heuristic * s, label='heuristic')
    ylim = axs.set_ylim(2e-3, 1.1)
    axs.set_xlim(1/s[-1], 1e3)
    axs.axvline(x=1/s_two, color='k', linestyle=':')
    axs.axvline(x=1/s_sticky, color='k', linestyle=':')
    axs.set_xlabel(r"Time, $r\tau$")
    axs.set_ylabel("SNR")
    ind = 27
    axs.text(1/s[ind], env_heuristic[ind] * s[ind],
             r"$\frac{\sqrt{N}(0.54)}{\sqrt{r\tau}}$",
             horizontalalignment='right', verticalalignment='top', fontsize=24)
    ind = 11
    axs.text(1/s[ind], env_heuristic[ind] * s[ind],
             r"$\frac{\sqrt{N}(M-1)}{r\tau}$",
             horizontalalignment='right', verticalalignment='top', fontsize=24)
    axs.text(1/s_two, ylim[0], r"$r\tau = 0.73$", fontsize=16,
             horizontalalignment='center', verticalalignment='bottom')
    axs.text(1/s_sticky, ylim[0] * 1.1, r"$r\tau = 0.22M^2$", fontsize=16,
             horizontalalignment='center', verticalalignment='bottom')
    axs.legend(loc="upper right")
    mplt.clean_axes(axs)
    return fig


# =============================================================================
# ## Diagnosis and treatment
# =============================================================================


def reoptim(senvelope_gen, envelope_gen, smodels_gen, sss, **options):
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


def cond_plot(sss, smodels_gen, senvelope_gen, envelope_gen):
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


def save_data(fname, nst, sss, envs, models, **options):
    """Save generated data
    """
    envelope_theory, envelope_gen, envelope_srl, senvelope_gen = envs
    models_gen, models_srl, smodels_gen = models
    np.savez_compressed(fname, sss=sss, envelope_theory=envelope_theory,
                        envelope_gen=envelope_gen, envelope_srl=envelope_srl,
                        senvelope_gen=senvelope_gen,
                        models_gen=models_gen, models_srl=models_srl,
                        smodels_gen=smodels_gen,
                        nst=nst, options=options)

def load_data(fname: str = 'optim.npz') -> ty.Tuple[la.lnarray, ...]:
    """Load previously generated data
    """
    with np.load(fname, allow_pickle=True) as saved:
        sss = saved['sss']
        envelope_theory = saved['envelope_theory']
        envelope_gen = saved['envelope_gen']
        envelope_srl = saved['envelope_srl']
        senvelope_gen = saved['senvelope_gen']
        models_gen = saved['models_gen']
        models_srl = saved['models_srl']
        smodels_gen = saved['smodels_gen']
        options = saved['options'][()]
        nst = saved['nst'][()]
    return (nst, sss, options, envelope_theory,
            envelope_gen, envelope_srl, senvelope_gen,
            models_gen, models_srl, smodels_gen)


# =============================================================================
if __name__ == "__main__":
    M = 10
    s = la.geomspace(1e-4, 10, 50)
    opts = {'repeats': 10, 'method': 'SLSQP'}
    env_th = cso.proven_envelope_laplace(s, M)
    env_gen, env_srl, mods_gen, mods_srl = optim_data(M, s, **opts)
    senv_gen, smods_gen = shift_data(M, s, cond_thresh=1e3, **opts)

    fig_th = theory_plot(s, env_th)
    fig_th.savefig("../../Notes/Figures/LenvProven.pdf")

    fig_num = optim_plot(s, env_th, env_gen, env_srl)
    fig_num.savefig("../../Notes/Figures/LenvNum.pdf")

    fig_shift = shift_plot(s, env_gen, senv_gen)
    fig_shift.savefig("../../Notes/Figures/LenvShift.pdf")

    fig_cond = cond_plot(s, smods_gen, senv_gen, env_gen)
    fig_cond.savefig("../../Notes/Figures/shift_cond.pdf")

    fig_heuristic = heuristic_plot(M, env_th, env_gen)
    fig_heuristic.savefig("../../Notes/Figures/LenvHeuristic.pdf")

    save_data('optim', M, s, (env_th, env_gen, env_srl, senv_gen),
              (mods_gen, mods_srl, smods_gen), **opts)
