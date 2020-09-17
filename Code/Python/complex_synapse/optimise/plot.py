"""Creating figures for Envelopes
"""
import typing as ty

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import numpy_linalg as la
import sl_py_tools.matplotlib_tricks as mplt

import complex_synapse.optimise as cso
import complex_synapse.optimise.shorten as sh

np.set_printoptions(precision=3, suppress=False, linewidth=90)
mplt.rc_fonts()
mplt.rc_colours()

Array = la.lnarray
Arrays = ty.Tuple[Array, ...]
OptDict = ty.Dict[str, ty.Any]
# =============================================================================
# ## All data
# =============================================================================


def all_data(nst: int, sss: Array, **opts) -> ty.Tuple[Arrays, Arrays]:
    """All of the data generation
    """
    env_th = cso.proven_envelope_laplace(sss, nst)
    env_gen, mods_gen = cso.optim_laplace_range(sss, nst, serial=False, **opts)
    env_srl, mods_srl = cso.optim_laplace_range(sss, nst, serial=True, **opts)
    env_shf, mods_shf = cso.optim_laplace_range(sss, nst, cond_lim=True,
                                                maker=cso.shifted_problem,
                                                **opts)
    return (env_gen, env_srl, env_shf, env_th), (mods_gen, mods_srl, mods_shf)


# =============================================================================
# ## Theory
# =============================================================================


def theory_plot(nst: int, sss: Array, env_th: Array) -> mpl.figure.Figure:
    """Theoretical envelope plot
    """
    fig, axs = plt.subplots()
    t_pts = np.array([1e-1, 0.9e1, 1.1e1, 1e4])
    axs.loglog(1/sss, env_th * sss, label='theory')
    axs.loglog(t_pts[:2], np.ones(2), 'k:')
    axs.loglog(t_pts[-2:], (nst-1) / t_pts[-2:], 'k:')
    axs.set_xlabel(r"Time, $r\tau$")
    axs.set_ylabel("SNR")
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
               envelope_srl: Array):
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


def shift_plot(sss: Array, envelope_gen: Array, senvelope_gen: Array):
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


def heuristic_plot(nst: int, sss: Array, envelope_theory: Array,
                   envelope_gen: Array):
    """Heuristic envelope plot
    """
    s_two = sh.s_star(2)
    s_sticky = sh.s_star(nst)
    env_heuristic = cso.heuristic_envelope_laplace(sss, nst)

    fig, axs = plt.subplots()
    axs.loglog(1/sss, envelope_theory * sss, label='theory')
    axs.loglog(1/sss, envelope_gen * sss, label='numeric')
    axs.loglog(1/sss, env_heuristic * sss, label='heuristic')
    ylim = axs.set_ylim(2e-3, 1.1)
    axs.set_xlim(1/sss[-1], 1e3)
    axs.axvline(x=1/s_two, color='k', linestyle=':')
    axs.axvline(x=1/s_sticky, color='k', linestyle=':')
    axs.set_xlabel(r"Time, $r\tau$")
    axs.set_ylabel("SNR")
    ind = 27
    axs.text(1/sss[ind], env_heuristic[ind] * sss[ind],
             r"$\frac{\sqrt{N}(0.54)}{\sqrt{r\tau}}$",
             horizontalalignment='right', verticalalignment='top', fontsize=24)
    ind = 11
    axs.text(1/sss[ind], env_heuristic[ind] * sss[ind],
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


def reoptim(senvelope_gen: Array, envelope_gen: Array, smodels_gen: Array,
            sss: Array, **options):
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
              envelope_gen: Array):
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


def save_data(fname: str, nst: int, sss: Array, options: OptDict,
              envs: Arrays, models: Arrays):
    """Save generated data
    """
    env_gen, env_srl, env_shf, env_theory = envs
    models_gen, models_srl, models_shf = models
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
        save_data('optim', nst, sss, opts, envs, mods)

    fig_th = theory_plot(nst, sss, envs[-1])
    fig_num = optim_plot(sss, envs[-1], envs[0], envs[1])
    fig_shift = shift_plot(sss, envs[0], envs[2])
    fig_cond = cond_plot(sss, mods[2], envs[2], envs[0])
    fig_heuristic = heuristic_plot(nst, sss, envs[-1], envs[0])

    if save_pdf:
        fig_th.savefig("../../Notes/Figures/LenvProven.pdf")
        fig_num.savefig("../../Notes/Figures/LenvNum.pdf")
        fig_shift.savefig("../../Notes/Figures/LenvShift.pdf")
        fig_cond.savefig("../../Notes/Figures/shift_cond.pdf")
        fig_heuristic.savefig("../../Notes/Figures/LenvHeuristic.pdf")


# =============================================================================
if __name__ == "__main__":
    main(False, False)
