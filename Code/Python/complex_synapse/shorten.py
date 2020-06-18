# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# -------------------------------------
"""Properties of optimal shortened models
"""
from typing import Tuple
import numpy as np
import scipy.optimize as sco

Arrays = Tuple[np.ndarray, ...]


# -------------------------------------
sol = sco.root_scalar(lambda x: np.sinh(x) - np.tanh(x/2) - x, bracket=(1, 2))
Y_STAR = sol.root
A_STAR = 2 * np.sinh(Y_STAR / 2)**2 / (Y_STAR * np.cosh(Y_STAR))


# -------------------------------------
def alpha_to_s(alpha: float) -> float:
    """Convert alpha to s"""
    return (alpha - 1)**2 / (2 * alpha)


def s_to_alpha(sval: float) -> (float, float):
    """Convert s to alpha two possibilities, reciprocals"""
    sval = sval + 1
    disc = np.sqrt(sval**2 - 1)
    return sval - disc, sval + disc


# -------------------------------------
def uniform(alpha: float, num: int) -> float:
    """Laplace-SNR for uniform serial model"""
    afm = alpha**(num / 2)
    sval = alpha_to_s(alpha)
    return 2 * (afm - 1)**2 / (num * sval * (afm**2 + 1))


# -------------------------------------
def num_star(alpha: float) -> float:
    """Optimal num for given alpha"""
    return 2 * Y_STAR / np.abs(np.log(alpha))


def alpha_star(num: int) -> (float, float):
    """Alpha where num is optimal"""
    alpha = np.exp(2 * Y_STAR / num)
    return 1 / alpha, alpha


def uni_star(alpha: float) -> float:
    """Heuristic envelope"""
    return A_STAR * np.abs(np.log(alpha)) / alpha_to_s(alpha)


# -------------------------------------
def short(eps: float, alpha: float, num: int) -> float:
    """Laplace-SNR for shortened serial model"""
    afm = alpha**(num / 2)
    afmm = alpha**((num-2) / 2)
    sval = alpha_to_s(alpha)
    fac = alpha**2 - alpha + 1
    numer = (afm - 1)**2
    denom = (afm**2 + 1)
    dnumer = numer - fac * (afmm - 1)**2
    ddenom = denom - fac * (afmm**2 + 1)
    fnumer = numer - eps * dnumer
    fdenom = denom - eps * ddenom
    return 2 * fnumer / ((num - 2*eps) * sval * fdenom)


def neg_short(eps: float, alpha: float, num: int) -> float:
    """negative Laplace-SNR for shortened serial model"""
    return - short(eps, alpha, num)


# -------------------------------------
def eps_star(alpha: float, num: int) -> float:
    """Laplace-SNR for shortened serial model"""
    if num == 2:
        return 0, 0
    afm = alpha**(num / 2)
    afmm = alpha**((num-2) / 2)
    fac = alpha**2 - alpha + 1
    numer = (afm - 1)**2
    denom = (afm**2 + 1)
    dnumer = numer - fac * (afmm - 1)**2
    ddenom = denom - fac * (afmm**2 + 1)
    fnumer = numer / dnumer
    fdenom = denom / ddenom
    disc = np.sqrt((num/2 - fnumer) * (fdenom - fnumer))
    return np.minimum(np.maximum(fnumer - disc, 0), 1)


# -------------------------------------
def lower(alpha: float, num: int) -> float:
    """derivative of short wrt eps at eps==1
    Note: lower t -> higher s -> (lower, higher) alpha[:],
    """
    afm = alpha**(num / 2)
    afmm = alpha**((num-2) / 2)
    fac = alpha**2 - alpha + 1
    numer = (afmm - 1)**2
    denom = (afmm**2 + 1)
    dnumer = (afm - 1)**2  # - fac * numer cancels
    ddenom = (afm**2 + 1)  # - fac * denom cancels
    return (num - 2)*(denom*dnumer - numer*ddenom) - 2*fac*numer*denom


# -------------------------------------
def upper(alpha: float, num: int) -> float:
    """derivative of short wrt eps at eps==0
    Note: higher t -> lower s -> (higher, lower) alpha[:],
    """
    afm = alpha**(num / 2)
    afmm = alpha**((num-2) / 2)
    fac = alpha**2 - alpha + 1
    numer = (afm - 1)**2
    denom = (afm**2 + 1)
    dnumer = - (afmm - 1)**2  # + numer cancels
    ddenom = - (afmm**2 + 1)  # + denom cancels
    return num*fac*(denom*dnumer - numer*ddenom) - 2*numer*denom


# -------------------------------------
def limits(num: int, debug: bool = False) -> float:
    """range of alpha where shortened model can be optimised
    Note: higher t -> lower s -> higher M
    """
    if num == 2:
        return (0.99, 0.01)
    x_0, x_1 = alpha_star(num - 2)[0], alpha_star(num)[0]
    lo_sol = sco.root_scalar(lower, args=(num,), x0=x_0, x1=x_1)
    if debug:
        print(lo_sol.flag, eps_star(lo_sol.root, num))
    x_0, x_1 = alpha_star(num)[0], alpha_star(num + 2)[0]
    hi_sol = sco.root_scalar(upper, args=(num,), x0=x_0, x1=x_1)
    if debug:
        print(hi_sol.flag, eps_star(hi_sol.root, num))
    return lo_sol.root, hi_sol.root


# -------------------------------------
def envelope(num: int, count: int, **kwds) -> Arrays:
    """Optimised shortened model"""
    lims = kwds.pop('lims', limits(num))
    fudge = kwds.pop('fudge', 0.0)
    lims = ((1-fudge) * lims[0], (1+fudge) * lims[1])
    alphas = np.geomspace(*lims, count)
    svals = alpha_to_s(alphas)
    env = uni_star(alphas)
    if num == 2:
        return svals, 1 / (1 + svals), env
    epss = eps_star(alphas, num)
    avals = short(epss, alphas, num)
    return svals, avals, env


# -------------------------------------
def envelopes(num_max: int, count: int) -> Arrays:
    """set of optimised shortened models"""
    siz = (num_max // 2, count)
    svals, avals, env = np.empty(siz), np.empty(siz), np.empty(siz)
    for i in range(num_max // 2):
        svals[i], avals[i], env[i] = envelope(2*i + 2, count)
    return svals, avals, env


# -------------------------------------
def env_ravel(svals: np.ndarray, env: np.ndarray) -> Arrays:
    """ravel and sort env"""
    new_s, new_env = svals.ravel(), env.ravel()
    inds = np.argsort(new_s)
    return new_s[inds], new_env[inds]
