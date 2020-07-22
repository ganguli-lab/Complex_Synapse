"""Properties of optimal shortened models
"""
from typing import Tuple

import numpy as np
import scipy.optimize as sco

from .shorten import Y_STAR

Arrays = Tuple[np.ndarray, ...]


# -------------------------------------
def beta_to_s(beta: float) -> float:
    """Convert alpha to s"""
    return np.cosh(beta) - 1


def s_to_beta(sval: float) -> float:
    """Convert s to beta. positive one of two possibilities"""
    return np.arccosh(sval + 1)


# -------------------------------------
def beta_star(num: int) -> float:
    """Value of beta at which num is optimal"""
    return 2 * Y_STAR / num


# =============================================================================


def components(beta: float, num: int) -> Tuple[float, ...]:
    """Components of numerator and denominator

    numerator = numer - eps * dnumer
    denominator = denom - eps * ddenom
    """
    denom = np.cosh(num * beta / 2)
    numer = denom - 1
    ddenom = np.cosh((num - 2) * beta / 2)
    dnumer = ddenom - 1
    return numer, denom, dnumer, ddenom


# -------------------------------------
def sticky(eps: float, beta: float, num: int) -> float:
    """Laplace-SNR for sticky serial model"""
    sval = beta_to_s(beta)
    numer, denom, dnumer, ddenom = components(beta, num)
    fnumer = numer - eps * dnumer
    fdenom = denom - eps * ddenom
    fmu = num - eps * (num - 2)
    return 2 * (1 - eps) * fnumer / (fmu * sval * fdenom)


def neg_sticky(eps: float, alpha: float, num: int) -> float:
    """negative Laplace-SNR for sticky serial model"""
    return - sticky(eps, np.log(alpha), num)


# -------------------------------------
def eps_stars(beta: float, num: int) -> float:
    """Optimal epsilon for sticky serial model, raw"""
    if num == 2:
        return 0, 0
    numer, denom, dnumer, ddenom = components(beta, num)
    fnumer = numer / dnumer
    fdenom = denom / ddenom
    fmu = num / (num - 2)
    if np.isscalar(fnumer):
        fnumer = fmu**2 if np.isclose(beta, 0) else fnumer
    else:
        fnumer[np.isclose(np.array(beta), 0)] = fmu**2
    disc = np.sqrt((fmu - 1)*(fdenom - 1)*(fmu - fnumer)*(fdenom - fnumer))
    quad = 1 + fnumer - fdenom - fmu
    linr = (fnumer - fdenom * fmu)
    return (linr - disc) / quad, (linr + disc) / quad


def eps_star(beta: float, num: int) -> float:
    """Optimal epsilon for sticky serial model, clipped"""
    epss = eps_stars(beta, num)[0]
    return np.minimum(np.maximum(epss, 0), 1)


def eps_star_star(beta: float, num: int) -> float:
    """Optimal epsilon for sticky serial model, clipped, other soln"""
    epss = eps_stars(beta, num)[1]
    return np.minimum(np.maximum(epss, 0), 1)


# -------------------------------------
def lower(beta: float, num: int) -> float:
    """derivative of sticky wrt eps at eps==0
    """
    numer, denom, dnumer, ddenom = components(beta, num)
    # d/d(eps) (1-eps)(n - eps dn)/(m - (m-2) eps)(d - eps dd)
    return num * (denom*dnumer - numer*ddenom) + 2 * numer * denom


# -------------------------------------
def limits(num: int, debug: bool = False) -> float:
    """range of beta where sticky model can be optimised
    Note: higher t -> lower s -> higher M
    """
    if num == 2:
        return (0.01, np.exp(-Y_STAR))
    x_0, x_1 = beta_star(num), beta_star(num + 2)
    lo_sol = sco.root_scalar(lower, args=(num,), x0=x_0, x1=x_1)
    if debug:
        print(lo_sol.flag, eps_star(lo_sol.root, num))
    return lo_sol.root, 0


# -------------------------------------
def envelope(num: int, count: int, **kwds) -> Arrays:
    """Optimised sticky model"""
    lims = kwds.pop('lims', limits(num))
    fudge = kwds.pop('fudge', 0.0)
    lims = ((1-fudge) * lims[0], (1+fudge) * lims[1])
    betas = np.linspace(*lims, count)
    svals = beta_to_s(betas)
    env = (num - 1) / (1 + (num - 1) * svals)
    if num == 2:
        return svals, 1 / (1 + svals), env
    epss = kwds.pop('eps', eps_star(betas, num))
    avals = sticky(epss, betas, num)
    return svals, avals, env


# -------------------------------------
def sticky_star(svals: np.ndarray, num: int) -> np.ndarray:
    """actual envelope"""
    betas = s_to_beta(svals)
    epss = eps_star(betas, num)
    avals = sticky(epss, betas, num)
    return avals


# -------------------------------------
def env_approx(svals: np.ndarray, num: int) -> np.ndarray:
    """approximate envelope"""
    return (num - 1) * (1 - np.sqrt(num * (num - 2) * svals))
