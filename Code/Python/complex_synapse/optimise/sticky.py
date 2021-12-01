"""Properties of optimal shortened models
"""
from typing import Tuple

import numpy as np
import scipy.optimize as _sco

from .shorten import Y_STAR, Values
# =============================================================================


def beta_to_s(beta: Values) -> Values:
    """Convert alpha to s"""
    return np.cosh(beta) - 1


def s_to_beta(sval: Values) -> Values:
    """Convert s to beta. positive one of two possibilities"""
    return np.arccosh(sval + 1)


# -------------------------------------
def beta_star(num: int) -> float:
    """Value of beta at which num is optimal"""
    return 2 * Y_STAR / num


# =============================================================================


def _components(beta: Values, numh: int) -> Tuple[Values, ...]:
    """Components of numerator and denominator

    numerator = numer - eps * dnumer

    denominator = denom - eps * ddenom

    Parameters
    ----------
    numh
        num/2

    Returns
    -------
    numer
        cosh(num*beta/2) - 1
    denom
        cosh(num*beta/2)
    dnumer
        cosh((num-2)*beta/2) - 1
    ddenom
        cosh((num-2)*beta/2)
    """
    denom = np.cosh(numh * beta)
    numer = denom - 1
    ddenom = np.cosh((numh - 1) * beta)
    dnumer = ddenom - 1
    return numer, denom, dnumer, ddenom


# -------------------------------------
def sticky(eps: Values, beta: Values, num: int) -> Values:
    """Laplace-SNR for sticky serial model"""
    sval = beta_to_s(beta)
    numer, denom, dnumer, ddenom = _components(beta, num / 2)
    fnumer = numer - eps * dnumer
    fdenom = denom - eps * ddenom
    fmu = num - eps * (num - 2)
    return 2 * (1 - eps) * fnumer / (fmu * sval * fdenom)


# -------------------------------------
def eps_stars(beta: Values, num: int) -> Tuple[Values, Values]:
    """Optimal epsilon for sticky serial model, raw"""
    if num == 2:
        return 0, 0
    numer, denom, dnumer, ddenom = _components(beta, num / 2)
    fnumer = numer / dnumer
    fdenom = denom / ddenom
    fmu = num / (num - 2)
    if np.isscalar(fnumer):
        fnumer = fmu**2 if np.isclose(0, beta) else fnumer
    else:
        fnumer[np.isclose(0, np.array(beta))] = fmu**2
    disc = np.sqrt((fmu - 1)*(fdenom - 1)*(fmu - fnumer)*(fdenom - fnumer))
    quad = 1 + fnumer - fdenom - fmu
    linr = (fnumer - fdenom * fmu)
    return (linr - disc) / quad, (linr + disc) / quad


def eps_star(beta: Values, num: int) -> Values:
    """Optimal epsilon for sticky serial model, clipped"""
    epss = eps_stars(beta, num)[0]
    return np.clip(epss, 0, 1)


def eps_star_star(beta: Values, num: int) -> Values:
    """Optimal epsilon for sticky serial model, clipped, other solution"""
    epss = eps_stars(beta, num)[1]
    return np.clip(epss, 0, 1)


# -------------------------------------
def _env_actual(betas: Values, svals: Values, numh: int) -> Values:
    """actual heuristic envelope

    Parameters
    ----------
    numh
        num/2
    """
    numer, _, dnumer, _ = _components(betas, numh)
    fdm = np.sqrt((numh - 1)*numer - numh*dnumer) + 1
    core = (numer - dnumer) / fdm**2
    return core / svals


def sticky_star(betas: Values, num: int) -> Values:
    """actual heuristic envelope"""
    svals = beta_to_s(betas)
    return _env_actual(betas, svals, num/2)


def sticky_star_s(svals: Values, num: int) -> Values:
    """actual heuristic envelope"""
    betas = s_to_beta(svals)
    return _env_actual(betas, svals, num/2)


# -------------------------------------
def env_approx(svals: Values, num: int) -> Values:
    """approximate heuristic envelope"""
    return (num - 1) * (1 - np.sqrt(num * (num - 2) * svals))


# =============================================================================


def lower(beta: float, num: int) -> float:
    """derivative of sticky wrt eps at eps==0
    """
    numer, denom, dnumer, ddenom = _components(beta, num / 2)
    # d/d(eps) (1-eps)(n - eps dn)/(m - (m-2) eps)(d - eps dd)
    return num * (denom*dnumer - numer*ddenom) + 2 * numer * denom


# -------------------------------------
def limits(num: int, debug: bool = False) -> Tuple[float, float]:
    """range of beta where sticky model can be optimised
    Note: higher t -> lower s -> higher M
    """
    if num == 2:
        return (0.01, np.exp(-Y_STAR))
    x_0, x_1 = beta_star(num), beta_star(num + 2)
    lo_sol = _sco.root_scalar(lower, args=(num,), x0=x_0, x1=x_1)
    if debug:
        print(lo_sol.flag, eps_star(lo_sol.root, num))
    return lo_sol.root, 0


# -------------------------------------
def envelope(num: int, count: int, **kwds) -> Tuple[np.ndarray, ...]:
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
