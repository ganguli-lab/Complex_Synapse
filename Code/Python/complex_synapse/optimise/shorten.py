"""Properties of optimal shortened models

alpha = exp(beta) from the notes
"""
from typing import Optional, Tuple, TypeVar

import numpy as np
import scipy.optimize as _sco

Values = TypeVar('Values', float, np.ndarray)
# =============================================================================
Y_STAR = _sco.root_scalar(lambda x: np.sinh(x) - np.tanh(x/2) - x,
                          bracket=(1, 2)).root
A_STAR = np.sinh(Y_STAR) / np.cosh(Y_STAR)**2


# -------------------------------------
def alpha_to_s(alpha: Values) -> Values:
    """Convert alpha to s"""
    return (alpha - 1)**2 / (2 * alpha)


def s_to_alpha(sval: Values) -> Tuple[Values, Values]:
    """Convert s to alpha, two possibilities, reciprocals"""
    sval = sval + 1
    disc = np.sqrt(sval**2 - 1)
    return sval - disc, sval + disc


def beta_to_s(beta: Values) -> Values:
    """Convert beta to s

    .. math:: s = 2 \\sinh^2(\\beta / 2)

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.

    Returns
    -------
    s : float|ndarray
        Parameter of Laplace transform.
    """
    return np.cosh(beta) - 1


def s_to_beta(sval: Values) -> Values:
    """Convert s to beta. positive one of two possibilities

    .. math:: s = 2 \\sinh^2(\\beta / 2)

    Parameters
    ----------
    s : float|ndarray
        Parameter of Laplace transform.

    Returns
    -------
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    """
    return np.arccosh(sval + 1)


# =============================================================================


def uniform(beta: Values, num: int) -> Values:
    """Laplace-SNR for uniform serial model

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    num : int
        Number of states

    Returns
    -------
    snr_laplace : float|ndarray
        Laplace transform of SNR curve.
    """
    sval = beta_to_s(beta)
    smp1 = np.cosh(beta * num / 2)
    return 2 * (smp1 - 1) / (num * sval * smp1)


# -------------------------------------
def num_star(beta: Values) -> Values:
    """Optimal num for given beta (continuous)

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.

    Returns
    -------
    num : int|float|ndarray
        Number of states
    """
    return 2 * Y_STAR / np.abs(beta)


def beta_star(num: int) -> float:
    """Beta where num is optimal

    Parameters
    ----------
    num : int|float|ndarray
        Number of states

    Returns
    -------
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    """
    return 2 * Y_STAR / num


def _uni_star(beta: Values, svals: Values) -> Values:
    """Heuristic envelope from varying M

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    s : float|ndarray
        Parameter of Laplace transform.

    Returns
    -------
    snr_laplace : float|ndarray
        Laplace transform of SNR curve.
    """
    return A_STAR * np.abs(beta) / svals


# -------------------------------------
def s_star(num: int) -> float:
    """s where num is optimal

    Parameters
    ----------
    num : int|float|ndarray
        Number of states

    Returns
    -------
    s : float|ndarray
        Parameter of Laplace transform.
    """
    beta = beta_star(num)
    return beta_to_s(beta)


def uni_star(beta: Values) -> Values:
    """Heuristic envelope from varying M for uniform model.

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.

    Returns
    -------
    snr_laplace : float|ndarray
        Laplace transform of SNR curve.
    """
    return _uni_star(beta, beta_to_s(beta))


def uni_star_s(svals: Values) -> Values:
    """Heuristic envelope from varying M for uniform model.

    Parameters
    ----------
    s : float|ndarray
        Parameter of Laplace transform.

    Returns
    -------
    snr_laplace : float|ndarray
        Laplace transform of SNR curve.
    """
    return _uni_star(s_to_beta(svals), svals)


# =============================================================================


# definition of alpha does matter for next 1
def _components(beta: Values, num: int, sval: Optional[Values] = None) -> Tuple[Values, ...]:
    """Components of numerator and denominator

    numerator = numer - eps * dnumer
    denominator = denom - eps * ddenom

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    num : int
        Number of states
    sval : Optional[float|ndarray], optional
        Parameter of Laplace transform, or calculate from beta by default None

    Returns
    -------
    components : Tuple[float|ndarray, ...]
        numer, denom, dnumer, ddenom
    """
    sval = beta_to_s(beta) if sval is None else sval
    denom = np.cosh(beta * num / 2)
    delta = np.cosh(beta * (num - 2) / 2)
    numer = denom - 1
    fac = 2 * sval + 1
    dnumer = numer - fac * (delta - 1)
    ddenom = denom - fac * delta
    return numer, denom, dnumer, ddenom


# -------------------------------------
def _short(eps: Values, beta: Values, sval: Values, num: int) -> Values:
    """Laplace-SNR for shortened serial model

    Parameters
    ----------
    eps : float|ndarray
        Transition rate out to end states.
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    sval : float|ndarray
        Parameter of Laplace transform.
    num : int
        Number of states

    Returns
    -------
    snr_laplace : float|ndarray
        Laplace transform of SNR curve.
    """
    if num == 2:
        return (1 - eps) / (1 - eps + sval)
    numer, denom, dnumer, ddenom = _components(beta, num, sval)
    fnumer = numer - eps * dnumer
    fdenom = denom - eps * ddenom
    return 2 * fnumer / ((num - 2*eps) * sval * fdenom)


def short(eps: Values, beta: Values, num: int) -> Values:
    """Laplace-SNR for shortened serial model

    Parameters
    ----------
    eps : float|ndarray
        Transition rate out to end states.
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    num : int
        Number of states

    Returns
    -------
    snr_laplace : float|ndarray
        Laplace transform of SNR curve.
    """
    return _short(eps, beta, beta_to_s(beta), num)


def short_s(eps: Values, sval: Values, num: int) -> Values:
    """Laplace-SNR for shortened serial model

    Parameters
    ----------
    eps : float|ndarray
        Transition rate out to end states.
    sval : float|ndarray
        Parameter of Laplace transform.
    num : int
        Number of states

    Returns
    -------
    snr_laplace : float|ndarray
        Laplace transform of SNR curve.
    """
    return _short(eps, s_to_beta(sval), sval, num)


# -------------------------------------
def eps_stars(beta: Values, num: int) -> Tuple[Values, Values]:
    """Optimal epsilon for shortened serial model, both solutions, unclipped.

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    num : int
        Number of states

    Returns
    -------
    eps : tuple[float|ndarray,float|ndarray]
        Transition rate out to end states.
    """
    if num == 2:
        return 0, 0
    numer, denom, dnumer, ddenom = _components(beta, num)
    fnumer = numer / dnumer
    fdenom = denom / ddenom
    disc2 = (num/2 - fnumer) * (fdenom - fnumer)
    # if np.any(disc2 < 0):
    #     print(f"{disc2=}")
    disc = np.sqrt((num/2 - fnumer) * (fdenom - fnumer))
    return fnumer - disc, fnumer + disc


def eps_star(beta: Values, num: int) -> Values:
    """Optimal epsilon for shortened serial model, 1st solution, clipped.

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    num : int
        Number of states

    Returns
    -------
    eps : float|ndarray
        Transition rate out to end states.
    """
    epss = eps_stars(beta, num)[0]
    return np.clip(epss, 0, 1)


def eps_star_star(beta: Values, num: int) -> Values:
    """Optimal epsilon for shortened serial model, 2nd solution, clipped.

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    num : int
        Number of states

    Returns
    -------
    eps : float|ndarray
        Transition rate out to end states.
    """
    epss = eps_stars(beta, num)[1]
    return np.clip(epss, 0, 1)


# -------------------------------------
def short_star(betas: Values, num: int) -> Values:
    """Actual heuristic envelope from optimal epsilon

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    num : int
        Number of states

    Returns
    -------
    snr_laplace : float|ndarray
        Laplace transform of SNR curve.
    """
    svals = beta_to_s(betas)
    if num == 2:
        return 1 / (1 + svals)
    epss = eps_star(betas, num)
    return _short(epss, betas, svals, num)


def short_star_s(svals: Values, num: int) -> Values:
    """Actual heuristic envelope from optimal epsilon

    Parameters
    ----------
    svals : float|ndarray
        Parameter of Laplace transform.
    num : int
        Number of states

    Returns
    -------
    snr_laplace : float|ndarray
        Laplace transform of SNR curve.
    """
    if num == 2:
        return 1 / (1 + svals)
    betas = s_to_beta(svals)
    epss = eps_star(betas, num)
    return _short(epss, betas, svals, num)


# -------------------------------------
def env_approx(svals: Values) -> Values:
    """approximate heuristic envelope

    Parameters
    ----------
    svals : float|ndarray
        Parameter of Laplace transform.

    Returns
    -------
    snr_laplace : float|ndarray
        Approximate Laplace transform of SNR curve.
    """
    return A_STAR * np.sqrt(2 / svals)


# =============================================================================


def lower(beta: float, num: int) -> float:
    """derivative of `short` wrt eps at eps==1
    Note: lower t -> higher s -> (higher, lower) beta[:],

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    num : int
        Number of states
    """
    # S(m), S(m)+1, S(m)-z(S(m-1)), S(m)+1-z(S(m-1)+1)
    numer, denom, dnumer, ddenom = _components(beta, num)
    # z(S(m-1))
    numer -= dnumer
    # z(S(m-1)+1)
    denom -= ddenom
    return (num - 2) * (denom * dnumer - numer * ddenom) - 2 * numer * denom


# -------------------------------------
def upper(beta: float, num: int) -> float:
    """derivative of `short` wrt eps at eps==0
    Note: higher t -> lower s -> (lower, higher) beta[:],

    Parameters
    ----------
    beta : float|ndarray
        Wavenumber across states for eta, etc.
    num : int
        Number of states
    """
    # S(m), S(m)+1, S(m)-z(S(m-1)), S(m)+1-z(S(m-1)+1)
    numer, denom, dnumer, ddenom = _components(beta, num)
    return num * (denom * dnumer - numer * ddenom) - 2 * numer * denom


# -------------------------------------
def limits(num: int, debug: bool = False) -> Tuple[float, float]:
    """range of beta where shortened model can be optimised
    Note: higher t -> lower s -> lower beta -> higher M

    Parameters
    ----------
    num : int
        Number of states
    debug : bool, optional
        Print diagnostic info, in case initial guesses fails to bracket roots,
        by default `False`.

    Returns
    -------
    beta1 : float|ndarray
        Wavenumber across states for eta, etc. Lower `t` side.
    beta2 : float|ndarray
        Wavenumber across states for eta, etc. higher `t` side.
    """
    if num == 2:
        return np.arccosh(100), beta_star(2)
    x_0, x_1 = 0.5 * beta_star(num + 2), 2 * beta_star(num - 2)
    if debug:
        print(f"{num=}-lo: f(x_0)={lower(x_0, num):.2f};"
              f" f(x_1)={lower(x_1, num):.2f};")
    lo_sol = _sco.root_scalar(lower, args=(num,), bracket=(x_0, x_1))
    if debug:
        print(f" {lo_sol.flag};"
              f" eps={eps_star(lo_sol.root, num):.2f};"
              f" f(x)={lower(lo_sol.root, num):.2f}")
        print(f"hi: f(x_0)={upper(x_0, num):.2f};"
              f" f(x_1)={upper(x_1, num):.2f};")
    hi_sol = _sco.root_scalar(upper, args=(num,), bracket=(x_0, x_1))
    if debug:
        print(f" {hi_sol.flag};"
              f" eps={eps_star(hi_sol.root, num):.2f};"
              f" f(x)={upper(hi_sol.root, num):.2f}")
    return lo_sol.root, hi_sol.root


# -------------------------------------
def envelope(num: int, count: int, **kwds) -> Tuple[np.ndarray, ...]:
    """Optimised shortened model, within limits

    Parameters
    ----------
    num : int
        Number of states
    count : int
        Number of points between limits where we optimise

    Returns
    -------
    svals : ndarray (count,)
        Parameter of Laplace transform.
    avals : ndarray (count,)
        Laplace transform of SNR curve. Actual maximum for given num.
    env : ndarray (count,)
        Laplace transform of SNR curve. Approximate maximum across nums.
    """
    lims = kwds.pop('lims', limits(num))
    fudge = kwds.pop('fudge', 0.0)
    lims = ((1-fudge) * lims[0], (1+fudge) * lims[1])
    betas = np.linspace(*lims, count)
    svals = beta_to_s(betas)
    env = uni_star(betas)
    if num == 2:
        return svals, 1 / (1 + svals), env
    epss = eps_star(betas, num)
    avals = _short(epss, betas, svals, num)
    return svals, avals, env


# -------------------------------------
def envelopes(num_max: int, count: int) -> Tuple[np.ndarray, ...]:
    """Set of optimised shortened models, within limits for all even nums

    Parameters
    ----------
    num_max : int
        Maximum number of states
    count : int
        Number of points between limits where we optimise

    Returns
    -------
    svals : ndarray (num_max/2,count)
        Parameter of Laplace transform.
    avals : ndarray (num_max/2,count)
        Laplace transform of SNR curve. Actual maximum for given num.
    env : ndarray (num_max/2,count)
        Laplace transform of SNR curve. Approximate maximum across nums.
    """
    siz = (num_max // 2, count)
    svals, avals, env = np.empty(siz), np.empty(siz), np.empty(siz)
    for i in range(num_max // 2):
        svals[i], avals[i], env[i] = envelope(2*i + 2, count)
    return svals, avals, env


# -------------------------------------
def gaps(num_max: int, count: int) -> Tuple[np.ndarray, ...]:
    """set of uniform models in the gaps between upper side for lower num and
    lower side of greater num.

    Parameters
    ----------
    num_max : int
        Maximum number of states
    count : int
        Number of points between limits where we compute

    Returns
    -------
    svals : ndarray (num_max/2 - 1,count)
        Parameter of Laplace transform.
    avals : ndarray (num_max/2 - 1,count)
        Laplace transform of SNR curve. Actual value for uniform model with
        a given num.
    env : ndarray (num_max/2 - 1,count)
        Laplace transform of SNR curve. Approximate maximum across nums.
    """
    gaps = np.full((num_max//2 - 1, 3, count), np.nan)
    lims = np.stack([limits(n) for n in range(2, num_max+2, 2)])
    for i, (high, low) in enumerate(zip(lims[:-1, 1], lims[1:, 0])):
        if high > low:
            gaps[i] = envelope(2*i + 2, count, lims=(low, high))
    return tuple(gaps.swapaxes(0, 1))


# -------------------------------------
def env_ravel(svals: np.ndarray, env: np.ndarray) -> Tuple[np.ndarray, ...]:
    """ravel and sort env
    """
    new_s, new_env = svals.ravel(), env.ravel()
    inds = np.argsort(new_s)
    return new_s[inds], new_env[inds]
