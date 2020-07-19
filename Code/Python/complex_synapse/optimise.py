# -*- coding: utf-8 -*-
"""Optimising synapse modelopts
"""
from functools import wraps
from numbers import Number
from typing import Callable, List, Optional, Tuple, TypeVar, Union

import numpy as np
import scipy.optimize as sco

import numpy_linalg as la
from sl_py_tools.arg_tricks import default
from sl_py_tools.containers import listify, map_join
from sl_py_tools.iter_tricks import (dcount, delay_warnings, denumerate, dzip,
                                     zenumerate)
from sl_py_tools.numpy_tricks.markov import params as mp

from . import builders as bld
from .synapse_opt import SynapseOptModel

Data = TypeVar('Data', Number, np.ndarray)
Func = Callable[[np.ndarray], Data]
Constraint = Union[sco.LinearConstraint, sco.NonlinearConstraint]
Problem = Tuple[Func[Number], Func[np.ndarray], Func[np.ndarray],
                sco.Bounds, List[Constraint]]
Maker = Callable[[SynapseOptModel, Number, bool, bool], Problem]
# =============================================================================
# Problem creation helper functions
# =============================================================================


def constraint_coeff(model: SynapseOptModel) -> la.lnarray:
    """Coefficient matrix for upper bound on off-diagonal row sums.

    Parameters
    ----------
    nst : int
        Number of states.
    npl : int
        Number of types of plasticity

    Returns
    -------
    coeffs : la.lnarray (2*nst, 2*nst(nst-1))
        matrix of coefficients s.t ``coeffs @ params <= 1``.
    """
    npl, nst = model.nplast, model.nstates
    rows, cols = la.arange(npl*nst), la.arange(nst-1)
    cols = cols + rows.c * (nst-1)
    coeffs = la.zeros((npl*nst, npl*nst*(nst-1)))
    coeffs[rows.c, cols] = 1
    return coeffs


def make_fn(model: SynapseOptModel, method: str, *args, **kwds) -> Func:
    """Make a loss function from model, method and value
    """
    if isinstance(method, str):
        method = getattr(model, method)

    @wraps(method)
    def loss_function(*params):
        """Loss function
        """
        model.set_params(params[0])
        return method(*args, *params[1:], **kwds)
    loss_function.model = model
    return loss_function


def cond_limit(model: SynapseOptModel, rate: Number,
               keep_feasible: bool = False, **kwds) -> sco.NonlinearConstraint:
    """Create a constraint on the condition number
    """
    kwds['svd'] = True
    cond_fn = make_fn(model, model.cond_fun, rate, **kwds)
    cond_jac = make_fn(model, model.cond_grad, rate, **kwds)
    return sco.NonlinearConstraint(cond_fn, 0, np.inf, cond_jac,
                                   keep_feasible=keep_feasible)


def conv_constraint(constraint: Constraint) -> List[dict]:
    """Convert constraint from trust-constr to SLSQP format
    """
    if isinstance(constraint, dict):
        return [constraint]

    slsqp_lb = {'type': 'ineq', 'args': ()}
    slsqp_ub = slsqp_lb.copy()
    if isinstance(constraint, sco.LinearConstraint):
        slsqp_lb['fun'] = lambda x: constraint.A @ x - constraint.lb
        slsqp_ub['fun'] = lambda x: constraint.ub - constraint.A @ x
        slsqp_lb['jac'] = lambda x: constraint.A
        slsqp_ub['jac'] = lambda x: - constraint.A
    elif isinstance(constraint, sco.NonlinearConstraint):
        slsqp_lb['fun'] = lambda x: constraint.fun(x) - constraint.lb
        slsqp_ub['fun'] = lambda x: constraint.ub - constraint.fun(x)
        if callable(constraint.jac):
            slsqp_lb['jac'] = constraint.jac
            slsqp_ub['jac'] = lambda x: - constraint.jac(x)
    else:
        raise TypeError(f"Unknown constraint type: {type(constraint)}")

    if not np.isfinite(constraint.lb):
        return [slsqp_ub]
    if not np.isfinite(constraint.ub):
        return [slsqp_lb]
    if constraint.lb == constraint.ub:
        slsqp_lb['type'] = 'eq'
        return [slsqp_lb]
    return [slsqp_lb, slsqp_ub]


def set_param_opts(opts: Optional[dict] = None):
    """Set options dict for SynapseOptModel Markov processes.
    """
    opts = default(opts, {})
    # set opts
    SynapseOptModel.CondThresh = opts.pop('CondThresh', 1e4)
    SynapseOptModel.type['serial'] = opts.pop('serial', False)
    SynapseOptModel.type['ring'] = opts.pop('ring', False)
    SynapseOptModel.type['uniform'] = opts.pop('uniform', False)
    if any(SynapseOptModel.type.values()):
        SynapseOptModel.directions = (1, -1)
    else:
        SynapseOptModel.directions = (0, 0)
    SynapseOptModel.directions = opts.pop('drn', SynapseOptModel.directions)


def get_model_opts(opts: Optional[dict] = None) -> dict:
    """Make options dict for SynapseOptModel instances.

    Returns
    -------
    modelopts : dict
        Options for `Synapse*Model` instances.
    """
    opts = default(opts, {})
    modelopts = opts.pop('modelopts', {})
    modelopts.setdefault('frac', opts.pop('frac', 0.5))
    modelopts.setdefault('binary', opts.pop('binary', True))
    modelopts.setdefault('npl', len(SynapseOptModel.directions))
    return modelopts


def make_model(opts: Optional[dict] = None, **kwds) -> SynapseOptModel:
    """Make options dict for Markov processes.

    `kwds` are added to `opts`, then all model related options are popped.

    Returns
    -------
    model : SynapseOptModel
        An instance to use in loss functions etc
    """
    opts = default(opts, {})
    opts.update(kwds)
    model = opts.pop('model', None)
    params = opts.pop('params', None)
    nst, npar = opts.pop('nst', None), opts.pop('npar', None)
    if model is not None:
        return model
    set_param_opts(opts)
    modelopts = get_model_opts(opts)
    if params is not None:
        return SynapseOptModel.from_params(params, **modelopts)
    paramopts = SynapseOptModel.type.copy()
    paramopts['drn'] = SynapseOptModel.directions[0]
    npl = modelopts.get('npl', 2)
    if nst is not None:
        npar = default(npar, npl * mp.num_param(nst, **paramopts))
    if npar is not None and not paramopts['uniform']:
        nst = default(nst, mp.num_state(npar // npl, **paramopts))
    if None in {nst, npar}:
        raise TypeError("Must specify one of [model, params, nst, npar]")
    return SynapseOptModel.rand(nst=nst, npar=npar, **modelopts)


# =============================================================================
# Problem creators
# =============================================================================


def make_problem(maker: Maker, rate: Number, **kwds) -> dict:
    """Make an optimize problem.
    """
    model = make_model(kwds)

    method = kwds.get('method', 'SLSQP')
    if method not in {'SLSQP', 'trust-constr'}:
        raise ValueError('method must be one of SLSQP, trust-constr')
    opts = {'keep_feasible': kwds.pop('keep_feasible', False),
            'inv': method == 'trust-constr',
            'svd': kwds.pop('cond', False)}

    x_init = model.get_params()
    fun, jac, hess, bounds, constraints = maker(model, rate, **opts)
    if not opts['inv']:
        constraints = map_join(conv_constraint, constraints)

    problem = {'fun': fun, 'x0': x_init, 'jac': jac, 'hess': hess,
               'bounds': bounds, 'constraints': constraints}
    if any(model.type.values()):
        del problem['constraints']
    problem.update(kwds)
    return problem


def normal_problem(model: SynapseOptModel, rate: Number, inv: bool = False,
                   keep_feasible: bool = False, **kwds) -> Problem:
    """Make an optimize problem.
    """
    fun = make_fn(model, model.laplace_fun, rate, inv=inv, **kwds)
    jac = make_fn(model, model.laplace_grad, rate, inv=inv, **kwds)
    hess = make_fn(model, model.laplace_hess, rate, **kwds) if inv else None

    con_coeff = constraint_coeff(model)
    bounds = sco.Bounds(0, 1, keep_feasible)
    diag = [sco.LinearConstraint(con_coeff, -np.inf, 1, keep_feasible)]
    if kwds.get('svd', False):
        diag.append(cond_limit(model, rate, keep_feasible, inv=inv))

    return fun, jac, hess, bounds, diag


# -----------------------------------------------------------------------------
# Shifting rate from function to constraints
# -----------------------------------------------------------------------------


def shifted_problem(model: SynapseOptModel, rate: Number, inv: bool = False,
                    keep_feasible: bool = False, **kwds) -> Problem:
    """Make an optimize problem with rate shifted to constraints.
    """
    if any(SynapseOptModel.type.values()):
        raise ValueError('Shifted problem cannot have special topology')

    fun = make_fn(model, model.area_fun, inv=True, **kwds)
    jac = make_fn(model, model.area_grad, inv=True, **kwds)
    hess = make_fn(model, model.area_hess, **kwds) if inv else None

    bounds = None
    cfun = make_fn(model, model.peq_min_fun, rate, **kwds)
    cjac = make_fn(model, model.peq_min_grad, rate, **kwds)
    chess = make_fn(model, model.peq_min_hess, rate, **kwds) if inv else None
    lims = sco.NonlinearConstraint(cfun, 0, np.inf, cjac, chess, keep_feasible)
    if kwds.get('svd', False):
        lims = [lims, cond_limit(model, None, keep_feasible, inv=True)]

    return fun, jac, hess, bounds, listify(lims)


# =============================================================================
# Optimisation helper functions
# =============================================================================


def update_laplace_problem(problem: dict):
    """Update an optimize problem with new x_init.
    """
    problem['x0'] = bld.RNG.random(problem['x0'].size)


def check_trust_constr(sol: np.ndarray, con: Constraint) -> (bool, np.ndarray):
    """Verify that solution satisfies a constraint for method trust-constr"""
    if isinstance(con, sco.LinearConstraint):
        vals = con.A @ sol
    elif isinstance(con, sco.NonlinearConstraint):
        vals = con.fun(sol)
    else:
        raise TypeError(f"Unknown constraint type: {type(con)}")
    if not np.isfinite(con.ub):
        return vals - con.lb, False
    if not np.isfinite(con.lb):
        return con.ub - vals, False
    if con.lb == con.ub:
        return vals - con.lb, True
    return np.r_[vals - con.lb, con.ub - vals], False


def constr_violation(prob: dict, result: sco.OptimizeResult) -> bool:
    """Verify that solution satisfies constraints"""
    maxcv = getattr(result, 'constr_violation', None)
    if maxcv is not None:
        return maxcv
    # must be SLSQP
    maxcv = 0
    solution = result.x
    bounds = prob['bounds']
    if bounds is not None:
        maxcv = max(maxcv, (solution - bounds.ub).max())
        maxcv = max(maxcv, (bounds.lb - solution).max())
    for cons in listify(prob.get('constraints', [])):
        kind, vals = cons['type'] == 'eq', cons['fun'](solution)
        maxcv = max(maxcv, np.fabs(vals).max() if kind else - vals.min())
    return maxcv


def verify_solution(prob: dict, result: sco.OptimizeResult) -> bool:
    """Verify that solution satisfies constraints"""
    # maxcv = constr_violation(prob, result)
    # return maxcv < prob.get('tol', 1e-3)
    itol = prob.get('tol', 1e-3)
    maxcv = getattr(result, 'constr_violation', None)
    if maxcv is not None:
        return maxcv < itol
    # must be SLSQP
    maxcv = 0
    solution = result.x
    bounds = prob['bounds']
    if bounds is not None:
        if (solution < bounds.lb).any() or (solution > bounds.ub).any():
            return False
    for cons in listify(prob.get('constraints', [])):
        vals, kind = cons['fun'](solution), cons['type'] == 'eq'
        fail = (not np.allclose(vals, 0)) if kind else (vals < -itol).any()
        if fail:
            return False
    return True


# =============================================================================
# Optimisation
# =============================================================================


def first_good(prob: dict) -> sco.OptimizeResult:
    """First solution that satisfies constraints"""
    max_tries = prob.pop('max_tries', 100)
    for _ in dcount('tries', max_tries):
        res = sco.minimize(**prob)
        if verify_solution(prob, res):
            break
        update_laplace_problem(prob)
    return res


def optim_laplace(rate: Number, nst: Optional[int] = None, *,
                  model: Optional[SynapseOptModel] = None,
                  maker: Maker = normal_problem, **kwds) -> sco.OptimizeResult:
    """Optimised model at one value of rate
    """
    repeats = kwds.pop('repeats', 0)
    prob = make_problem(maker, rate, nst=nst, model=model, **kwds)
    res = first_good(prob)
    for _ in dcount('repeats', repeats, disp_step=1):
        update_laplace_problem(prob)
        new_res = sco.minimize(**prob)
        if verify_solution(prob, new_res) and new_res.fun < res.fun:
            res = new_res
    if not verify_solution(prob, res):
        res.fun = np.nan
    return res


def optim_laplace_range(rates: np.ndarray, nst: int,
                        **kwds) -> (la.lnarray, la.lnarray):
    """Optimised model at many values of rate
    """
    model = make_model(kwds, nst=nst)
    snr = la.empty_like(rates)
    models = la.empty((len(rates), model.nparam))
    with delay_warnings():
        for i, rate in denumerate('rate', rates):
            res = optim_laplace(rate, model=model, **kwds)
            snr[i] = - res.fun
            models[i] = res.x
    return snr, models


def reoptim_laplace_range(inds: np.ndarray, rates: np.ndarray,
                          snr: np.ndarray, models: np.ndarray,
                          **kwds) -> (la.lnarray, la.lnarray):
    """Optimised model at many values of rate
    """
    model = make_model(kwds, params=models[inds[0]])
    with delay_warnings():
        for ind, rate in dzip('rate', inds, rates[inds]):
            res = optim_laplace(rate, model=model, **kwds)
            snr[ind] = - res.fun
            models[ind] = res.x
    return snr, models


def check_cond_range(rates: np.ndarray, models: np.ndarray,
                     **kwds) -> la.lnarray:
    """Inverse condition numbers of optiomised models

    Parameters
    ----------
    rates : np.ndarray (S,)
        inverse time (Laplace parameter)
    models : np.ndarray (S, P)
        optimised models

    Returns
    -------
    rcnd : np.ndarray (S,)
        inverse condition number, worst of :math:`Z(0), Z(s)`
    """
    cnd = la.empty_like(rates)
    model = make_model(kwds, params=models[0])
    for i, rate, params in zenumerate(rates, models):
        model.set_params(params)
        cnd[i] = model.cond(rate, rate_max=True)
    return cnd


# =============================================================================
# Theory
# =============================================================================


def proven_envelope_laplace(rate: Data, nst: int) -> Data:
    """Theoretical envelope for Laplace transform"""
    return (nst - 1) / (rate * (nst - 1) + 1)
