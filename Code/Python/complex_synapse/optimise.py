# -*- coding: utf-8 -*-
"""Optimising synapse modelopts
"""
from functools import wraps
from numbers import Number
from typing import Callable, Union, TypeVar, Tuple, Optional

import numpy as np
import scipy.optimize as sco

import numpy_linalg as la
from sl_py_tools.containers import listify
from sl_py_tools.arg_tricks import default
from sl_py_tools.iter_tricks import dcount, delay_warnings, denumerate, dzip
from sl_py_tools.numpy_tricks import markov_param as mp

from . import builders as bld
from .synapse_opt import SynapseOptModel

Data = TypeVar('Data', Number, np.ndarray)
Func = Callable[[np.ndarray], Data]
Constraint = Union[sco.LinearConstraint, sco.NonlinearConstraint]
Problem = Tuple[Func[Number], Func[np.ndarray], Func[np.ndarray],
                sco.Bounds, Constraint]
Maker = Callable[[Number, int], Problem]
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


def make_loss_function(model: SynapseOptModel, method: str,
                       *args, **kwds) -> Func:
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


def get_param_opts(opts: dict = None, **kwds) -> (int, int):
    """Make options dict for Markov processes.

    Returns
    -------
    nst : int or None
        total number of states. Use `npar` if `None`.
    npar : int or None
        total number of parameters. Use `nst` if `None` (default).
    """
    opts = default(opts, {})
    opts.update(kwds)
    # set opts
    SynapseOptModel.RCondThresh = opts.pop('RCondThresh', 1e-4)
    SynapseOptModel.type['serial'] = opts.pop('serial', False)
    SynapseOptModel.type['ring'] = opts.pop('ring', False)
    SynapseOptModel.type['uniform'] = opts.pop('uniform', False)
    if any(SynapseOptModel.type.values()):
        SynapseOptModel.directions = (1, -1)
    else:
        SynapseOptModel.directions = (0, 0)
    SynapseOptModel.directions = opts.pop('drn', SynapseOptModel.directions)
    # get opts
    nst, npar = None, None
    paramopts = SynapseOptModel.type.copy()
    paramopts['drn'] = SynapseOptModel.directions[0]
    if 'nst' in opts:
        nst = opts.pop('nst')
        npar = opts.pop('npar', 2 * mp.num_param(nst, **paramopts))
    if 'npar' in opts and not paramopts['uniform']:
        npar = opts.pop('npar')
        nst = opts.pop('nst', mp.num_state(npar // 2, **paramopts))
    return nst, npar


def get_model_opts(opts: dict = None, **kwds) -> dict:
    """Make options dict for SynapseOptModel.
    """
    opts = default(opts, {})
    opts.update(kwds)
    modelopts = opts.pop('modelopts', {})
    modelopts.setdefault('frac', opts.pop('frac', 0.5))
    modelopts.setdefault('binary', opts.pop('binary', True))
    modelopts.setdefault('npl', 2)
    return modelopts


def make_model(opts: dict = None, **kwds) -> SynapseOptModel:
    """Make options dict for Markov processes.

    Returns
    -------
    model : SynapseOptModel
        An instance to use in loss functions etc
    """
    opts.update(kwds)
    model = opts.pop('model', None)
    if model is not None:
        return SynapseOptModel
    modelopts = get_model_opts(opts)
    nst, npar = get_param_opts(opts)
    params = opts.pop('params', None)
    if params is not None:
        return SynapseOptModel.from_params(params, **modelopts)
    return SynapseOptModel.rand(nst=nst, npar=npar, **modelopts)


# =============================================================================
# Problem creators
# =============================================================================


def make_problem(maker: Maker, rate: Number, **kwds) -> dict:
    """Make an optimize problem.
    """
    model = make_model(kwds)

    feasible = kwds.pop('keep_feasible', False)
    method = kwds.get('method', 'SLSQP')
    if method not in {'SLSQP', 'trust-constr'}:
        raise ValueError('method must be one of SLSQP, trust-constr')

    x_init = model.get_params()
    fun, jac, hess, bounds, constraint = maker(model, rate, method, feasible)

    problem = {'fun': fun, 'x0': x_init, 'jac': jac, 'hess': hess,
               'bounds': bounds, 'constraints': [constraint]}
    if model.type['serial'] or model.type['ring']:
        del problem['constraints']
    problem.update(kwds)
    return problem


def make_normal_problem(model: SynapseOptModel, rate: Number, method: str,
                        keep_feasible: bool) -> Problem:
    """Make an optimize problem.
    """
    fun = make_loss_function(model, model.laplace_fun, rate)
    jac = make_loss_function(model, model.laplace_grad, rate)
    hess = None

    con_coeff = constraint_coeff(model)
    bounds = sco.Bounds(0, 1, keep_feasible)

    if method == 'trust-constr':
        hess = make_loss_function(model, model.laplace_hess, rate)
        constraint = sco.LinearConstraint(con_coeff, 0, 1, keep_feasible)
    else:
        constraint = {'type': 'ineq', 'args': (),
                      'fun': lambda x: 1 - con_coeff @ x,
                      'jac': lambda x: -con_coeff}

    return fun, jac, hess, bounds, constraint


# -----------------------------------------------------------------------------
# Shifting rate from function to constraints
# -----------------------------------------------------------------------------


def make_shifted_problem(model: SynapseOptModel, rate: Number, method: str,
                         keep_feasible: bool) -> Problem:
    """Make an optimize problem with rate shifted to constraints.
    """
    if SynapseOptModel.type['serial'] or SynapseOptModel.type['ring']:
        raise ValueError('Shifted problem cannot have special topology')

    fun = make_loss_function(model, model.laplace_fun, None, True)
    jac = make_loss_function(model, model.laplace_grad, None, True)
    hess = None

    bounds = sco.Bounds(0, 1 + rate, keep_feasible)
    ofun = make_loss_function(model, model.peq_min_fun, rate)
    ojac = make_loss_function(model, model.peq_min_grad, rate)

    if method == 'trust-constr':
        hess = make_loss_function(model, model.laplace_hess, None, True)
        ohess = make_loss_function(model, model.peq_min_hess, rate)
        constraint = sco.NonlinearConstraint(ofun, 0, 1, ojac, ohess,
                                             keep_feasible)
    else:
        constraint = {'type': 'ineq', 'args': (), 'fun': ofun, 'jac': ojac}

    return fun, jac, hess, bounds, constraint


# =============================================================================
# Optimisation helper functions
# =============================================================================


def update_laplace_problem(problem: dict):
    """Update an optimize problem with new x_init.
    """
    problem['x0'] = bld.RNG.random(problem['x0'].size)


def check_trust_constr(sol: np.ndarray, con: Constraint) -> (str, np.ndarray):
    """Verify that solution satisfies a constraint for method trust-constr"""
    if isinstance(con, sco.LinearConstraint):
        vals = con.A @ sol
    elif isinstance(con, sco.NonlinearConstraint):
        vals = con.fun(sol)
    else:
        raise TypeError(f"Unknown constraint type: {type(con)}")
    if con.lb == con.ub:
        return 'eq', vals - con.lb
    return 'ineq', np.r_[vals - con.lb, con.ub - vals]


def verify_solution(prob: dict, result: sco.OptimizeResult) -> bool:
    """Verify that solution satisfies constraints"""
    solution = result.x
    bounds = prob['bounds']
    if (solution < bounds.lb).any() or (solution > bounds.ub).any():
        return False
    for cons in listify(prob.get('constraints', [])):
        if isinstance(cons, dict):
            kind, vals = cons['type'], cons['fun'](solution)
        else:
            kind, vals = check_trust_constr(solution, cons)
        fail = (not np.allclose(vals, 0)) if kind == 'eq' else (vals < 0).any()
        if fail:
            return False
    return True


# =============================================================================
# Optimisation
# =============================================================================


def optim_laplace(rate: Number, nst: Optional[int] = None, *,
                  model: Optional[SynapseOptModel] = None,
                  maker: Maker = make_normal_problem,
                  **kwds) -> sco.OptimizeResult:
    """Optimised model at one value of rate
    """
    repeats = kwds.pop('repeats', 0)
    prob = make_problem(maker, rate, nst=nst, model=model, **kwds)
    res = sco.minimize(**prob)
    for _ in dcount('repeats', repeats, disp_step=1):
        update_laplace_problem(prob)
        new_res = sco.minimize(**prob)
        if verify_solution(prob, new_res) and new_res.fun < res.fun:
            res = new_res
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


def check_rcond_range(rates: np.ndarray, models: np.ndarray,
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
    rcnd = la.empty_like(rates)
    synmodel = make_model(kwds, params=models[0])
    with delay_warnings():
        for i, rate, params in denumerate('rate', rates, models):
            synmodel.set_params(params)
            rcnd[i] = synmodel.rcond(rate, rate_max=True)
    return rcnd


# =============================================================================
# Theory
# =============================================================================


def proven_envelope_laplace(rate: Data, nst: int) -> Data:
    """Theoretical envelope for Laplace transform"""
    return (nst - 1) / (rate * (nst - 1) + 1)
