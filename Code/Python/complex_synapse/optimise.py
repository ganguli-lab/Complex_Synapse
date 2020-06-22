# -*- coding: utf-8 -*-
"""Optimising synapse modelopts
"""
from numbers import Number
import numpy as np
import scipy.optimize as sco
from sl_py_tools.iter_tricks import dcount, denumerate, delay_warnings, dzip
from sl_py_tools.numpy_tricks import markov_param as mp
from sl_py_tools.arg_tricks import default
import numpy_linalg as la
from .synapse_opt import SynapseOptModel
from . import builders as bld

# =============================================================================
# Optimisation helper functions
# =============================================================================


def constraint_coeff(nst: int, npl: int = 2) -> la.lnarray:
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
    rows, cols = la.arange(npl*nst), la.arange(nst-1)
    cols = cols + rows.c * (nst-1)
    coeffs = la.zeros((npl*nst, npl*nst*(nst-1)))
    coeffs[rows.c, cols] = 1
    return coeffs


def make_loss_function(model: SynapseOptModel, method: str, *args, **kwds):
    """Make a loss function from model, method and value
    """
    if isinstance(method, str):
        method = getattr(model, method)

    def loss_function(*params):
        """Loss function
        """
        model.set_params(params[0])
        return method(*args, *params[1:], **kwds)
    return loss_function


def get_param_opts(opts: dict = None, **kwds) -> dict:
    """Make options dict for Markov processes.
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
    paramopts = SynapseOptModel.type.copy()
    paramopts['drn'] = SynapseOptModel.directions
    return paramopts


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


def make_laplace_problem(rate: Number, nst: int, **kwds) -> (dict, dict):
    """Make an optimize problem.
    """
    modelopts = get_model_opts(kwds)
    get_param_opts(kwds)
    keep_feasible = kwds.pop('keep_feasible', False)
    method = kwds.get('method', 'SLSQP')

    model = SynapseOptModel.rand(nst, **modelopts)
    fun = make_loss_function(model, model.laplace_grad, rate)
    x_init = model.get_params()
    hess = make_loss_function(model, model.laplace_hess, rate)

    con_coeff = constraint_coeff(nst, modelopts['npl'])
    bounds = sco.Bounds(0, 1, keep_feasible)
    if method not in {'SLSQP', 'trust-constr'}:
        raise ValueError('method must be one of SLSQP, trust-constr')
    if method == 'trust-constr':
        constraint = sco.LinearConstraint(con_coeff, 0, 1, keep_feasible)
    else:
        hess = None
        constraint = {'type': 'ineq', 'args': (),
                      'fun': lambda x: 1 - con_coeff @ x,
                      'jac': lambda x: -con_coeff}

    problem = {'fun': fun, 'x0': x_init, 'jac': True, 'hess': hess,
               'bounds': bounds, 'constraints': constraint}
    if model.type['serial'] or model.type['ring']:
        del problem['constraints']
    problem.update(kwds)
    return problem


def update_laplace_problem(problem: dict):
    """Update an optimize problem with new x_init.
    """
    problem['x0'] = bld.RNG.random(problem['x0'].size)


def optim_laplace(rate: Number, nst: int, **kwds) -> (SynapseOptModel,
                                                      sco.OptimizeResult):
    """Optimised model at one value of rate
    """
    repeats = kwds.pop('repeats', 0)
    prob = make_laplace_problem(rate, nst, **kwds)
    res = sco.minimize(**prob)
    for _ in dcount('repeats', repeats, disp_step=1):
        update_laplace_problem(prob)
        new_res = sco.minimize(**prob)
        if new_res.fun < res.fun:
            res = new_res
    return res


def optim_laplace_range(rates: np.ndarray, nst: int, **kwds) -> (la.lnarray,
                                                                 la.lnarray):
    """Optimised model at many values of rate
    """
    popts = get_param_opts(kwds)
    drn = popts['drn']
    popts['drn'] = drn[0]
    snr = la.empty_like(rates)
    models = la.empty((len(rates), 2 * mp.num_param(nst, **popts)))
    popts['drn'] = drn
    with delay_warnings():
        for i, rate in denumerate('rate', rates):
            res = optim_laplace(rate, nst, **kwds, **popts)
            snr[i] = - res.fun
            models[i] = res.x
    return snr, models


def reoptim_laplace_range(inds: np.ndarray, rates: np.ndarray,
                          models: np.ndarray, snr: np.ndarray, **kwds):
    """Optimised model at many values of rate
    """
    popts = get_param_opts(kwds)
    drn = popts['drn']
    popts['drn'] = drn[0]
    uniform = popts.pop('uniform', False)
    nst = mp.num_state(models.shape[1] // 2, **popts)
    popts['drn'] = drn
    popts['uniform'] = uniform
    with delay_warnings():
        for ind, rate in dzip('rate', inds, rates[inds]):
            res = optim_laplace(rate, nst, **kwds, **popts)
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
    popts = get_param_opts(kwds)
    mopts = get_model_opts(kwds)
    drn = popts['drn']
    popts['drn'] = drn[0]
    popts.pop('uniform', False)
    nst = mp.num_state(models.shape[1] // 2, **popts)
    # popts['drn'] = drn
    synmodel = SynapseOptModel.rand(nst, **mopts)
    with delay_warnings():
        for i, rate, params in denumerate('rate', rates, models):
            synmodel.set_params(params)
            rcnd[i] = min(synmodel.rcond(), synmodel.rcond(rate))
    return rcnd
