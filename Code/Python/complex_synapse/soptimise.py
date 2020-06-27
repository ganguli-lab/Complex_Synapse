# -*- coding: utf-8 -*-
"""Optimising synapse models
"""
from numbers import Number
import numpy as np
import scipy.optimize as sco
from sl_py_tools.iter_tricks import dcount, denumerate, delay_warnings, dzip
from sl_py_tools.numpy_tricks import markov_param as mp
import numpy_linalg as la
from .synapse_opt import SynapseOptModel
from . import optimise as opt

# =============================================================================
# Optimisation helper functions
# =============================================================================


def make_laplace_problem(rate: Number, nst: int, **kwds) -> (dict, dict):
    """Make an optimize problem.
    """
    modelopts = opt.get_model_opts(kwds)
    opt.get_param_opts(kwds)
    keep_feasible = kwds.pop('keep_feasible', False)
    method = kwds.get('method', 'SLSQP')
    if SynapseOptModel.type['serial'] or SynapseOptModel.type['ring']:
        raise ValueError('Shifted problem cannot have special topology')

    model = SynapseOptModel.rand(nst, **modelopts)
    fun = opt.make_loss_function(model, model.laplace_grad, None, True)
    x_init = model.get_params()
    hess = opt.make_loss_function(model, model.laplace_hess, None, True)

    bounds = sco.Bounds(0, 1 + rate, keep_feasible)
    ofun = opt.make_loss_function(model, model.peq_min_fun, rate)
    ojac = opt.make_loss_function(model, model.peq_min_grad, rate)
    if method not in {'SLSQP', 'trust-constr'}:
        raise ValueError('method must be one of SLSQP, trust-constr')
    if method == 'trust-constr':
        ohess = opt.make_loss_function(model, model.peq_min_hess, rate)
        constraint = sco.NonlinearConstraint(ofun, 0, 1, ojac, ohess,
                                             keep_feasible)
    else:
        hess = None
        constraint = {'type': 'ineq', 'args': (), 'fun': ofun, 'jac': ojac}

    problem = {'fun': fun, 'x0': x_init, 'jac': True, 'hess': hess,
               'bounds': bounds, 'constraints': constraint}
    problem.update(kwds)
    return problem


def optim_laplace(rate: Number, nst: int, **kwds) -> sco.OptimizeResult:
    """Optimised model at one value of rate
    """
    repeats = kwds.pop('repeats', 0)
    prob = make_laplace_problem(rate, nst, **kwds)
    res = sco.minimize(**prob)
    for _ in dcount('repeats', repeats, disp_step=1):
        opt.update_laplace_problem(prob)
        new_res = sco.minimize(**prob)
        if new_res.fun < res.fun:
            res = new_res
    return res


def optim_laplace_range(rates: np.ndarray, nst: int, **kwds) -> (la.lnarray,
                                                                 la.lnarray):
    """Optimised model at many values of rate
    """
    popts = opt.get_param_opts(kwds)
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
    popts = opt.get_param_opts(kwds)
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
