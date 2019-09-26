# -*- coding: utf-8 -*-
"""Optimising synapse modelopts
"""
from numbers import Number
import numpy as np
import scipy.optimize as sco
from sl_py_tools.iter_tricks import dcount, denumerate, delay_warnings
from .builders import la, mp
from .synapse_opt import SynapseOptModel

# =============================================================================
# %%* Optimisation helper functions
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

    def loss_function(params):
        """Loss function
        """
        model.set_params(params)
        return method(*args, **kwds)
    return loss_function


def make_laplace_problem(rate: Number, nst: int, **kwds) -> (dict, dict):
    """Make an optimize problem.
    """
    modelopts = kwds.pop('modelopts', {})
    modelopts.setdefault('frac', kwds.pop('frac', 0.5))
    modelopts.setdefault('binary', kwds.pop('binary', True))
    modelopts.setdefault('npl', 2)
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
    return problem, modelopts


def update_laplace_problem(problem: dict):
    """Update an optimize problem with new x_init.
    """
    problem['x0'] = la.random.rand(problem['x0'].size)


def optim_laplace(rate: Number, nst: int, **kwds) -> (SynapseOptModel,
                                                      sco.OptimizeResult):
    """Optimised model at one value of rate
    """
    res: sco.OptimizeResult
    SynapseOptModel.type['serial'] = kwds.pop('serial', False)
    if SynapseOptModel.type['serial']:
        SynapseOptModel.directions = (1, -1)
    else:
        SynapseOptModel.directions = (0, 0)
    repeats = kwds.pop('repeats', 0)
    prob, mopts = make_laplace_problem(rate, nst, **kwds)
    res = sco.minimize(**prob)
    for i in dcount('repeats', repeats, disp_step=1):
        update_laplace_problem(prob)
        new_res = sco.minimize(**prob)
        if new_res.fun < res.fun:
            res = new_res
    return SynapseOptModel.from_params(res.x, **mopts), res


def optim_laplace_range(rates: np.ndarray, nst: int, **kwds) -> (la.lnarray,
                                                                 la.lnarray):
    """Optimised model at many values of rate
    """
    snr = la.empty_like(rates)
    serl = kwds.get('serial', False)
    models = la.empty((len(rates), 2*mp.num_param(nst, serial=serl, drn=serl)))
    with delay_warnings():
        for i, rate in denumerate('rate', rates):
            res = optim_laplace(rate, nst, **kwds)[1]
            snr[i] = - res.fun
            models[i] - res.x
    return snr, models
