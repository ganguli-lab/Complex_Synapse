# -*- coding: utf-8 -*-
"""Optimising synapse modelopts
"""
from __future__ import annotations

from functools import wraps
from numbers import Number
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, TypeVar, Union

import numpy as np
import scipy.optimize as sco

import numpy_linalg as la
import sl_py_tools.containers as _cn
# import sl_py_tools.dict_tricks as _dt
# import sl_py_tools.arg_tricks as _ag
import sl_py_tools.iter_tricks as _it
import sl_py_tools.numpy_tricks.markov as _ma
import sl_py_tools.numpy_tricks.markov.params as _mp
import sl_py_tools.options_classes as _opt

import complex_synapse.builders as _bld
import complex_synapse.optimise.shorten as _sh
import complex_synapse.optimise.sticky as _st
# import complex_synapse.options as _opt
import complex_synapse.synapse_base as _sb
from complex_synapse.optimise.synapse_opt import SynapseOptModel
# =============================================================================
# Class for specifying types of synapses
# =============================================================================


# pylint: disable=too-many-ancestors
# =============================================================================
# Better LinearConstraint
# =============================================================================


class LinearConstraint(sco.LinearConstraint):
    """Extending Linear constraint to have callables like NonlinearConstraint

    Parameters
    ----------
    A : {array_like, sparse matrix}, shape (m, n)
        Matrix defining the constraint.
    lb, ub : array_like
        Lower and upper bounds on the constraint. Each array must have the
        shape (m,) or be a scalar, in the latter case a bound will be the same
        for all components of the constraint. Use ``np.inf`` with an
        appropriate sign to specify a one-sided constraint.
        Set components of `lb` and `ub` equal to represent an equality
        constraint. Note that you can mix constraints of different types:
        interval, one-sided or equality, by setting different components of
        `lb` and `ub` as  necessary.
    keep_feasible : array_like of bool, optional
        Whether to keep the constraint components feasible throughout
        iterations. A single value set this property for all components.
        Default is False. Has no effect for equality constraints.

    Methods
    -------
    fun
        The function defining the constraint.
        The signature is ``fun(x) -> array_like, shape (m,)``.
    jac
        Method of computing the Jacobian matrix (an m-by-n matrix,
        where element (i, j) is the partial derivative of f[i] with
        respect to x[j]). It must have the following signature:
        ``jac(x) -> ndarray, shape (m, n)``.
        Default is '2-point'.
    hess
        Method for computing the Hessian matrix.
        It must return the Hessian matrix of ``dot(fun, v)`` and
        must have the following signature:
        ``hess(x, v) -> ndarray, shape (n, n)``.
        Here ``v`` is ndarray with shape (m,) containing Lagrange multipliers.
    """
    def fun(self, x: np.ndarray) -> np.ndarray:
        """The function defining the constraint.
        The signature is ``fun(x) -> array_like, shape (m,)``.
        """
        return self.A @ x

    def jac(self, x: np.ndarray) -> np.ndarray:
        """Method of computing the Jacobian matrix (an m-by-n matrix,
        where element (i, j) is the partial derivative of f[i] with
        respect to x[j]). It must have the following signature:
        ``jac(x) -> ndarray, shape (m, n)``.
        Default is '2-point'.
        """
        return self.A

    def hess(self, x: np.ndarray) -> np.ndarray:
        """Method for computing the Hessian matrix.
        It must return the Hessian matrix of ``dot(fun, v)`` and
        must have the following signature:
        ``hess(x, v) -> ndarray, shape (n, n)``.
        Here ``v`` is ndarray with shape (m,) containing Lagrange multipliers.
        """
        return np.zeros((self.A.shape[1], self.A.shape[1], self.A.shape[0]))


class ModelOptions(_opt.Options):
    """Class that contains model definition options.

    The individual options can be accessed as object instance attributes
    (e.g. `obj.name`) or as dictionary items (e.g. `obj['name']`) for both
    getting and setting.
    It can be unpacked to give keywords for `SynapseParamModel` `classmethods`.

    Parameters
    ----------
    topology : TopologyOptions, optional keyword
        Topology specifying options. By default `TopologyOptions()`.
    binary : bool, optional keyword
        Is the weight vector binary? Otherwise it's linear. Default: `False`
    rng : np.random.Generator, optional keyword
        Source of random numbers. By default, `builders.RNG`.
    frac : ArrayLike, optional keyword
        Fraction of events of each plasticity type, An extra element is added
        if it is subnormalised. By default `0.5`.
    npl : int, optional keyword
        Total number of plasticity types. By default `topology.npl`
    nst : int or None, optional keyword
        Total number of states. Calculate from `npar` if `None` (default).
    npar : int or None, optional keyword
        Total number of parameters. Calculate from `nst` if `None` (default).

    All parameters are optional keywords. Any dictionary passed as positional
    parameters will be popped for the relevant items. Keyword parameters must
    be valid keys, otherwise a `KeyError` is raised.
    """
    map_attributes: _opt.Attrs = ('topology',)
    prop_attributes: _opt.Attrs = ('frac', 'npl', 'nst', 'npar')
    key_last: _opt.Attrs = ('frac', 'npl', 'nst', 'npar')
    # topology specifying options
    topology: _ma.TopologyOptions
    # Is the weight vector binary?
    binary: bool
    # source of random numbers
    rng: np.random.Generator
    # used by properties
    _frac: la.lnarray
    _nst: Optional[int]
    _npar: Optional[int]

    def __init__(self, *args, **kwds) -> None:
        self.topology = _ma.TopologyOptions()
        self.binary = False
        self.rng = _bld.RNG
        self._frac = la.array([0.5, 0.5])
        self._nst = None
        self._npar = None
        super().__init__(*args, **kwds)

    def check_complete(self) -> bool:
        """Check if model's options have been completely set

        Returns
        -------
        done : bool
            True if `nst` and `npr` both have values.

        Raises
        ------
        TypeError
            If `nst` and `npar` are both `None`
        """
        if None in {self.nst, self.npar}:
            raise TypeError("Must specify one of [nst, npar]")

    @property
    def frac(self) -> la.lnarray:
        """Probability of each plasticity type.
        """
        return self._frac

    @frac.setter
    def frac(self, value: Optional[_sb.ArrayLike]) -> None:
        """Set the probability of each plasticity type.

        Does nothing if `value` is `None`. Adds an element to the end if
        subnormalised.
        """
        if value is None:
            return
        self._frac = _sb.append_frac(value, self.topology.npl)
        self._frac = _sb.trim_frac(self._frac, self.topology.npl)

    @property
    def npl(self) -> int:
        """Number of plasticity types
        """
        return self._frac.shape[-1]

    @npl.setter
    def npl(self, value: Optional[int]) -> None:
        """Set the number of plasticity types.

        Does nothing if `value` is `None`. Removes end elements of `frac`
        if shortening. Appends 0s if lengthening.
        """
        if value is None:
            return
        self.topology.npl = value
        if value <= self.npl:
            self._frac = _sb.trim_frac(self._frac, value)
            return
        extra = la.zeros(self._frac.shape[:-1] + (value - self.npl,))
        self._frac = np.concatenate((self._frac, extra), axis=-1)

    @property
    def nst(self) -> int:
        """Number of states
        """
        return self._nst

    @nst.setter
    def nst(self, value: Optional[int]) -> None:
        """Set the number of states.

        Does nothing if `value` is `None`.
        """
        if value is None:
            return
        self._nst = value
        self._npar = self.npl * _mp.num_param(value, **self.topology)

    @property
    def npar(self) -> int:
        """Number of parameters
        """
        return self._npar

    @npar.setter
    def npar(self, value: Optional[int]) -> None:
        """Set the number of states.

        Does nothing if `value` is `None`.
        """
        if value is None:
            return
        self._npar = value
        self._nst = _mp.num_state(value // self.npl, **self.topology)

# =============================================================================
# Class for specifying optimiser options
# =============================================================================


class ProblemOptions(_opt.Options):
    """Class for specifying problem creator options

    The individual options can be accessed as object instance attributes
    (e.g. `obj.name`) or as dictionary items (e.g. `obj['name']`) for both
    getting and setting.

    Parameters
    ----------
    keep_feasible : bool
        Ensure constraints are kept at intermediate steps? By default `False`.
    inv : bool
        Store inverted fundamental matrices? This will be set `True` by
        other options when needed, so there is no need to set it manually.
    cond : bool
        Include condition number when checking if model is well behaved? This
        will be set to the inverse of `cond_lim`, so there is usually no need
        to set it manually.
    hess : bool
        Use Hessian of loss function/constraints? By default `False`.
        It is automatically set `False` when parent sets `method` to 'SLSQP'.
    cond_lim : bool
        Use the condition number as a constraint? By default `False`.
    cond_thresh : float
        Upper bound on condition number. By default `1e4`.

    All parameters are optional keywords. Any dictionary passed as positional
    parameters will be popped for the relevant items. Keyword parameters must
    be valid keys, otherwise a `KeyError` is raised.
    """
    prop_attributes: _opt.Attrs = ('hess', 'cond_lim')
    key_last: _opt.Attrs = ('hess', 'cond_lim')
    # Ensure constraints are kept at intermediate steps?
    keep_feasible: bool
    # Include condition number when checking if model is well behaved?
    cond: bool
    # Store inverted fundamental matrices?
    inv: bool
    # Use Hessian of loss function/constraints?
    _hess: bool
    # Use the condition number as a constraint?
    _cond_lim: bool
    # Upper bound on condition number
    cond_thresh: float

    def __init__(self, *args, **kwds) -> None:
        self.keep_feasible = False
        self.inv = False
        self.cond = True
        self._hess = False
        self._cond_lim = False
        # One param we don't store here
        self.cond_thresh = 1e4
        super().__init__(*args, **kwds)

    def for_fns(self) -> Dict[str, bool]:
        """keywords for a Maker function"""
        return {'inv': self.inv, 'cond': self.cond, 'svd': self._cond_lim}

    def lin(self) -> Dict[str, bool]:
        """keywords for a `LinearConstraint`"""
        return {'keep_feasible': self.keep_feasible}

    def nlin(self) -> Dict[str, bool]:
        """keywords for a `NonlinearConstraint`"""
        # could add _diff_rel_step=None, finite_diff_jac_sparsity=None
        return {'keep_feasible': self.keep_feasible}

    @property
    def hess(self) -> bool:
        """Do we use the hessian?
        """
        return self._hess

    @hess.setter
    def hess(self, value: Optional[bool]) -> None:
        """Choose whether to use the hessian.

        Does nothing if `value` is `None`.
        """
        if value is None:
            return
        self._hess = value
        if value:
            self.inv = True

    @property
    def cond_lim(self) -> bool:
        """Do we use a condition number constraint?
        """
        return self._cond_lim

    @cond_lim.setter
    def cond_lim(self, value: Optional[bool]) -> None:
        """Choose whether to use a condition number constraint.

        Does nothing if `value` is `None`.
        """
        if value is None:
            return
        self._cond_lim = value
        self.cond = not value


class OptimOptions(_opt.MasterOptions, fallback='extra'):
    """Class for specifying optimiser options

    The individual options can be accessed as object instance attributes
    (e.g. `obj.name`) or as dictionary items (e.g. `obj['name']`) for both
    getting and setting. When subscripting with an unknown key, it will search
    through `model_opts`, `problem_opts` and `extra`. Setting completely
    unknown keys will store them in `extra`.

    Parameters
    ----------
    model_opts : ModelOptions
        Options for constructing models, by default `ModelOptions()`.
    problem_opts : ProblemOptions
        Options for creating problem, by default `ProblemOptions()`.
    extra : Dict[str, Any]
        Extra keyword arguments for `scipy.optimize.minimize`, by default `{}`.
    max_tries : int
        Maximum number of times to attempt optimisation until a good solution
        is found. By default 100
    repeats : int
        Number of times each optimisation is repeated after the first good
        solution. By default 0
    method : str
        Method to use in `scipy.optimize.minimize`, which must be one of
        `'SLSQP'` or `'trust-constr'`. By default `'SLSQP'`.
    maker : Maker
        The function that creates the problem, usually `normal_problem` or
        `shifted_problem`. By default `normal_problem`.

    All parameters are optional keywords. Any dictionary passed as positional
    parameters will be popped for the relevant items.
    """
    map_attributes: _opt.Attrs = ('model_opts', 'problem_opts', 'extra')
    prop_attributes: _opt.Attrs = ('method',)
    key_last: _opt.Attrs = (ModelOptions.key_last + ProblemOptions.key_last
                             + ('method', 'maker'))
    model_opts: ModelOptions
    problem_opts: ProblemOptions
    extra: Dict[str, Any]
    max_tries: int
    repeats: int
    _method: str
    _maker: Maker

    def __init__(self, *args, **kwds) -> None:
        self.model_opts = ModelOptions()
        self.problem_opts = ProblemOptions()
        self.extra = {}
        self.max_tries = 100
        self.repeats = 0
        self._method = "SLSQP"
        self._maker = normal_problem
        super().__init__(*args, **kwds)

    @property
    def shifted(self) -> bool:
        """Is it a shifted problem?

        Returns
        -------
        shifted : bool
            Is problem maker `shifted_problem`?
        """
        return self._maker is shifted_problem

    @shifted.setter
    def shifted(self, value: Optional[bool]):
        """Is it a shifted problem?

        Parameters
        ----------
        value : bool
            Should the problem maker be `shifted_problem`?
        """
        if value is None:
            return
        if value:
            self._maker = shifted_problem
        else:
            self._maker = normal_problem

    def maker(self, model: SynapseOptModel, rate: Optional[Number]
              ) -> Problem:
        """Call the `Maker` function

        Parameters
        ----------
        model : SynapseOptModel
            The model to be used to compute function and gradients.
        rate : Optional[Number]
            Rate parameter to loss fuctions

        Returns
        -------
        problem : Tuple[Func, Func, Func|None, Bounds, List[Constraint]]
            The problem's ingredients: (fun, jac, hess, bounds, lims)
        """
        return self._maker(model, rate, self.problem_opts)

    def set_maker(self, value: Optional[Maker]) -> None:
        """Choose the function that creates the problem.

        Does nothing if `value` is `None`.
        """
        if value is None:
            return
        self._maker = value
        if value is shifted_problem:
            self.problem_opts.hess = False

    @property
    def method(self) -> str:
        """Method to use in `scipy.optimize.minimize`.
        """
        return self._method

    @method.setter
    def method(self, value: Optional[str]) -> None:
        """Choose method to use in `scipy.optimize.minimize`.

        Does nothing if `value` is `None`.
        """
        if value is None:
            return
        if value not in {'SLSQP', 'trust-constr'}:
            raise ValueError('method must be one of "SLSQP", "trust-constr"')
        self._method = value
        if value == 'SLSQP':
            self.problem_opts.inv = True


# =============================================================================
# Class for holding optimisation problems
# =============================================================================
# options, model, fun, jac, hess, bounds, constraints, x_init, rate


class OptimProblem:
    """Class for holding optimisation problems

    Parameters
    ----------
    rate : float
        The value at which we evaluate the loss function, etc.
    params : ArrayLike
        The parameters to initialise the model.
    opts : OptimOptions
        Options for customising the optimiser.

    Attributes
    ----------
    model : SynapseOptModel|None
        The model object used for calculating loss functions etc.
    x_init : la.lnarray|None
        The parameters to initialise the model.
    fns : Tuple[Callable[ndarray], float|ndarray]] * 3
        The loss function, its Jacobian, and its Hessian.
    bounds : Bounds
        Bounds on the model parameters.
    constraints : List[LinearConstraint|NonlinearConstraint|Dict[str, Any]]
        Constraints on the model parameters.
    """
    # The model object used for calculating loss functions etc.
    model: Optional[SynapseOptModel]
    # The parameters to initialise the model
    x_init: Optional[la.lnarray]
    # The loss function, its Jacobian, and its Hessian.
    fns: Tuple[Optional[Func], Optional[Func], Optional[Func]]
    # Bounds on the model parameters
    bounds: Optional[sco.Bounds]
    # Constraints on the model parameters
    constraints: List[Constraint]
    # The value at which we evaluate the loss function, etc.
    _rate: Optional[float]
    # Options for customising the optimiser.
    opts: OptimOptions

    def __init__(self, **kwds) -> None:
        self.model = None
        self.fns = (None, None, None)
        self.bounds, self.constraints = None, []
        self._rate = kwds.pop('rate', None)
        self.x_init = kwds.pop('params', None)
        if self.x_init is not None:
            self.x_init = la.asarray(self.x_init)
        self.opts = kwds.pop('opts', OptimOptions())
        self.opts.update(kwds)
        self._make_model()

    def _make_model(self) -> None:
        """Create a model"""
        if self.x_init is None:
            self.model = SynapseOptModel.rand(**self.opts.model_opts)
            self.x_init = self.model.get_params()
        else:
            self.model = SynapseOptModel.from_params(
                self.x_init, **self.opts.model_opts)

    def _make_problem(self) -> None:
        """Create the loss function, gradient, constraints, etc."""
        if self.model is None:
            self._make_model()
        self.x_init = self.model.get_params()
        prob = self.opts.maker(self.model, self._rate)
        self.fns = prob[:3]
        self.bounds, self.constraints = prob[3:]
        if self.model.topology.constrained:
            self.constraints = []
        elif self.opts.method == 'SLSQP':
            self.constraints = conv_constraints(self.constraints)

    def for_scipy(self) -> Dict[str, Any]:
        """Make an `scipy.optimize.minimize` problem.
        """
        if self.fns[0] is None:
            self._make_problem()
        return {'fun': self.fns[0], 'jac': self.fns[1], 'hess': self.fns[2],
                'bounds': self.bounds, 'constraints': self.constraints,
                'x0': self.x_init, **self.opts.extra}

    def update_init(self, value: Optional[la.lnarray] = None) -> None:
        """Change initial guess"""
        if value is None:
            self.x_init = self.opts.model_opts.rng.random(self.x_init.size)
        else:
            self.x_init = value

    def update_rate(self, value: Optional[float]) -> None:
        """Change the value at which we evaluate the loss function, etc."""
        self.fns = (None, None, None)
        self.bounds, self.constraints = None, []
        if value is not None:
            self._rate = value

    def verify_solution(self, result: sco.OptimizeResult) -> bool:
        """Verify that solution satisfies constraints"""
        # maxcv = constr_violation(prob, result)
        # return maxcv < prob.get('tol', 1e-3)
        if not np.isfinite(result.fun):
            return False
        itol = self.opts.extra.get('tol', 1e-3)
        maxcv = getattr(result, 'constr_violation', None)
        if maxcv is not None:
            return maxcv < itol
        # must be SLSQP
        if _fail_bnd(result.x, self.bounds, itol):
            return False
        for constr in conv_constraints(self.constraints):
            if _fail_cons(result.x, constr, itol):
                return False
        return True


# =============================================================================
# Problem creation helper functions
# =============================================================================


def make_fn(model: SynapseOptModel, method: str, *args,
            opts: Optional[ProblemOptions] = None, **kwds) -> Func:
    """Make a loss function from model, method and value

    Parameters
    ----------
    model : SynapseOptModel
        The object used for computations
    method : str|Callable[[np.ndarray,...], Number]
        The method of `model` that gives us the loss function, or the method.
    opts : ProblemOptions|None, optional
        Object containing options for the whole problem, by default None

    Returns
    -------
    fun : Callable[[np.ndarray], Number]
        The loss function, wrapped so that it only needs model parameters.
    The `model` is also returned as `fun.model`.
    """
    if opts is not None:
        kwds.update(opts.for_fns())

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


def make_fn_set(model: SynapseOptModel, method: str, opts: ProblemOptions,
                *args, **kwds) -> Tuple[Func, ...]:
    """Make a loss function, jacobian, hessian from model, method and value

    Parameters
    ----------
    model : SynapseOptModel
        The object used for computations
    method : str|Callable[[ndarray,...], Number]
        The method of `model` that gives us the loss function, or the method.
    opts : ProblemOptions|None, optional
        Object containing options for the whole problem, by default `None`.

    Returns
    -------
    fun : Callable[[ndarray(M)], Number|ndarray(N)]
        The loss function, wrapped so that it only needs model parameters.
    jac : Callable[[ndarray(M)], ndarray(M)|(M,N)]
        Jacobian of `fun`, wrapped so that it only needs model parameters.
    hess : Callable[[ndarray(M)], ndarray(M,M)|(M,M,N)]
        Hessian of `fun`, wrapped so that it only needs model parameters.
    The `model` is also returned as `fun.model`.
    """
    fun = make_fn(model, method + '_fun', *args, opts=opts, **kwds)
    jac = make_fn(model, method + '_grad', *args, opts=opts, **kwds)
    if opts.hess and hasattr(model, method + '_hess'):
        hess = make_fn(model, method + '_hess', *args, opts=opts, **kwds)
    else:
        hess = None
    return fun, jac, hess


def cond_limit(model: SynapseOptModel, rate: float, opts: ProblemOptions,
               ) -> sco.NonlinearConstraint:
    """Create a constraint on the condition number

    Parameters
    ----------
    model : SynapseOptModel
        The object used for computations
    rate : float
        Value of Laplace transorm parameter.
    opts : ProblemOptions
        Object containing options for the whole problem.

    Returns
    -------
    constraint : sco.NonlinearConstraint
        The constraint.
    """
    cond_fn, cond_jac, _ = make_fn_set(model, 'cond', opts, rate,
                                       cond_thresh=opts.cond_thresh)
    return sco.NonlinearConstraint(cond_fn, 0, np.inf, cond_jac, **opts.nlin())


def constraint_coeff(model: SynapseOptModel) -> la.lnarray:
    """Coefficient matrix for upper bound on off-diagonal row sums.

    Parameters
    ----------
    model: SynapseOptModel
        Model whose shape determines `coeffs`.

    Returns
    -------
    coeffs : la.lnarray (2*nst, 2*nst(nst-1))
        matrix of coefficients s.t ``coeffs @ params <= 1``.
    """
    npl, nst = model.nplast, model.nstate
    return type(model).constraint_coeff(npl, nst)


def conv_constraints(constraints: Sequence[Constraint]) -> List[SLSQPCons]:
    """Convert sequence of constraints from trust-constr to SLSQP format

    Parameters
    ----------
    constraints : Sequence[dict|LinearConstraint|sco.NonlinearConstraint]
        The constraints in their original form

    Returns
    -------
    constraint_dicts : List[dict]
        Equivalent constraints in SLSQP format, in a list.

    Raises
    ------
    TypeError
        If elements of `constraints` are not of one of the listed types.
        Note that it is the version of `LinearConstraint` from this module.
    """
    return _cn.map_join(conv_constraint, constraints)


def conv_constraint(constraint: Constraint) -> List[SLSQPCons]:
    """Convert one constraint from trust-constr to SLSQP format

    Parameters
    ----------
    constraint : dict|LinearConstraint|sco.NonlinearConstraint
        The constraint in its original form

    Returns
    -------
    constraint_dicts : List[dict]
        Equivalent constraints in SLSQP format, in a list. If `constraint.lb`
        and `constraint.ub` are both finite and different, the list will
        contain two constraints.

    Raises
    ------
    TypeError
        If `constraint` is not of one of the listed types. Note that it is the
        version of `LinearConstraint` from this module.
    """
    if not isinstance(constraint,
                      (sco.NonlinearConstraint, LinearConstraint, dict)):
        raise TypeError(f"Unknown constraint type: {type(constraint)}")
    elif isinstance(constraint, dict):
        return [constraint]

    slsqp_lb = {'type': 'ineq', 'args': ()}
    slsqp_ub = slsqp_lb.copy()
    slsqp_lb['fun'] = lambda x: constraint.fun(x) - constraint.lb
    slsqp_ub['fun'] = lambda x: constraint.ub - constraint.fun(x)
    if callable(constraint.jac):
        slsqp_lb['jac'] = constraint.jac
        slsqp_ub['jac'] = lambda x: - constraint.jac(x)

    if not np.isfinite(constraint.lb):
        return [slsqp_ub]
    if not np.isfinite(constraint.ub):
        return [slsqp_lb]
    if constraint.lb == constraint.ub:
        slsqp_lb['type'] = 'eq'
        return [slsqp_lb]
    return [slsqp_lb, slsqp_ub]


# =============================================================================
# Problem creators
# =============================================================================


def normal_problem(model: SynapseOptModel, rate: float,
                   opts: ProblemOptions) -> Problem:
    """Make an optimize problem.

    Parameters
    ----------
    model : SynapseOptModel
        The object used for computations
    rate : float
        Value of Laplace transorm parameter.
    opts : ProblemOptions
        Object containing options for the whole problem.

    Returns
    -------
    fun : Callable[[ndarray(M)], Number|ndarray(N)]
        The loss function, wrapped so that it only needs model parameters.
    jac : Callable[[ndarray(M)], ndarray(M)|(M,N)]
        Jacobian of `fun`, wrapped so that it only needs model parameters.
    hess : Callable[[ndarray(M)], ndarray(M,M)|(M,M,N)]
        Hessian of `fun`, wrapped so that it only needs model parameters.
    bounds : sco.Bounds
        Bounds on each transition rate
    diag : list[LinearConstraint|sco.NonNonlinearConstraint]
        Bound on off-diagonal row sums.
    """
    fun, jac, hess = make_fn_set(model, 'laplace', opts, rate)

    con_coeff = constraint_coeff(model)
    bounds = sco.Bounds(0, 1, **opts.lin())
    diag = [LinearConstraint(con_coeff, -np.inf, 1, **opts.lin())]
    if opts.cond_lim:
        diag.append(cond_limit(model, rate, opts))

    return fun, jac, hess, bounds, diag


# -----------------------------------------------------------------------------
# Shifting rate from function to constraints
# -----------------------------------------------------------------------------


def shifted_problem(model: SynapseOptModel, rate: float,
                    opts: ProblemOptions) -> Problem:
    """Make an optimize problem with rate shifted to constraints.

    Parameters
    ----------
    model : SynapseOptModel
        The object used for computations
    rate : float
        Value of Laplace transorm parameter.
    opts : ProblemOptions
        Object containing options for the whole problem.

    Returns
    -------
    fun : Callable[[ndarray(M)], Number|ndarray(N)]
        The loss function, wrapped so that it only needs model parameters.
    jac : Callable[[ndarray(M)], ndarray(M)|(M,N)]
        Jacobian of `fun`, wrapped so that it only needs model parameters.
    hess : Callable[[ndarray(M)], ndarray(M,M)|(M,M,N)]
        Hessian of `fun`, wrapped so that it only needs model parameters.
    bounds : sco.Bounds|None
        Bounds on each transition rate.
    diag : list[sco.NonlinearConstraint]
        Constraints on transition rates.
    """
    if model.topology.constrained:
        raise ValueError('Shifted problem cannot have special topology')

    fun, jac, hess = make_fn_set(model, 'area', opts, rate)

    bounds = None
    cfun, cjac, chess = make_fn_set(model, 'peq_min', opts, rate)
    lims = sco.NonlinearConstraint(cfun, 0, np.inf, cjac, chess, **opts.nlin())
    if opts.cond_lim:
        lims = [lims, cond_limit(model, None, opts)]

    return fun, jac, hess, bounds, _cn.listify(lims)


# =============================================================================
# Karush-Kuhn-Tucker multipliers
# =============================================================================


def _grad_to_kkt(grad: np.ndarray, model: _sb.SynapseParam) -> np.ndarray:
    """Fix diagonals of kkt from gradient wrt off-diagonal elements.

    Parameters
    ----------
    grad : np.ndarray (...,2n(n-1))
        Gradient of loss function wrt off-diagonal elements.
    model : SynapseParam (2,n,n)
        Model with transition rates.

    Returns
    -------
    kkt : np.ndarray (...,2,n,n)
        grad with diagonals `-rowsum(plast*grad)`, subtracted from row.
    """
    kkt = _ma.params.params_to_mat(grad, **model.topology.directed())
    kkt[..., np.diag_indices(model.nstate, 2)] = 0
    kkt_diag = -(model.plast * kkt).sum(axis=-1, keepdims=True)
    kkt += kkt_diag
    return kkt


def _all_kkt_normal(prob: OptimProblem) -> np.ndarray:
    """Calculate Karush-Kuhn-Tucker multiplier for every element.

    Parameters
    ----------
    prob : OptimProblem
        The problem.

    Returns
    -------
    kkt : np.ndarray (2,n,n)
        KKT multipliers.
    """
    model = prob.model
    grad = prob.fns[1](model.get_params())
    return _grad_to_kkt(grad, model)


def _all_kkt_shifted(prob: OptimProblem) -> np.ndarray:
    """Calculate Karush-Kuhn-Tucker multiplier for every element.

    Parameters
    ----------
    prob : OptimProblem
        The problem.

    Returns
    -------
    kkt : np.ndarray (2,n,n)
        KKT multipliers.
    """
    model = prob.model
    constr = conv_constraint(prob.constraints[0])[0]
    # (2,n,n) -> (2*n**2,)
    rhs = _all_kkt_normal(prob).ravel()
    # (2n**2,2n(n-1))
    cjac = constr['jac'](model.get_params())
    # (2n**2,2,n,n) -> (2n**2,2n**2)
    cfs = _grad_to_kkt(cjac, model).ravelaxes(-3).T
    return la.solve(cfs, rhs).reshape(model.nplast, model.nstate, model.nstate)


def kkt_multipliers(prob: OptimProblem) -> np.ndarray:
    """Calculate Karush-Kuhn-Tucker multiplier for every element.

    Parameters
    ----------
    prob : OptimProblem
        The solved problem.

    Returns
    -------
    kkt : np.ndarray (2,n,n)
        KKT multipliers
    """
    if prob.opts.shifted:
        return _all_kkt_shifted(prob)
    return _all_kkt_normal(prob)



# =============================================================================
# Optimisation helper functions
# =============================================================================


def _fail_bnd(sol: np.ndarray, bnds: Optional[sco.Bounds], tol: float) -> bool:
    """Check if solution fails SLSQP constraint"""
    if bnds is None:
        return False
    return (sol < bnds.lb - tol).any() or (sol > bnds.ub + tol).any()


def _fail_cons(soln: np.ndarray, cons: SLSQPCons, itol: float) -> bool:
    """Check if solution fails SLSQP constraint"""
    vals = cons['fun'](soln)
    if cons['type'] == 'eq':
        return not np.allclose(0, vals)
    return (vals < -itol).any()

# =============================================================================
# Optimisation
# =============================================================================


def first_good(prob: OptimProblem) -> sco.OptimizeResult:
    """First solution that satisfies constraints

    Parameters
    ----------
    prob : OptimProblem
        Object containing the whole problem. Updated in place.

    Returns
    -------
    res : sco.OptimizeResult
        Object containing the results of optimisation.
    """
    res = sco.OptimizeResult()
    for _ in _it.dcount('tries', prob.opts.max_tries):
        try:
            res = sco.minimize(**prob.for_scipy())
        except np.linalg.LinAlgError:
            pass
        else:
            if prob.verify_solution(res):
                break
        prob.update_init()
    return res


def optim_laplace(rate: float, prob: Optional[OptimProblem] = None,
                  **kwds) -> sco.OptimizeResult:
    """Optimised model at one value of rate

    Parameters
    ----------
    rate : float
        Value of Laplace transorm parameter.
    prob : OptimProblem|None, optional
        Object containing the (previous) whole problem. If `None`, we build a
        new one using `rate` and extra keyword arguments. By default `None`.

    Returns
    -------
    res : sco.OptimizeResult
        Object containing the results of optimisation.
    """
    if prob is None:
        prob = OptimProblem(rate=rate, **kwds)
    else:
        prob.update_rate(rate)
    res = first_good(prob)
    if not prob.verify_solution(res):
        res.fun = 0
    for _ in _it.dcount('repeats', prob.opts.repeats, disp_step=1):
        prob.update_init()
        new_res = sco.minimize(**prob.for_scipy())
        if prob.verify_solution(new_res) and new_res.fun < res.fun:
            res = new_res
    if not prob.verify_solution(res):
        res.fun = np.nan
    return res


def optim_laplace_range(rates: np.ndarray, nst: int,
                        **kwds) -> Tuple[la.lnarray, la.lnarray]:
    """Optimised model at many values of rate

    Parameters
    ----------
    rates : np.ndarray (T,)
        Parameters of Laplace transform at which we maximise.
    nst : int
        Number of states

    Returns
    -------
    snr : la.lnarray (T,)
        Maximum snr (Laplace transformed) at each value of `rates`.
    models : la.lnarray (T,Q)
        Parameters of models that achieve `snr`.
    """
    prob = OptimProblem(nst=nst, **kwds)
    # model = make_model(kwds, nst=nst)
    snr = la.empty_like(rates)
    models = la.empty((len(rates), prob.opts.model_opts.npar))
    with _it.delay_warnings():
        for i, rate in _it.denumerate('rate', rates):
            res = optim_laplace(rate, prob)
            snr[i] = - res.fun
            models[i] = res.x
    return snr, models


def reoptim_laplace_range(inds: np.ndarray, rates: np.ndarray,
                          snr: np.ndarray, models: np.ndarray,
                          **kwds) -> Tuple[la.lnarray, la.lnarray]:
    """Reoptimised model at many values of rate

    Parameters
    ----------
    inds : np.ndarray[int]
        Indices of models to reoptimise.
    rates : np.ndarray
        Parameters of Laplace transform at which we maximise.
    snr : np.ndarray
        Old maximum snr (Laplace transformed) at each value of `rates`.
    models : np.ndarray
        Old parameters of models that achieve `snr`.

    Returns
    -------
    snr : la.lnarray (T,)
        New maximum snr (Laplace transformed) at each value of `rates`.
    models : la.lnarray (T,Q)
        New parameters of models that achieve `snr`.
    """
    prob = OptimProblem(params=models[inds[0]], **kwds)
    # model = make_model(kwds, params=models[inds[0]])
    with _it.delay_warnings():
        for ind, rate in _it.dzip('rate', inds, rates[inds]):
            res = optim_laplace(rate, prob)
            if - res.fun > snr[ind]:
                snr[ind] = - res.fun
                models[ind] = res.x
    return snr, models


def check_cond_range(rates: np.ndarray, models: np.ndarray,
                     **kwds) -> la.lnarray:
    """Condition numbers of optimised models

    Parameters
    ----------
    rates : np.ndarray (S,)
        inverse time (Laplace parameter)
    models : np.ndarray (S, P)
        optimised models

    Returns
    -------
    cond : np.ndarray (S,)
        condition number, worst of :math:`Z(0), Z(s)`
    """
    cnd = la.empty_like(rates)
    prob = OptimProblem(params=models[0], **kwds)
    model = prob.model
    for i, rate, params in _it.zenumerate(rates, models):
        model.set_params(params)
        cnd[i] = model.cond(rate, rate_max=True)
    return cnd


# =============================================================================
# Theory
# =============================================================================


def proven_envelope_laplace(rate: Data, nst: int) -> Data:
    """Theoretical envelope for Laplace transform.

    Parameters
    ----------
    rates : Number|ndarray
        Rate parameter of Laplace transform.
    nst : int
        Number of states.

    Returns
    -------
    envelope : Data
        Putative maximum A(s) for all models.
    """
    return (nst - 1) / (rate * (nst - 1) + 1)


def equlibrium_envelope_laplace(rate: Data, nst: int, gam: float = 2.) -> Data:
    """Theoretical envelope for Laplace transform assuming detailed balance.

    Parameters
    ----------
    rates : Number|ndarray
        Rate parameter of Laplace transform.
    nst : int
        Number of states.
    gam : float, optional
        Value to use in the square root bound, by default 2.

    Returns
    -------
    envelope : Data
        Putative maximum A(s) for all models.
    """
    if nst < 3:
        return (nst - 1) / (rate * (nst - 1) + 1)
    s_two, s_sticky, i_a = 1./gam, gam/(nst-1)**2, gam/(nst-1)
    env = 0.5 / np.sqrt(s_two * rate)
    env = np.where(rate > s_two, 1 / (rate + s_two), env)
    env = np.where(rate < s_sticky, i_a / (rate + s_sticky), env)
    return env


def heuristic_envelope_laplace(rate: Data, nst: int) -> Data:
    """Heuristic envelope for Laplace transform.

    Parameters
    ----------
    rate : Number|ndarray
        Rate parameter of Laplace transform.
    nst : int
        Number of states.

    Returns
    -------
    envelope : Data
        Putative maximum A(s) for all models.
    """
    rate = la.array(rate)
    s_two = rate >= _sh.s_star(2)
    s_sticky = rate < _sh.s_star(nst)
    s_short = np.logical_not(np.logical_or(s_two, s_sticky))
    env_two = _sh.short_star_s(rate[s_two], 2)
    env_short = _sh.uni_star_s(rate[s_short])
    env_sticky = _st.sticky_star_s(rate[s_sticky], nst)
    return np.concatenate((env_sticky, env_short, env_two))


# =============================================================================
# Type hints
# =============================================================================
Data = TypeVar('Data', Number, np.ndarray, float, covariant=True)
Func = Callable[[np.ndarray], Data]
SLSQPCons = Dict[str, Any]
Constraint = Union[LinearConstraint, sco.NonlinearConstraint, SLSQPCons]
Problem = Tuple[Func[float], Func[np.ndarray], Func[np.ndarray],
                Optional[sco.Bounds], List[Constraint]]
Maker = Callable[[SynapseOptModel, Optional[float], ProblemOptions], Problem]
