# -*- coding: utf-8 -*-
"""Testing the square root bound
"""
from __future__ import annotations

from typing import Tuple, List, Optional, Type, cast

import numpy as np
import scipy.optimize as sco
import scipy.linalg as sla

import numpy_linalg as la
import sl_py_tools.iter_tricks as _it
import sl_py_tools.numpy_tricks.markov as _ma
import sl_py_tools.numpy_tricks.markov.params as _mp

import complex_synapse.builders as bld
import complex_synapse.optimise as cso
import complex_synapse.synapse_base as sb
import complex_synapse.synapse_mem as sm
from complex_synapse.optimise.optimise import Func, Constraint, Problem

# =============================================================================


def check_sqrt(model: sm.SynapseModel) -> float:
    """Check if sqrt bound is met

    Parameters
    ----------
    model : SynapseModel
        Synapse model object for computation

    Returns
    -------
    sqrt : float
        max_a I_a sqrt(tau_a)
    """
    taus, inits = model.spectrum()
    assert (taus.imag**2 < 1e-6).all()
    assert (inits.imag**2 < 1e-6).all()
    return (inits.real * np.sqrt(taus.real / 2)).max()[()]


def rand_sqrt(nst: int, **kwds) -> Tuple[float, SynapseSqrt]:
    """Test sqrt bound with a random model.

    Parameters
    ----------
    nst : int
        Number of states.

    Returns
    -------
    sqrt : float
        max_a I_a sqrt(tau_a)
    model : SynapseSqrt
        The random model.
    """
    kwds.setdefault('binary', True)
    model = SynapseSqrt.rand(nst, npl=2, **kwds)
    model.symmetrise()
    return check_sqrt(model), model


def rand_balanced(nst: int, frac: float = 0.5) -> la.lnarray:
    """A random transition matrix with detailed balance

    Parameters
    ----------
    nst : int
        Number of states, M.
    frac : float, optional
        Fraction of events that are potentiating, by default 0.5

    Returns
    -------
    trans_mat : la.lnarray (M,M)
        Transition matrix.
    """
    flux = bld.RNG.random(nst**2)
    flux /= flux.sum()
    return ptom_balanced(flux[:-1], frac)


def ptom_balanced(param: la.lnarray, frac: float = 0.5) -> la.lnarray:
    """A transition matrix with detailed balance from its parameters.

    Parameters
    ----------
    param : la.lnarray (M^2 - 1,)
        Parameters: fluxes, except last element.
    frac : float, optional
        Fraction of events that are potentiating, by default 0.5

    Returns
    -------
    mat : la.lnarray (M,M)
        Transition matrix.
    """
    nst = int(np.sqrt(param.size + 1))
    frac = sb.append_frac(frac, 2)
    # coeffs, lbs = constraint(nst, frac)
    # print(coeffs @ param - lbs)
    param = la.r_[param.copy(), 1 - param.sum()].reshape((nst, nst))
    peq, pneq = param.sum(-1), param.sum(-2)
    fpot = (param - np.diagflat(peq)) * frac[0]
    fdep = (param.t - np.diagflat(pneq)) * frac[0]
    mrk = np.stack((fpot, fdep)) / peq.c / frac.s
    return mrk


def mtop_balanced(mat: la.lnarray, frac: float = 0.5,
                  peq: Optional[la.lnarray] = None) -> la.lnarray:
    """Parameters of a transition matrix with detailed balance.

    Parameters
    ----------
    mat : la.lnarray (M,M)
        Transition matrix.
    frac : float, optional
        Fraction of events that are potentiating, by default 0.5
    peq : lnarray (M,)|None, optional
        Steady-state distribution, calculate from mat if None (default).

    Returns
    -------
    param : la.lnarray (M^2 - 1,)
        Parameters: fluxes, except last element.
    """
    if peq is None:
        mrk = frac * mat[0] + (1-frac) * mat[1]
        peq = _ma.calc_peq(mrk)
    mat += np.eye(mat.shape[-1])
    flux = peq.c * mat
    param = flux[0].ravel()
    return param[:-1]


def _row_cnst(row: int, nst: int, frac: float = 0.5) -> la.lnarray:
    """Constraint for one diagonal of M^dep to be positive

    Parameters
    ----------
    row : int
        Row whose sum we compute.
    nst : int
        Number of states, M.
    frac : float, optional
        Fraction of events that are potentiating, by default 0.5

    Returns
    -------
    coeffs : lnarray (M^2 - 1,)
        Coefficients of parameters for constraint.
    """
    fpp, fmm = sb.append_frac(frac, 2)
    coeff = la.zeros((nst, nst))
    coeff[:, row] = -fpp
    coeff[row, :] = fmm
    if row == nst - 1:
        coeff -= fmm
    return coeff.ravel()[:-1]


def constraint(nst: int, frac: float = 0.5):
    """Constraint to enforce non-negative diagonals.

    Parameters
    ----------
    nst : int
        Number of states.
    frac : float, optional
        Fraction of events that are potentiating, by default 0.5

    Returns
    -------
    coeffs : lnarray (nst+1,nst**2-1)
        Coefficients of linear constraints on parameters.
    lbs : lnarray (nst+1,)
        Lower bounds on `coeffs @ params`.
    """
    frac = sb.append_frac(frac, 2)
    rows = [_row_cnst(i, nst, frac) for i in range(nst)]
    rows.append(-la.ones(nst**2 - 1))
    coeffs = np.stack(rows)
    lbs = la.zeros(nst+1)
    lbs[-2] = -frac[1]
    lbs[-1] = -1
    return coeffs, lbs


def valid(params: la.lnarray, frac: float = 0.5) -> bool:
    """Are these parameter values valid?

    Parameters
    ----------
    param : la.lnarray (M^2 - 1,)
        Parameters: fluxes, except last element.
    frac : float, optional
        Fraction of events that are potentiating, by default 0.5

    Returns
    -------
    is_valid : bool
        Are the parameters and their sums in range?
    """
    if (params < 1e-5).any():
        return False
    nst = int(np.sqrt(params.size + 1))
    coeffs, lbs = constraint(nst, frac)
    return (coeffs @ params > lbs).all()


def get_valid(nst: int, frac: float = 0.5, rng: np.random.Generator = bld.RNG):
    """Get a valid set of parameter values

    Parameters
    ----------
    nst : int
        Number of states, M.
    frac : float, optional
        Fraction of events that are potentiating, by default 0.5
    rng : np.random.Generator, optional
        Our random number generator, by default la.random.default_rng().

    Returns
    -------
    param : la.lnarray (M^2 - 1,)
        Parameters: fluxes, except last element.
    """
    params = rng.random(nst**2 - 1) * 2 / nst**2
    for _ in _it.undcount('tries', 100, disp_step=10):
        if valid(params, frac):
            return params
        params = rng.random(nst**2 - 1) * 2 / nst**2
    return params


def get_serial(nst: int, frac: float = 0.5, rng: np.random.Generator = bld.RNG
               ) -> la.lnarray:
    """Get a valid set of parameter values for a serial model

    Parameters
    ----------
    nst : int
        Number of states, M.
    frac : float, optional
        Fraction of events that are potentiating, by default 0.5
    rng : np.random.Generator, optional
        Our random number generator, by default la.random.default_rng().

    Returns
    -------
    param : la.lnarray (2*M - 2,)
        Parameters: fluxes, except last element.
    """
    serial = np.diagflat(la.full(nst-1, 1/nst), 1).ravel()[:-1]
    params = (rng.random(nst**2 - 1) - 0.5) * 2 / nst**2
    for _ in _it.undcount('tries', 100, disp_step=10):
        if valid(params + serial, frac):
            return params + serial
        params = (rng.random(nst**2 - 1) - 0.5) * 2 / nst**2
    return serial


# =============================================================================


class SynapseSqrt(sb.SynapseParam, sm.SynapseModel):
    """class for testing sqrt bound"""
    _param: la.lnarray
    _peq: la.lnarray
    _pneq: la.lnarray
    _flux: la.lnarray
    _ind: Optional[int]
    _changed: bool
    _saved: Tuple[la.lnarray, ...]


    def __init__(self, *args, **kwds) -> None:
        super().__init__(*args, **kwds)
        self._param = la.zeros((self.nstate, self.nstate))
        self._flux = la.zeros((self.nstate, self.nstate))
        self._peq = la.zeros(self.nstate)
        self._pneq = la.zeros(self.nstate)
        self._ind = None
        self._changed = True
        self._saved = (la.array(0), la.array(0), la.array(0))

    def biggest(self) -> float:
        """max(I_a * sqrt(t_a / 2)), the quamtity we optimise."""
        # taus, inits = self.spectrum()
        # if (taus.imag**2 > 1e-6).any() or (inits.imag**2 > 1e-6).any():
        #     self.symmetrise()
        qas, init_ab = self._eig()[:2]
        sqrts = np.diagonal(init_ab) / np.sqrt(2 * qas)
        self._ind = sqrts.argmax()[()]
        return - sqrts[self._ind]

    def biggest_grad(self) -> la.lnarray:
        """gradient of max(I_a * sqrt(t_a / 2))"""
        if self._ind is None:
            self.biggest()
        ind = self._ind
        qas, init_ab, dq_ab, di_aa = self._eig()
        i_by_q = init_ab / (qas.c - qas)
        i_by_q[np.diag_indices(qas.size)] = np.diagonal(init_ab) / qas
        grad = dq_ab[..., ind] @ i_by_q[ind, :]
        grad -= dq_ab[..., ind, :] @ i_by_q[:, ind]
        grad -= dq_ab[..., ind, ind] * i_by_q[ind, ind] / 2
        grad += di_aa[..., ind]
        grad = grad.ravel()
        return - (grad[:-1] - grad[-1]) / np.sqrt(2 * qas[ind])

    def _eig(self) -> Tuple[la.lnarray, ...]:
        """Do eigen-decomp, etc"""
        if not self._changed:
            return self._saved
        self._changed = False
        self._ind = None
        frac = self.frac[0]
        evals, evc = sla.eigh(-self._flux,
                              np.diagflat(self._peq))
        evals, evc = la.asarray(evals), la.asarray(evc)
        uda_w = 2 * frac * evc.T @ (self.weight * self._peq)
        init_ab = uda_w.c * ((self._pneq - self._peq) @ evc)
        dq_mn_ab = evc[:, None, :, None] * evc[:, None] * frac
        dq_mn_ab += evc[:, :, None] * evc[:, None, None] * frac
        dq_mn_ab += evc[:, None, :, None] * evc[:, None, None] * (evals - frac)
        dq_mn_ab -= evc[:, :, None] * evc[:, None] * frac
        di_mn_aa = uda_w * (evc - evc[:, None])
        self._saved = (evals, init_ab, dq_mn_ab, di_mn_aa)
        return self._saved

    def symmetrise(self) -> None:
        """Symmetrise wrt time reversal.
        """
        self += self.time_rev()
        self /= 2

    def set_params(self, params: np.ndarray, *args, **kwds) -> None:
        # inherit docstring
        # super().set_params(params, *args, **kwds)
        # self.symmetrise()
        if np.allclose(params, self.get_params()):
            return
        self._changed = True
        nst = int(np.sqrt(params.size + 1))
        self._param = la.r_[params, 1 - params.sum()].reshape((nst, nst))
        self._peq, self._pneq = self._param.sum(-1), self._param.sum(-2)
        fpot = (self._param - np.diagflat(self._peq)) * self.frac[0]
        fdep = (self._param.t - np.diagflat(self._pneq)) * self.frac[0]
        factor = self._peq.c * self.frac.s
        self.plast = np.stack((fpot, fdep)) / factor
        self._flux = fpot + fdep
        # self.plast = ptom_balanced(params, self.frac[0])

    def get_params(self, ravel: bool = True) -> la.lnarray:
        # inherit docstring
        # return mtop_balanced(self.plast, self.frac[0], self.peq())
        return self._param.ravel()[:-1]

    def current(self) -> la.lnarray:
        """Antisymmetric part of flux"""
        if self.topology.serial:
            return 0
        if self.topology.ring:
            param = self.get_params(ravel=False)
            return np.diff(param.prod(axis=-1))
        flux = self.peq().c * self.markov()
        curr = flux - flux.t
        return _mp.mat_to_params(curr, **self.directed(None))

    @classmethod
    def rand(cls: Type[sb.Syn], nst: int, *args, **kwargs) -> sb.Syn:
        # inherit docstring
        kwargs['npar'] = nst**2 - 1
        kwargs['nst'] = nst
        serial = kwargs.pop('serialish', False)
        rng = kwargs.pop('rng', bld.RNG)
        frac = kwargs.get('frac', 0.5)
        if serial:
            params = get_serial(nst, frac, rng)
        else:
            params = get_valid(nst, frac, rng)
        return super().from_params(params, *args, **kwargs)
        # return super().rand(nst, *args, **kwargs)

    @staticmethod
    def constraint_coeff(npl: int, nst: int) -> la.lnarray:
        """Coefficient matrix for upper bound on off-diagonal row sums.

        Parameters
        ----------
        npl : int
            Number of types of plasticity
        nst: int
            Number of states.

        Returns
        -------
        coeffs : la.lnarray (nst**2-1,)
            matrix of coefficients s.t ``coeffs @ params <= 1``.
        """
        assert npl == 2
        return la.ones(nst**2 - 1)


# =============================================================================


class OptimProblem(cso.OptimProblem):
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
    model: Optional[SynapseSqrt]
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
    opts: cso.OptimOptions

    def __init__(self, **kwds) -> None:
        kwds.setdefault('npl', 2)
        kwds.setdefault('rate', None)
        kwds.setdefault('model_cls', SynapseSqrt)
        kwds.setdefault('maker', normal_problem)
        super().__init__(**kwds)

    # def _make_model(self) -> None:
    #     """Create a model object"""
    #     if self.x_init is None:
    #         self.model = SynapseSqrt.rand(**self.opts.model_opts)
    #     else:
    #         self.model = SynapseSqrt.from_params(self.x_init,
    #                                              **self.opts.model_opts)

    def update_init(self, value: Optional[la.lnarray] = None,
                    serialish: bool = False) -> None:
        """Change initial guess (random)

        Parameters
        ----------
        value : la.lnarray|None, optional
            New parameters for the initial guess, by default None
        serialish : bool, optional
            Do we initialise favouring serial transitions? By default False
        """
        if value is None:
            if serialish:
                self.x_init = get_serial(self.model.nstate, self.model.frac[0],
                                         self.opts.model_opts.rng)
            else:
                self.x_init = get_valid(self.model.nstate, self.model.frac[0],
                                        self.opts.model_opts.rng)
        else:
            self.x_init = value


# =============================================================================


def balance_cons(model: SynapseSqrt, opts: cso.ProblemOptions,
                 ) -> sco.NonlinearConstraint:
    """Create a constraint for detailed balance.

    Parameters
    ----------
    model : SynapseSqrt
        The object used for computations
    opts : ProblemOptions|None, optional
        Object containing options for the whole problem, by default None

    Returns
    -------
    constraint : NonlinearConstraint
        The constraint object in trust-constr format.
    """
    cond_fn = cso.optimise.make_fn(model, 'current')
    return sco.NonlinearConstraint(cond_fn, 0, 0, **opts.nlin())


def normal_problem(model: SynapseSqrt, rate: Optional[float],
                   opts: cso.ProblemOptions) -> Problem:
    """Make an optimize problem.

    Parameters
    ----------
    model : SynapseOptModel
        The object used for computations
    rate : float|None
        Value of Laplace transorm parameter. Should be None.
    opts : ProblemOptions
        Object containing options for the whole problem.

    Returns
    -------
    fun : Callable[[ndarray(M)], Number]
        The loss function, wrapped so that it only needs model parameters.
    jac : Callable[[ndarray(M)], ndarray(M)]
        Jacobian of `fun`, wrapped so that it only needs model parameters.
    hess : Callable[[ndarray(M)], ndarray(M,M)]
        Hessian of `fun`, wrapped so that it only needs model parameters.
    bounds : sco.Bounds
        Bounds on each transition rate
    diag : list[LinearConstraint|sco.NonNonlinearConstraint]
        Bound on off-diagonal row sums.
    """
    assert rate is None
    fun = cso.optimise.make_fn(model, 'biggest')
    jac = cso.optimise.make_fn(model, 'biggest_grad')
    hess = None

    con_coeff, lbs = constraint(model.nstate, model.frac[0])
    # con_coeff = la.ones(model.nstate**2 - 1)
    bounds = sco.Bounds(1e-5, np.inf, **opts.lin())
    diag = [cso.optimise.LinearConstraint(con_coeff, lbs, np.inf, **opts.lin())]
    # diag.append(balance_cons(model, opts))
    if opts.cond_lim:
        diag.append(cso.optimise.cond_limit(model, None, opts))

    return fun, jac, hess, bounds, diag


def optim_sqrt(nst: int, prob: Optional[OptimProblem] = None,
               **kwds) -> sco.OptimizeResult:
    """Optimised model

    Parameters
    ----------
    nst : int
        Number of states
    prob : OptimProblem|None, optional
        Object containing the (previous) whole problem. If `None`, we build a
        new one using extra keyword arguments. By default `None`.

    Returns
    -------
    res : sco.OptimizeResult
        Object containing the results of optimisation.
    """
    kwds.setdefault('repeats', 10)
    kwds.setdefault('nst', nst)
    serialish = kwds.pop('serialish', False)
    step = kwds.pop('disp_step', 1)
    if prob is None:
        prob = OptimProblem(**kwds)
    res = cso.optimise.first_good(prob)
    if not prob.verify_solution(res):
        res.fun = 0
    for _ in _it.dcount('repeats', prob.opts.repeats, disp_step=step):
        prob.update_init(serialish=serialish)
        try:
            new_res = sco.minimize(**prob.for_scipy())
        except np.linalg.LinAlgError:
            pass
        else:
            if prob.verify_solution(new_res) and new_res.fun < res.fun:
                res = new_res
    # if not prob.verify_solution(res):
    #     res.fun = np.nan
    return res
