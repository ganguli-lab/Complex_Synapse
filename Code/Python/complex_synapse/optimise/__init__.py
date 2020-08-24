# -*- coding: utf-8 -*-
"""
Code for studying complex synapses
"""
from . import optimise, sticky, shorten
from .synapse_opt import SynapseOptModel
from .optimise import (make_problem, normal_problem, shifted_problem,
                       verify_solution, update_laplace_problem, optim_laplace,
                       optim_laplace_range, reoptim_laplace_range,
                       check_cond_range, proven_envelope_laplace,
                       heuristic_envelope_laplace, make_model)

assert any((True, SynapseOptModel, optimise, sticky, shorten))
assert any((True, make_model, make_problem, normal_problem, shifted_problem,
            verify_solution, update_laplace_problem, optim_laplace,
            optim_laplace_range, reoptim_laplace_range, check_cond_range,
            proven_envelope_laplace, heuristic_envelope_laplace))
