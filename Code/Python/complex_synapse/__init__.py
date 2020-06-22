# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from . import builders
from .synapse_memory_model import SynapseMemoryModel
from .synapse_opt import SynapseOptModel
from . import optimise
from . import soptimise
from .optimise import make_loss_function, make_laplace_problem
assert any((builders, SynapseMemoryModel, SynapseOptModel))
assert any((make_loss_function, make_laplace_problem, optimise, soptimise))
