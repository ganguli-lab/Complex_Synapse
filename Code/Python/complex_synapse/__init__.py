# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from . import builders
from .synapse_memory_model import SynapseMemoryModel
from .synapse_opt import SynapseOptModel
from . import optimise
assert any((True, builders, SynapseMemoryModel, SynapseOptModel, optimise))
