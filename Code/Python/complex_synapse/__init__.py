# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from . import builders
from . import markov
from .synapse_memory_model import SynapseMemoryModel
from .synapse_opt import SynapseOptModel
assert any((builders, markov, SynapseMemoryModel, SynapseOptModel))
