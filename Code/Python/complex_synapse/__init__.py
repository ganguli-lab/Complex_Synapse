# -*- coding: utf-8 -*-
"""
Code for studying complex synapses
"""
from . import builders, optimise, sticky, shorten, identify
from .synapse_mem import SynapseMemoryModel
from .synapse_opt import SynapseOptModel
assert any((True, builders, SynapseMemoryModel, SynapseOptModel, optimise,
            sticky, shorten, identify))
