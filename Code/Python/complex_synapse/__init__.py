# -*- coding: utf-8 -*-
"""
Code for studying complex synapses
"""
from . import builders, identify, optimise, synapse_base, synapse_mem, sqrt
from .synapse_mem import SynapseMemoryModel

assert any((True, SynapseMemoryModel, builders, optimise, identify, sqrt,
            synapse_mem, synapse_base))
