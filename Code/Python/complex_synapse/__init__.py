# -*- coding: utf-8 -*-
"""
Code for studying complex synapses
"""
from . import builders, optimise, identify, synapse_mem, synapse_base, graphs
from .synapse_mem import SynapseMemoryModel
assert any((True, SynapseMemoryModel, builders, optimise, identify,
            synapse_mem, synapse_base, graphs))
