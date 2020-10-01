# -*- coding: utf-8 -*-
"""
Code for studying complex synapses
"""
from . import builders, examples, identify, optimise, synapse_base, synapse_mem
from .synapse_mem import SynapseMemoryModel

assert any((True, SynapseMemoryModel, builders, optimise, identify, examples,
            synapse_mem, synapse_base))
