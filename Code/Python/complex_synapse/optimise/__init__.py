# -*- coding: utf-8 -*-
"""
Code for studying complex synapses
"""
from . import optimise, sticky, shorten, synapse_opt
from .synapse_opt import SynapseOptModel
assert any((True, SynapseOptModel, optimise, sticky, shorten, synapse_opt))
