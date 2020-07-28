"""Fitting synapse models to experiment
"""
from .synapse_id import SynapseIdModel
from .plast_seq import PlasticitySequence, SimPlasticitySequence

assert any((True, SynapseIdModel, PlasticitySequence, SimPlasticitySequence))
