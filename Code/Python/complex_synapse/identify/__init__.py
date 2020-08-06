"""Fitting synapse models to experiment
"""
from .synapse_id import SynapseIdModel
from .plast_seq import PlasticitySequence, SimPlasticitySequence
from . import baum_welch

assert any((True, SynapseIdModel, PlasticitySequence, SimPlasticitySequence,
            baum_welch))
