"""Fitting synapse models to experiment
"""
from .synapse_id import SynapseIdModel
from .plast_seq import PlasticitySequence, SimPlasticitySequence
from . import baum_welch
from .fit_synapse import SynapseFitter

assert any((True, SynapseIdModel, PlasticitySequence, SimPlasticitySequence,
            baum_welch, SynapseFitter))
