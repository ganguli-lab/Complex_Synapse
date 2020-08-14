"""Fitting synapse models to experiment
"""
from .synapse_id import SynapseIdModel
from .plast_seq import PlasticitySequence, SimPlasticitySequence
from .fit_synapse import SynapseFitter, GroundedFitter
from . import baum_welch
from .baum_welch import BaumWelchFitter, GroundedBWFitter
from .video import FitterVideo

assert any((True, SynapseIdModel, PlasticitySequence, SimPlasticitySequence,
            baum_welch, SynapseFitter, GroundedFitter, BaumWelchFitter,
            GroundedBWFitter, FitterVideo))
