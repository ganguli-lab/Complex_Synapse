"""Fitting synapse models to experiment
"""
from .synapse_id import SynapseIdModel
from .plast_seq import PlasticitySequence, SimPlasticitySequence
from .fit_synapse import SynapseFitter, GroundedFitter
from . import baum_welch, fit_synapse, video, synapse_id, plast_seq
from .baum_welch import BaumWelchFitter, GroundedBWFitter
from .video import FitterVideo

assert any((True, baum_welch, fit_synapse, video, synapse_id, plast_seq,
            SynapseIdModel, PlasticitySequence, SimPlasticitySequence,
            SynapseFitter, GroundedFitter, BaumWelchFitter, GroundedBWFitter,
            FitterVideo))
