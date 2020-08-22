"""Fitting synapse models to experiment
"""
from . import baum_welch, fit_synapse, synapse_id
from .synapse_id import SynapseIdModel
from .plast_seq import PlasticitySequence, SimPlasticitySequence
from .fit_synapse import SynapseFitter, GroundedFitter, SynapseFitOptions
from .baum_welch import BaumWelchFitter, GroundedBWFitter
from .video import FitterVideo, VideoOptions
from . import video, plast_seq

assert any((True, baum_welch, fit_synapse, video, synapse_id, plast_seq,
            SynapseIdModel, PlasticitySequence, SimPlasticitySequence,
            SynapseFitter, GroundedFitter, BaumWelchFitter, GroundedBWFitter,
            SynapseFitOptions, FitterVideo, VideoOptions))
