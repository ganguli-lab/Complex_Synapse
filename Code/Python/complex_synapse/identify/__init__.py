"""Fitting synapse models to experiment
"""
from . import baum_welch, fit_synapse, synapse_id
from .fit_synapse import print_callback
from .synapse_id import SynapseIdModel
from .plast_seq import PlasticitySequence, SimPlasticitySequence
from .fit_synapse import SynapseFitter, GroundedFitter, SynapseFitOptions
from .baum_welch import BaumWelchFitter, GroundedBWFitter, BaumWelchOptions
from .record import RecordingCallback, FitterReplay, GroundedFitterReplay
from .video import FitterVideo, VideoOptions, VideoLabels, VideoLayout
# from . import video, plast_seq

assert any((True, baum_welch, fit_synapse, synapse_id, print_callback,
            SynapseIdModel, PlasticitySequence, SimPlasticitySequence,
            SynapseFitOptions, SynapseFitter, GroundedFitter,
            BaumWelchOptions, BaumWelchFitter, GroundedBWFitter,
            FitterVideo, VideoOptions, VideoLabels, VideoLayout,
            RecordingCallback, FitterReplay, GroundedFitterReplay))
