"""Fitting synapse models to experiment
"""
from . import baum_welch, synapse_id
from .fit_synapse import print_callback
from .synapse_id import SynapseIdModel
from .plast_seq import PlasticitySequence, SimPlasticitySequence
from .fit_synapse import SynapseFitter, GroundedFitter, SynapseFitOptions
from .baum_welch import BaumWelchFitter, GroundedBWFitter, BaumWelchOptions
from .record import RecordingCallback, FitterReplay, GroundedFitterReplay
from .video import FitterPlots, VideoOptions, VideoLabels, VideoLayout, animate
from ..options import ImageOptions, AnimationOptions
# from . import fit_synapse, video, plast_seq

assert any((True, baum_welch, synapse_id, print_callback, animate,
            SynapseIdModel, PlasticitySequence, SimPlasticitySequence,
            SynapseFitter, GroundedFitter, BaumWelchFitter, GroundedBWFitter,
            FitterReplay, GroundedFitterReplay, RecordingCallback, FitterPlots,
            SynapseFitOptions, BaumWelchOptions, VideoOptions,
            VideoLabels, VideoLayout, ImageOptions, AnimationOptions,
            ))
