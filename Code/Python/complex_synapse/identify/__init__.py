# -*- coding: utf-8 -*-
"""Fitting synapse models to experiment
"""
from . import baum_welch, synapse_id
from .baum_welch import BaumWelchFitter, BaumWelchOptions, GroundedBWFitter
from .fit_synapse import (GroundedFitter, SynapseFitOptions, SynapseFitter,
                          print_callback)
from .plast_seq import PlasticitySequence, SimPlasticitySequence
from .record import FitterReplay, GroundedFitterReplay, RecordingCallback
from .synapse_id import SynapseIdModel
from .video import FitterPlots, VideoLabels, VideoLayout, VideoOptions, animate

# from . import fit_synapse, video, plast_seq

assert any((True, baum_welch, synapse_id, print_callback, animate,
            SynapseIdModel, PlasticitySequence, SimPlasticitySequence,
            SynapseFitter, GroundedFitter, BaumWelchFitter, GroundedBWFitter,
            FitterReplay, GroundedFitterReplay, RecordingCallback, FitterPlots,
            SynapseFitOptions, BaumWelchOptions, VideoOptions,
            VideoLabels, VideoLayout,
            ))
