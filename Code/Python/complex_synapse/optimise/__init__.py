# -*- coding: utf-8 -*-
"""
Code for studying complex synapses
"""
from . import optimise, sticky, shorten, plot, video
from .synapse_opt import SynapseOptModel, SynapseParamModel, TopologyOptions
from .optimise import ModelOptions, ProblemOptions, OptimOptions, OptimProblem
from .optimise import (normal_problem, shifted_problem, optim_laplace,
                       optim_laplace_range, reoptim_laplace_range,
                       check_cond_range, proven_envelope_laplace,
                       heuristic_envelope_laplace)
from .video import (VideoLabels, VideoLayout, VideoOptions,
                    GraphPlots, ModelPlots, EnvelopeFig, animate)
from ..options import ImageOptions

assert any((True, optimise, sticky, shorten, plot, video))
assert any((True, SynapseOptModel, SynapseParamModel, TopologyOptions))
assert any((True, ModelOptions, ProblemOptions, OptimOptions, OptimProblem))
assert any((True, normal_problem, shifted_problem, optim_laplace,
            optim_laplace_range, reoptim_laplace_range, check_cond_range,
            proven_envelope_laplace, heuristic_envelope_laplace))
assert any((True, ImageOptions, VideoLabels, VideoLayout, VideoOptions,
            GraphPlots, ModelPlots, EnvelopeFig, animate))
