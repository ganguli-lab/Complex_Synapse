# -*- coding: utf-8 -*-
"""Optimising complex synapses
"""
from . import optimise, plot, shorten, sticky, video
from .optimise import (ModelOptions, OptimOptions, OptimProblem,
                       ProblemOptions, check_cond_range,
                       heuristic_envelope_laplace, normal_problem,
                       optim_laplace, optim_laplace_range,
                       proven_envelope_laplace, reoptim_laplace_range,
                       shifted_problem, equlibrium_envelope_laplace)
from .plot import GraphPlots, load_data, save_data
from .synapse_opt import SynapseOptModel
from .video import (EnvelopeFig, ModelPlots, VideoLabels, VideoLayout,
                    VideoOptions, animate)

assert any((True, optimise, sticky, shorten, plot, video))
assert any((True, SynapseOptModel, OptimProblem))
assert any((True, GraphPlots, ModelPlots, EnvelopeFig))
assert any((True, ModelOptions, ProblemOptions, OptimOptions, VideoLabels,
            VideoLayout, VideoOptions))
assert any((True, load_data, save_data, animate, normal_problem,
            shifted_problem, optim_laplace, optim_laplace_range,
            reoptim_laplace_range, check_cond_range, proven_envelope_laplace,
            heuristic_envelope_laplace, equlibrium_envelope_laplace))
