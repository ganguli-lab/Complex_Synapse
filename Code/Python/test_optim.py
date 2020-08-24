# -*- coding: utf-8 -*-
"""Test script for optimisation
"""
import numpy_linalg as la
import complex_synapse as cs
from sl_py_tools.time_tricks import time_with

s = la.geomspace(1e-4, 10, 5)
# m = cs.SynapseOptModel.rand(nst=10, npl=2)
# print(m.peq_min_fun(0.1))
with time_with():
    envelope, models = cs.optimise.optim_laplace_range(
        s, 10, repeats=10, cond=True, CondThresh=1e3,
        maker=cs.optimise.normal_problem, method='trust-constr')
