# -*- coding: utf-8 -*-
"""Test script for optimisation
"""
# %%
import numpy as np
import matplotlib as mpl
import matplotlib.animation as mpla
import matplotlib.pyplot as plt
import numpy_linalg as la
import complex_synapse as cs
import complex_synapse.optimise as cso
import complex_synapse.graphs as csg
import complex_synapse.optimise.plot as csp
import sl_py_tools.matplotlib_tricks as mpt
from sl_py_tools.time_tricks import time_with
np.set_printoptions(precision=2, suppress=True)
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# %%
m = cs.SynapseMemoryModel.rand(nst=10, npl=2)
# print(m.peq_min_fun(0.1))
# %%
s = la.geomspace(1e-4, 10, 5)
with time_with():
    envelope, models = cso.optim_laplace_range(
        s, 10, repeats=10, cond_lim=True, cond_thresh=1e3, serial=False,
        maker=cso.normal_problem, method='SLSQP')
# %%
nst, sss, opts, envs, mods = csp.load_data()
# %%
obj = cso.video.EnvelopeFig(sss, envs[3], envs[1], mods[1])
obj.update(40)
# %%
serial = cs.optimise.SynapseParamModel.rand(
    6, serial=True, directions=(1, -1), binary=True)
graph = csg.param_model_to_graph(serial)
nodes, edges = csg.draw_graph(graph)
# %%
