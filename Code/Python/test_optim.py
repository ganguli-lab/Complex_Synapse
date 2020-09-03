# -*- coding: utf-8 -*-
"""Test script for optimisation
"""
# %%
import numpy as np
import numpy_linalg as la
import complex_synapse as cs
import complex_synapse.optimise as cso
from sl_py_tools.time_tricks import time_with
np.set_printoptions(precision=2, suppress=True)
# %%
m = cs.SynapseMemoryModel.rand(nst=10, npl=2)
# print(m.peq_min_fun(0.1))
# %%
s = la.geomspace(1e-4, 10, 5)
with time_with():
    envelope, models = cso.optim_laplace_range(
        s, 10, repeats=10, cond_lim=True, cond_thresh=1e3,
        maker=cso.normal_problem, method='SLSQP')
# %%
nst, sss, opts, envs, mods = cso.plot.load_data()
fig, obj = cso.video.make_plots(sss, envs[1::2], mods[1][0], 0)
# %%
fig, obj = cso.video.make_plots(sss, envs[1::2], mods[1][45], 45)
# %%
fig, obj = cso.video.make_plots(sss, envs[1::2], mods[1][45], 45)
obj.update_plots(sss, mods[1][25], 25)

# %%
