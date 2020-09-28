# -*- coding: utf-8 -*-
"""Test script for optimisation
"""
# %%
from complex_synapse.optimise.video import serial_eqp
import numpy as np
import matplotlib as mpl
import matplotlib.animation as mpla
import matplotlib.pyplot as plt
import numpy_linalg as la
import complex_synapse as cs
import complex_synapse.optimise as cso
import complex_synapse.graphs as csg
import sl_py_tools.matplotlib_tricks as mpt
from sl_py_tools.time_tricks import time_with
np.set_printoptions(precision=2, suppress=True)
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# %%
nst, sss, opts, envs, mods = cso.plot.load_data()
obj = cso.video.EnvelopeFig(sss, envs[3], envs[1], mods[1], vmax=0.5)
ani = cso.video.animate(obj)
plt.show()
# %%
