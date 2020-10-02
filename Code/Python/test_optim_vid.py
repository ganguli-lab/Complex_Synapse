# -*- coding: utf-8 -*-
"""Test script for optimisation
"""
# %%
import numpy as np
import matplotlib.pyplot as plt
import sl_py_tools.matplotlib_tricks as mpt
import sl_py_tools.display_tricks as dt
import complex_synapse.optimise as cso
np.set_printoptions(precision=2, suppress=True, threshold=90)
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# %%
nst, sss, opts, envs, mods = cso.plot.load_data()
obj = cso.video.EnvelopeFig(sss, envs[1], mods[1], vmax=0.5)
ani = cso.video.animate(obj)
plt.show()
# %%
ani.save('~/Documents/videos/envelope.mp4',
         progress_callback=dt.FormattedTempDisplay("frame {:3d}/{:3d}"))
