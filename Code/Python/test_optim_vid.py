# -*- coding: utf-8 -*-
"""Test script for optimisation
"""
# %%
from pathlib import Path
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
folder = Path('~/Documents/videos').expanduser()
# %%
plt.show()
# %%
# ani.save(str(folder / 'envelope.mp4'),
#          progress_callback=dt.FormattedTempDisplay("frame {:3d}/{:3d}"))
# %%
# writer = mpt.FileSeqWriter(fps=2)
# writer.setup(obj.fig, str(folder / 'envelope/test_.pdf'), None, 2)
# ani.save(folder / 'envelope/test_.pdf',
#          progress_callback=dt.FormattedTempDisplay("frame {:3d}/{:3d}"))
