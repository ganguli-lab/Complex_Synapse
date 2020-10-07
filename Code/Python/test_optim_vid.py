# -*- coding: utf-8 -*-
"""Test script for optimisation
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sl_py_tools.matplotlib_tricks as mpt
import sl_py_tools.display_tricks as dt
import complex_synapse.optimise as cso
# -----------------------------------------------------------------------------
np.set_printoptions(precision=2, suppress=True, threshold=90)
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
fmt = 'pdf'
# -----------------------------------------------------------------------------
nst, sss, opts, envs, mods = cso.plot.load_data()
vid = cso.video.EnvelopeFig(sss, envs[1], mods[1], vmax=0.5)
ani = cso.video.animate(vid)
disp = dt.FormattedTempDisplay("frame: {:2d}/{:2d}")
folder = Path('~/Documents/videos').expanduser()
# -----------------------------------------------------------------------------
if fmt == 'mp4':
    writer, fname = 'ffmpeg', 'envelope.mp4'
elif fmt == 'pdf':
    writer, fname = 'pdf_pages', 'envelope.pdf'
elif fmt == 'png':
    writer = mpt.FileSeqWriter(fps=2, ndigit=3)
    fname = 'envelope/test_.png'
else:
    writer, fname = None, None
# -----------------------------------------------------------------------------
if writer and fname:
    ani.save(folder / fname, writer=writer, progress_callback=disp)
else:
    plt.show()
