# -*- coding: utf-8 -*-
"""Test script for synapse identification
"""
from pathlib import Path
# import logging

# import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# import numpy_linalg as la
import sl_py_tools.display_tricks as dt
import sl_py_tools.matplotlib_tricks as mpt

import complex_synapse as cs
import complex_synapse.identify as idfy

# -----------------------------------------------------------------------------
np.set_printoptions(precision=3, suppress=True, threshold=90)
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
fmt = 'pdf'
# logging.basicConfig(filename='save_fitter.log', filemode='w', level=logging.INFO)
# logging.captureWarnings(True)
# -----------------------------------------------------------------------------
opt = idfy.BaumWelchOptions(disp_step=10, verbosity=2 + 6 + 9)
true_model = idfy.SynapseIdModel.build(cs.builders.build_serial, 6, jmp=0.7)
# -----------------------------------------------------------------------------
opt.verbosity = 9
with np.load('test_fit.npz') as saved_file:
    saved = {**saved_file}
vid = idfy.FitterPlots(np.s_[:100, 0], transpose=False)
old_fit = idfy.GroundedFitterReplay(saved, callback=vid, opt=opt)
# opt.max_it = 10
ani = idfy.animate(old_fit, blit=False)
disp = dt.FormattedTempDisplay("frame: {:3d}/{:3d}")
folder = Path('~/Documents/videos').expanduser()
# -----------------------------------------------------------------------------
if fmt == 'mp4':
    writer, fname = 'ffmpeg', 'identification.mp4'
elif fmt == 'pdf':
    writer, fname = 'pdf_pages', 'identification.pdf'
elif fmt == 'png':
    writer = mpt.FileSeqWriter(fps=2, ndigit=3)
    fname = 'identification/test_.png'
else:
    writer, fname = None, None
# -----------------------------------------------------------------------------
if writer and fname:
    ani.save(folder / fname, writer=writer, progress_callback=disp)
else:
    plt.show()
