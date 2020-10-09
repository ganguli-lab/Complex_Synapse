# -*- coding: utf-8 -*-
"""Test script for synapse identification
"""
from pathlib import Path
# import logging

# import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import numpy_linalg as la
import sl_py_tools.display_tricks as dt
import sl_py_tools.matplotlib_tricks as mpt

import complex_synapse as cs
import complex_synapse.identify as idfy

# -----------------------------------------------------------------------------
np.set_printoptions(precision=3, suppress=True, threshold=90)
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
fmt = 'mp4'
nfr = 25
# logging.basicConfig(filename='save_fitter.log', filemode='w', level=logging.INFO)
# logging.captureWarnings(True)
# -----------------------------------------------------------------------------
opt = idfy.BaumWelchOptions(disp_step=10, verbosity=2 + 6 + 9)
true_model = idfy.SynapseIdModel.build(cs.builders.build_serial, 6, jmp=0.7)
# -----------------------------------------------------------------------------
opt.verbosity = 9
with la.load('test_fit.npz') as saved_file:
    saved = dict(saved_file)
vid = idfy.FitterPlots(np.s_[:100, 0], transpose=False)
old_fit = idfy.GroundedFitterReplay(saved, callback=vid, opt=opt, max_it=nfr)
ani = idfy.animate(old_fit, blit=False)
disp = dt.FormattedTempDisplay("frame: {:2d}/{:2d}")
folder = Path('~/Documents/videos/identification').expanduser()
# -----------------------------------------------------------------------------
fname = 'identification.' + fmt
if fmt == 'mp4':
    writer = 'ffmpeg'
elif fmt == 'pdf':
    writer = 'pdf_pages'
elif fmt == 'png':
    writer = mpt.FileSeqWriter(fps=2, ndigit=2)
else:
    writer, fname = None, None
# -----------------------------------------------------------------------------
if writer and fname:
    ani.save(folder / fname, writer=writer, progress_callback=disp)
else:
    plt.show()
