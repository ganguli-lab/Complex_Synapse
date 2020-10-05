# -*- coding: utf-8 -*-
"""Test script for synapse identification
"""
# %%
from pathlib import Path
# import logging
import numpy as np
import matplotlib as mpl
import matplotlib.animation as mpla
import matplotlib.pyplot as plt
import numpy_linalg as la
import complex_synapse as cs
import complex_synapse.identify as idfy
import sl_py_tools.matplotlib_tricks as mpt
# import sl_py_tools.iter_tricks as it
import sl_py_tools.display_tricks as dt
np.set_printoptions(precision=3, suppress=True, threshold=90)
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# logging.basicConfig(filename='save_fitter.log', filemode='w', level=logging.INFO)
# logging.captureWarnings(True)
# %%
opt = idfy.BaumWelchOptions(disp_step=10, verbosity=2 + 6 + 9)
true_model = idfy.SynapseIdModel.build(cs.builders.build_serial, 6, jmp=0.7)
# %%
opt.verbosity = 9
with np.load('test_fit.npz') as saved_file:
    saved = {**saved_file}
vid = idfy.FitterPlots(np.s_[:100, 0], transpose=False)
old_fit = idfy.GroundedFitterReplay(saved, callback=vid, opt=opt)
opt.max_it = 10
folder = Path('~/Documents/videos').expanduser()
# %%
# vid(old_fit, 0)
# vid.fig.savefig(str(folder / 'identification_frame.pdf'), format='pdf', dpi=72)
# %%
# ani = idfy.animate(old_fit, blit=False)
# plt.show()
# %%
# writer = mpt.FileSeqWriter(fps=2)
# fname = str(folder / 'identification/test_.pdf')
# writer.setup(vid.fig, fname, ndigit=3)
# %%
# writer = 'ffmpeg'
# fname = str(folder / 'identification.mp4')
# %%
writer = 'pdf_pages'
fname = str(folder / 'identification.pdf')
# %%
ani = idfy.animate(old_fit, blit=False)
ani.save(fname, writer=writer, progress_callback=dt.FormattedTempDisplay(
    "frame: {:3d}/{:3d}"))
# %%
# writer = mpt.PdfPagesWriter(fps=2)
# old_fit.init()
# with writer.saving(vid.fig, fname, dpi=None):
#     for i in it.dcount('frame', old_fit.opt.max_it):
#         old_fit.step(i)
#         writer.grab_frame()
# %%
if __name__ != "__main__":
    # %%
    sim = true_model.simulate(100, 10)
    fit_model = idfy.SynapseIdModel.rand(6, binary=True)
    fit_model.normalise()
    # %%
    opt.verbosity = 8
    rec = idfy.RecordingCallback(np.s_[:100, 0], idfy.print_callback)
    fitter = idfy.GroundedBWFitter(sim, fit_model, true_model, rec, opt=opt)
    fitter.run()
    np.savez_compressed('test_fit', **rec.info)
    # %%
    fitter = idfy.GroundedBWFitter(sim, fit_model, true_model, opt=opt)
    vid = idfy.FitterPlots(np.s_[:100, 0], transpose=False)
    vid(fitter, 0)
    # %%
    fitter = idfy.GroundedBWFitter(sim, fit_model, true_model, opt=opt)
    fitter.run()
    # %%
    fitter = idfy.GroundedBWFitter(sim, fit_model, true_model, opt=opt)
    vid = idfy.FitterPlots(np.s_[:100, 0])
    fitter.run(vid)
    # %%
    vo = idfy.VideoOptions()
    vo
    # %%
