"""Test script for synapse identification
"""
# %%
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy_linalg as la
import complex_synapse as cs
import complex_synapse.identify as idfy
import sl_py_tools.matplotlib_tricks as mpt
np.set_printoptions(precision=3, suppress=True, threshold=90)
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# %%
true_model = idfy.SynapseIdModel.build(cs.builders.build_serial, 6, jmp=0.7)
sim = true_model.simulate(100, 10)
fit_model = idfy.SynapseIdModel.rand(6, binary=True)
fit_model.normalise()
fitter = idfy.baum_welch.GroundedBWFitter(sim, fit_model, true_model)
# %%
fitter.opt.verbosity = 2 + 6 + 9
fitter.opt.disp_step = 10
# %%
rec = idfy.RecordingCallback(np.s_[:100, 0], idfy.print_callback)
fitter.run(rec)
# %%
old_fit = idfy.GroundedFitterReplay(sim, rec.info, verbosity=17)
vid = idfy.FitterVideo(np.s_[:100, 0], transpose=False)
old_fit.run(vid)
# %%
vid = idfy.FitterVideo(np.s_[:100, 0], transpose=False)
vid(fitter, 0)
# %%
fitter.run()
# %%
vid = idfy.FitterVideo(fitter, np.s_[:100, 0],
                       'C:/Users/subhy/Documents/videos/Fit/fit_{:03d}.png')
fitter.run(vid)
# %%
vo = idfy.VideoOptions()
vo
# %%
