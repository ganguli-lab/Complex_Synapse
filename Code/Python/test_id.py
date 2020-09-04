"""Test script for synapse identification
"""
# %%
import numpy as np
import matplotlib as mpl
import matplotlib.animation as mpla
import matplotlib.pyplot as plt
import numpy_linalg as la
import complex_synapse as cs
import complex_synapse.identify as idfy
import sl_py_tools.matplotlib_tricks as mpt
import sl_py_tools.iter_tricks as it
np.set_printoptions(precision=3, suppress=True, threshold=90)
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# %%
opt = idfy.BaumWelchOptions(disp_step = 10, verbosity = 2 + 6 + 9)
true_model = idfy.SynapseIdModel.build(cs.builders.build_serial, 6, jmp=0.7)
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
opt.verbosity = 17
with np.load('test_fit.npz') as saved_file:
    saved = {**saved_file}
vid = idfy.FitterPlots(np.s_[:100, 0], transpose=False)
old_fit = idfy.GroundedFitterReplay(saved, callback=vid, opt=opt)
old_fit.init()
ani = mpla.FuncAnimation(vid.fig, old_fit.step, init_func=old_fit.init,
                         frames=it.erange(old_fit.opt.max_it),
                         interval=500, repeat_delay=1000, blit=False)
plt.show()
# %%
fitter = idfy.GroundedBWFitter(sim, fit_model, true_model, opt=opt)
# %%
vid = idfy.FitterPlots(np.s_[:100, 0], transpose=False)
vid(fitter, 0)
# %%
fitter.run()
# %%
vid = idfy.FitterPlots(fitter, np.s_[:100, 0],
                       'C:/Users/subhy/Documents/videos/Fit/fit_{:03d}.png')
fitter.run(vid)
# %%
vo = idfy.VideoOptions()
vo
# %%
