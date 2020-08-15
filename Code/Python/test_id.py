# %%
import numpy as np
import numpy_linalg as la
import matplotlib as mpl
import matplotlib.pyplot as plt
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
vid = idfy.FitterVideo(fitter, np.s_[:, 0], 'fit_{}.png')
# %%
fitter.opt['verbose'] = 1
fitter.run(vid)
# %%
