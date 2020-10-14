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
from sl_py_tools.time_tricks import time_with
np.set_printoptions(precision=3, suppress=True, threshold=90)
np.seterr(all='warn')
mpt.rc_colours()
mpt.rc_fonts('sans-serif')
# logging.basicConfig(filename='save_fitter.log', filemode='w', level=logging.INFO)
# logging.captureWarnings(True)

def compare(x: la.lnarray, y: la.lnarray):
    return np.abs(x - y).max()[()]

# %%
opt = idfy.BaumWelchOptions(disp_step=10, verbosity=2 + 6 + 9)
true_model = idfy.SynapseIdModel.build(cs.builders.build_serial, 6, jmp=0.7)
# %%
sim = true_model.simulate(20, 400)
fit_model = idfy.SynapseIdModel.rand(6, binary=True)
fit_model.normalise()
# %%
with time_with(subsec=True, absolute=False):
    update, initial = fit_model.updaters()
    sim_data = sim.move_t_axis(-1).plast_type, sim.move_t_axis(-1).readouts
    a, b, e = idfy.baum_welch._calc_bw_abe_c(update, initial, *sim_data)
    p, init = idfy.baum_welch._calc_model_c(update, *sim_data, (a, b, e))
with time_with(subsec=True, absolute=False):
    update, initial, plast_type = idfy.baum_welch._get_updaters(fit_model, sim)
    aa, bb, ee = idfy.baum_welch._calc_bw_abe(update, initial)
    pp, iinit = idfy.baum_welch._calc_model(update, plast_type, (aa, bb, ee))
print(compare(a, aa))
print(compare(p, pp))
# %%
if __name__ != "__main__":
    # %%
    opt.verbosity = 8
    rec = idfy.RecordingCallback(idfy.print_callback)
    fitter = idfy.GroundedBWFitter(sim, fit_model, true_model, rec, opt=opt)
    fitter.run()
    # %%
    np.savez_compressed('test_fit', **rec.info)
    # %%
    opt.verbosity = 9
    with la.load('test_fit.npz') as saved_file:
        saved = {**saved_file}
    # %%
    vid = idfy.FitterPlots(np.s_[0, :100], transpose=False)
    old_fit = idfy.GroundedFitterReplay(saved, callback=vid, opt=opt)
    # opt.max_it = 10
    # %%
    opt.disp_step = 5
    opt.verbosity = 8
    rec = idfy.RecordingCallback(idfy.print_callback)
    fitter = idfy.GroundedBWFitter.rerun(saved, callback=rec, opt=opt)
    fitter.run()
    # %%
    vid(old_fit, 0)
    vid(old_fit, 1)
    # folder = Path('~/Documents/videos').expanduser()
    # vid.fig.savefig(str(folder / 'identification_frame.pdf'), format='pdf')
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
