"""Creating figures for Examples
"""
import typing as ty

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import numpy_linalg as la
import sl_py_tools.arg_tricks as ag
import sl_py_tools.matplotlib_tricks as mplt

import complex_synapse as cs

np.set_printoptions(precision=3, suppress=False, linewidth=90)
mplt.rc_fonts()
mplt.rc_colours()

def ser_casc_data(nst: int, jmp: float, time: la.lnarray
                  ) -> ty.Tuple[la.lnarray, la.lnarray, la.lnarray]:
    """Data for example plot

    Parameters
    ----------
    nst : int
        Number of states
    jmp : float
        Model parameter
    time : la.lnarray
        Vector of times to evaluate SNR

    Returns
    -------
    time : la.lnarray
        Vector of times at which SNR is evaluated
    serial_snr : la.lnarray
        SNR for serial model at each time
    cascade_snr : la.lnarray
        SNR for cascade model at each time
    """
    nst = ag.default(nst, 10)
    time = ag.default(time, 1 / la.geomspace(1e-4, 10, 50))
    jmp = ag.default(jmp, 0.7)

    serial = cs.SynapseMemoryModel.build(cs.builders.build_serial,
                                         nst, jmp=jmp)
    cascade = cs.SynapseMemoryModel.build(cs.builders.build_cascade,
                                          nst, jmp=jmp)


    serial_snr = serial.snr_exp_ave(time)
    cascade_snr = cascade.snr_exp_ave(time)

    return time, serial_snr, cascade_snr


def ser_casc_plot(time: la.lnarray, serial_snr: la.lnarray,
                  cascade_snr: la.lnarray) -> mpl.figure.Figure:
    """Example plot

    Parameters
    ----------
    time : la.lnarray
        Vector of times at which SNR is evaluated
    serial_snr : la.lnarray
        SNR for serial model at each time
    cascade_snr : la.lnarray
        SNR for cascade model at each time
    """
    fig, axs = plt.subplots()
    axs.loglog(time, serial_snr, label='serial')
    axs.loglog(time, cascade_snr, label='cascade')
    axs.set_xlabel(r"Time, $r\tau$")
    axs.set_ylabel("SNR")
    axs.set_xlim(1e-1, 1e4)
    axs.legend(loc="lower left")
    mplt.clean_axes(axs)
    return fig


if __name__ == "__main__":
    j, n, t = 0.7, 10, 1 / la.geomspace(1e-4, 10, 50)
    t, snr_serial, snr_cascade = ser_casc_data(n, j, t)
    fig_sc = ser_casc_plot(t, snr_serial, snr_cascade)
    fig_sc.savefig("../../Notes/Figures/serial_vs_cascade.pdf")
