"""Baum welch forward backward variables"""
from numpy import ufunc

alpha_beta: ufunc
alpha_beta__doc__ = """Calculate BW forward/backward variables

Parameters
-----------
updaters : la.lnarray, (R,P,M,M)
    Plasticity matrices multiplied by readout indicators of 'to' state,
    `Prob(i(t+1)=j, w(t+1)=r|i(t)=i, mu(t)=p)`.
initial : la.lnarray, (R,M)
    Initial state distribution multiplied by readout indicators of state,
    `Prob(w(0)=r, i(0)=i)`.
plast_type : ArrayLike, (T-1,E), int[0:P]
    id of plasticity type after each time-step.
readouts : ArrayLike, (T,E), int[0:R]
    id of readout from synapse at each time-step.

Returns
-------
alpha : la.lnarray, (T,M)
    Normalised Baum-Welch forward variable.
beta : la.lnarray, (T,M)\n"
    Scaled Baum-Welch backward variable.
eta : la.lnarray, (T,)\n"
    Norm of Baum-Welch forward variable.
"""
