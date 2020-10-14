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
beta : la.lnarray, (T,M)
    Scaled Baum-Welch backward variable.
eta : la.lnarray, (T,)
    Norm of Baum-Welch forward variable.
"""

plast_init: ufunc
plast_init__doc__ = """One Baum-Welch/Rabiner-Juang update of the model

Parameters
-----------
updaters : la.lnarray, (R,P,M,M) float[0:1]
    Plasticity matrices multiplied by readout probability given 'to' state.
plast_type : la.lnarray, (T-1,E), int[0:P]
    id of plasticity type after each time-step
alpha : la.lnarray, (T,M) float[0:1]
    Normalised Baum-Welch forward variable
beta : la.lnarray, (T,M) float
    Scaled Baum-Welch backward variable
eta : la.lnarray, (T,) float[1:]
    Norm of Baum-Welch forward variable

Returns
-------
plast : array_like, (P,M,M), float[0:1]
    new estimate of transition probability matrix.
initial : array_like, (M,) float[0:1]
    new estimate of distribution of initial state,
    not assumed to be the steady-state distribution.
"""

plast_steady: ufunc
plast_steady__doc__ = """One Baum-Welch/Rabiner-Juang update of the model

Parameters
-----------
updaters : la.lnarray, (R,P,M,M) float[0:1]
    Plasticity matrices multiplied by readout probability given 'to' state.
plast_type : la.lnarray, (T-1,E), int[0:P]
    id of plasticity type after each time-step
alpha : la.lnarray, (T,M) float[0:1]
    Normalised Baum-Welch forward variable
beta : la.lnarray, (T,M) float
    Scaled Baum-Welch backward variable
eta : la.lnarray, (T,) float[1:]
    Norm of Baum-Welch forward variable

Returns
-------
plast : array_like, (P,M,M), float[0:1]
    new estimate of transition probability matrix.
initial : array_like, (M,) float[0:1]
    new estimate of distribution of initial state,
    assumed to be the steady-state distribution.
"""
