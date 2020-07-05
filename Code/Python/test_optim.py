# -*- coding: utf-8 -*-
import numpy_linalg as la
import complex_synapse as cs

s = la.geomspace(1e-4, 10, 50)
m = cs.SynapseOptModel.rand(nst=10, npl=2)
# print(m.peq_min_fun(0.1))
envelope, models = cs.soptimise.optim_laplace_range(s, 10, repeats=10,
                                                    serial=False)
