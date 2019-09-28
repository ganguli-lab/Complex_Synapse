# -*- coding: utf-8 -*-
import numpy_linalg as la
import complex_synapse as cs

s = la.geomspace(1e-4, 10, 50)
envelope, models = cs.optimise.optim_laplace_range(s, 10, repeats=10,
                                                   serial=True)
