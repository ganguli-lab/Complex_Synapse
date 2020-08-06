# %%
import numpy as np
import numpy_linalg as la
import complex_synapse as cs
import complex_synapse.identify as idfy
np.set_printoptions(precision=3, suppress=True, threshold=90)
# %%
true_model = idfy.SynapseIdModel.build(cs.builders.build_serial, 6, jmp=0.7)
sim = true_model.simulate(100, 10)
true_nllike = idfy.baum_welch.likelihood(true_model, sim)
fit_model = idfy.SynapseIdModel.rand(6, binary=True)
fit_model.normalise()
nllike = idfy.baum_welch.likelihood(fit_model, sim)
distance = (true_model - fit_model).norm()
print(f"-log P(readouts|true) = {true_nllike:.3f}")
print(f"-log P(readouts|fit) = {nllike:.3f} ||true-fit|| = {distance:.3f}")
# %%
for _ in range(100):
    prob, nllike = idfy.baum_welch.update_model(fit_model, sim)
    distance = (true_model - fit_model).norm()
    print(f"-log P(readouts|fit) = {nllike:.3f} ||true-fit|| = {distance:.3f}")
# %%
