import numpy as np

def calculate_UHF_S2(uhf_parameters, S, Na, Nb, log=False):
    C_a = uhf_parameters.expansion().alpha.matrix()
    C_b = uhf_parameters.expansion().beta.matrix()

    Sigma = C_a.T.conjugate() @ S @ C_b

    Sz = (Na - Nb) // 2

    S2 = Sz * (Sz + 1) + Nb

    spin_contamination = 0

    for i in range(Na):
        for j in range(Nb):
            spin_contamination -= np.abs(Sigma[i, j]) ** 2.0

    S2 += spin_contamination
    return S2
