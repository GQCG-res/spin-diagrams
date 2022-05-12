import gqcpy

def create_H2_PES(R_values):
    H_left = gqcpy.Nucleus(1, 0.0,0.0,0.0)

    molecules = []
    for R in R_values:
        H_right = gqcpy.Nucleus(1, R,0.0,0.0)  # The molecule is aligned on the x-axis.
        molecules.append(gqcpy.Molecule([H_left, H_right]))

    return molecules


def create_ring_PES(n, R_values):
    molecules = []

    for R in R_values:
        molecules.append(gqcpy.Molecule.HRingFromDistance(n, R))

    return molecules
