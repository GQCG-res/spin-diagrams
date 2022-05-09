# -*- coding: utf-8 -*-
"""File containing short-hand functions for quantum chemical methods."""
import gqcpy
import numpy as np


def do_magnetic_UHF_calculation(molecule, basis_set, B, Sz, C_initial=None, threshold=1.0e-08, max_iter=10000, log=False, stability=True):
    """
        @param Sz:              The spin projection to target.
        @param C_initial:       An initial guess for the UHF expansion coefficients.

        @return:
            - The UHF energy.
            - The converged UHF parameters.
    """

    # Set up the molecular magnetic Hamiltonian in the AO basis.
    spin_orbital_basis = gqcpy.LondonUSpinOrbitalBasis(molecule, basis_set, B)
    hamiltonian_AO = spin_orbital_basis.quantize(gqcpy.FQMolecularPauliHamiltonian(molecule, B))

    S = spin_orbital_basis.quantize(gqcpy.OverlapOperator())


    # Do an UHF calculation, calculate the constant spin Zeeman contribution and calculate the total energy.
    N = molecule.numberOfElectrons()

    N_alpha = int(0.5 * N + Sz)
    N_beta = int(0.5 * N - Sz)
    print("N_alpha, N_beta: {}, {}".format(N_alpha, N_beta))

    K_alpha = spin_orbital_basis.numberOfSpinors() // 2
    K_beta = K_alpha
    if log: print("K_alpha, K_beta: {}, {}".format(K_alpha, K_beta))

    Oa, Va = N_alpha, K_alpha - N_alpha
    Ob, Vb = N_beta, K_beta - N_beta
    if log: print("(Oa, Va): ({},{})".format(Oa, Va))
    if log: print("(Ob, Vb): ({},{})".format(Ob, Vb))

    if C_initial is None:
        environment = gqcpy.UHFSCFEnvironment_cd.WithCoreGuess(N_alpha, N_beta, hamiltonian_AO, S)
    else:
        environment = gqcpy.UHFSCFEnvironment_cd(N_alpha, N_beta, hamiltonian_AO, S, C_initial)

    solver = gqcpy.UHFSCFSolver_cd.Plain(threshold=threshold, maximum_number_of_iterations=max_iter)

    qc_structure = gqcpy.UHF_cd.optimize(solver, environment)
    UHF_parameters = qc_structure.groundStateParameters()

    C = UHF_parameters.expansion()
    hamiltonian_MO = hamiltonian_AO.transformed(C)
    stability_matrices = UHF_parameters.calculateStabilityMatrices(hamiltonian_MO)
    internal_stability = stability_matrices.isInternallyStable()
    if stability:
        while internal_stability is False:
            if log: print("Following internal instability.")
            U = stability_matrices.instabilityRotationMatrix(Oa, Ob, Va, Vb)
        
            C.rotate(U)

            environment2 = gqcpy.UHFSCFEnvironment_cd(N_alpha, N_beta, hamiltonian_AO, S, C)

            qc_structure = gqcpy.UHF_cd.optimize(solver, environment2)
            UHF_parameters = qc_structure.groundStateParameters()
            C2 = UHF_parameters.expansion()

            # Check for stability of the GHF solution.
            hamiltonian_MO2 = hamiltonian_AO.transformed(C2)
            stability_matrices2 = UHF_parameters.calculateStabilityMatrices(hamiltonian_MO2)
            if log: stability_matrices2.printStabilityDescription()

            internal_stability = stability_matrices2.isInternallyStable()

    electronic_energy = qc_structure.groundStateEnergy()
    nuclear_repulsion = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value()

    UHF_energy = (electronic_energy + nuclear_repulsion).real  # Omit the imaginary component.

    return UHF_energy, UHF_parameters


def do_magnetic_GHF_calculation(molecule, basis_set, B, mu=None, log=False, threshold=1.0e-08, max_iter=10000, solver_algorithm='Plain', stability=True):
    """
        @param stability            A boolean indicating if stability analysis should be performed.
        @param mu                   A 3-vector of Lagrange multipliers, signaling to do a modified GHF calculation.
                                    The modification consists of (-mu . S), where is the spin vector operator.
        @param log                  A boolean indicating if logging output should be given.

        @return:
            - The GHF energy.
            - The converged GHF parameters.
    """

    N = molecule.numberOfElectrons()

    # Set up the molecular magnetic Pauli Hamiltonian in the AO basis.
    spinor_basis = gqcpy.LondonGSpinorBasis(molecule, basis_set, B)
    hamiltonian_AO = spinor_basis.quantize(gqcpy.FQMolecularPauliHamiltonian(molecule, B))

    if mu is not None:
        # Apply the modification potential to the Hamiltonian.
        Spin_op = spinor_basis.quantize(gqcpy.ElectronicSpinOperator())
        S_parameters = Spin_op.allParameters()

        W = gqcpy.ScalarGSQOneElectronOperator_cd(mu[0] * S_parameters[0] + mu[1] * S_parameters[1] + mu[2] * S_parameters[2])
        hamiltonian_AO -= W


    # Do a GHF calculation and calculate the total energy.
    S_AO = spinor_basis.quantize(gqcpy.OverlapOperator())

    # Initialize the GHF SCF environment with an initial guess that does not consist of spin-orbitals, but consists of true spinors that are not eigenvectors of S_z.
    environment = gqcpy.GHFSCFEnvironment_cd.WithComplexlyTransformedCoreGuess(N, hamiltonian_AO, S_AO)

    if solver_algorithm == 'Plain':
        solver = gqcpy.GHFSCFSolver_cd.Plain(threshold=threshold, maximum_number_of_iterations=max_iter)
    elif solver_algorithm == 'DIIS':
        solver = gqcpy.GHFSCFSolver_cd.DIIS(threshold=threshold, maximum_number_of_iterations=max_iter)
    else:
        raise ValueError("An unspecified solver was supplied.")

    qc_structure = gqcpy.GHF_cd.optimize(solver, environment)
    GHF_parameters = qc_structure.groundStateParameters()

    # Check for stability of the GHF solution.
    M = spinor_basis.numberOfSpinors()

    C = GHF_parameters.expansion()
    hamiltonian_MO = hamiltonian_AO.transformed(C)
    stability_matrices = GHF_parameters.calculateStabilityMatrices(hamiltonian_MO)
    internal_stability = stability_matrices.isInternallyStable()
    
    if stability:
        while internal_stability is False:
            if log: print("Following internal instability.")
            U = stability_matrices.instabilityRotationMatrix(N, M-N)

            C.rotate(U)
            environment2 = gqcpy.GHFSCFEnvironment_cd(N, hamiltonian_AO, S_AO, C)

            qc_structure = gqcpy.GHF_cd.optimize(solver, environment2)
            GHF_parameters = qc_structure.groundStateParameters()
            C2 = GHF_parameters.expansion()

            # Check for stability of the GHF solution.
            hamiltonian_MO2 = hamiltonian_AO.transformed(C2)
            stability_matrices2 = GHF_parameters.calculateStabilityMatrices(hamiltonian_MO2)
            if log: stability_matrices2.printStabilityDescription()

            internal_stability = stability_matrices2.isInternallyStable()

    # Calculate the GHF total energy.
    electronic_energy = qc_structure.groundStateEnergy()
    nuclear_repulsion = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value()

    GHF_energy = (electronic_energy + nuclear_repulsion).real  # Omit the imaginary component, since it's zero.

    if mu is not None:
        # Apply the modification correction to the energy.

        S_AO = GHF_parameters.calculateExpectationValueOf(gqcpy.ElectronicSpinOperator(), S_AO)
        E_modification = mu.dot(S_AO).real

        GHF_energy += E_modification  # + because we apply the modification -(mu.S)

    return GHF_energy, GHF_parameters



def do_magnetic_FCI_calculation(molecule, basis_set, B, mu=None, n=1, basis="GHF", GHF_threshold=1.0e-03):
    """
        @param basis:   "GHF" or "Lowdin", indicating a GHF or Löwdin orbital basis.
        @param n:       The number of states to be found.
        @param mu       A 3-vector of Lagrange multipliers, signaling to do a modified FCI calculation.
                        The modification consists of (-mu . S), where is the spin vector operator.

        @return:
            - The ground state FCI energy.
            - The converged FCI ground state parameters.
            - The orbitals used to define the GHF reference.
    """

    N = molecule.numberOfElectrons()

    # Set up the basis of London ASOs.
    spinor_basis = gqcpy.LondonGSpinorBasis(molecule, basis_set, B)
    if basis == "GHF":

        # Set up the molecular magnetic Hamiltonian in the AO basis.
        hamiltonian = spinor_basis.quantize(gqcpy.FQMolecularPauliHamiltonian(molecule, B))

        # Optimize the GHF wave function model.
        (GHF_energy, GHF_parameters) = do_magnetic_GHF_calculation(molecule, basis_set, B, threshold=GHF_threshold)
        C = GHF_parameters.expansion()

        spinor_basis.transform(C)
        hamiltonian.transform(C)

    elif basis == "Lowdin":
        spinor_basis.lowdinOrthonormalize()
        hamiltonian = spinor_basis.quantize(gqcpy.FQMolecularPauliHamiltonian(molecule, B))

        C = spinor_basis.expansion()

    else: raise ValueError("The orbital basis must be either 'GHF' or 'Löwdin'.")


    if mu is not None:
        # Apply the modification potential to the Hamiltonian.
        Spin_op = spinor_basis.quantize(gqcpy.ElectronicSpinOperator())
        S_parameters = Spin_op.allParameters()

        W = gqcpy.ScalarGSQOneElectronOperator_cd(mu[0] * S_parameters[0] + mu[1] * S_parameters[1] + mu[2] * S_parameters[2])
        hamiltonian -= W

    # Set up the FCI ONV basis and diagonalize the dense Hamiltonian matrix.
    M = spinor_basis.numberOfSpinors()
    N = molecule.numberOfElectrons()

    full_onv_basis = gqcpy.SpinUnresolvedONVBasis(M, N)
    onv_basis = gqcpy.SpinUnresolvedSelectedONVBasis(full_onv_basis)

    environment = gqcpy.CIEnvironment.Dense_cd(hamiltonian, onv_basis)
    solver = gqcpy.EigenproblemSolver.Dense_cd()
    qc_structure = gqcpy.CI_cd(onv_basis, n).optimize(solver, environment)

    nuclear_repulsion = gqcpy.NuclearRepulsionOperator(molecule.nuclearFramework()).value()
    if n == 1:  # We only want the ground state.
        FCI_parameters = qc_structure.groundStateParameters()
        FCI_energy = qc_structure.groundStateEnergy().real + nuclear_repulsion

        if mu is not None:
            # Apply the modification correction to the energy.

            D = FCI_parameters.calculate1DM()
            S_vector = np.array(Spin_op.calculateExpectationValue(D)).real
            E_modification = mu.dot(S_vector).real

            FCI_energy += E_modification  # + because we apply the modification -(mu.S)

        return FCI_energy, FCI_parameters, C

    else:

        if mu is not None: raise ValueError("Multiple states in combination with a multiplier is unsupported.")

        return [qc_structure.energy(i).real + nuclear_repulsion for i in range(n)], [qc_structure.parameters(i) for i in range(n)], C
