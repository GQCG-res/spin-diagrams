"""This file contains functionality related to constraining spin expectation values for the GHF wave function model."""

import gqcpy

import numpy as np
import scipy.optimize

from prototypes.methods import *
from prototypes.optimization import *


def find_constrained_GHF_state(molecule, basis_set, B, S_z_target, threshold=1.0e-08, max_iter=10000, bracket=None, tolerance=1.0e-08, log=False, method='brentq', solver='Plain', stability=True):
    """
        Attempt to find a GHF state with the requested spin expectation value, using a modified (mu) algorithm.

        @param S_z_target                   The target S_z expectation value.
        @param bracket                      A bracket (tuple) in which the target multiplier (mu) is to be found.

        @param stability                    A boolean indicating if stability analysis should be performed.
        @param tolerance                    The tolerance for the requested spin expectation value.
        @param log                          A boolean indicating if logging output should be given.
        @param method                       The scipy.optimize.root_scalar method that should be used to find the root mu.

        @return:
            - The constrained GHF energy.
            - The constrained GHF state.
            - The 'optimal' Lagrange multiplier.
    """

    def f(mu):
        """
            The function whose root is to be found / should be minimized:
                f(mu) = <S_z>(mu) - S_target
        """

        #if log: print("\tEvaluating f(mu) with mu = {}".format(mu))

        # Do a GHF calculation with the modification dictated by the Lagrange multiplier.
        mu_vector = mu * np.array([0, 0, 1])
        (E_constrained, GHF_parameters) = do_magnetic_GHF_calculation(molecule, basis_set, B, threshold=threshold, max_iter=max_iter, mu=mu_vector, log=log, solver_algorithm=solver, stability=stability)

        # Calculate the S_z expectation value for the modified GHF parameters.
        spinor_basis = gqcpy.LondonGSpinorBasis(molecule, basis_set, B)
        S_AO = spinor_basis.overlap()

        S_calculated_vector = GHF_parameters.calculateExpectationValueOf(gqcpy.ElectronicSpinOperator(), S_AO)
        S_z_calculated = S_calculated_vector[2].real

        #if log: print("\tFound S_z expectation value of {}".format(S_z_calculated))

        return S_z_calculated - S_z_target

    # When trying to reach the maximum or minimum S_z value that the molecule can sustain, a root-finding algorithm won't work, so we use a line search algorithm instead.
    # In the minimum case, we expect that a Lagrange multiplier = -B will work.
    # In the maximum case, we expect that a very large Lagrange multiplier will work.
    S_z_max = molecule.numberOfElectrons() / 2
    S_z_min = (molecule.numberOfElectrons() % 2) / 2

    targeting_S_z_max = (S_z_target == S_z_max)
    targeting_S_z_min = (S_z_target == S_z_min)
    targeting_extremal_S_z = (targeting_S_z_max or targeting_S_z_min)

    if targeting_extremal_S_z:

        # Do a line search algorithm, but set up the bracket smartly.
        if targeting_S_z_min:
            if log: print("Finding minimal S_z value.")
            if bracket is None: bracket = (-0.1, 0.1)

            result = scipy.optimize.minimize_scalar(lambda mu: f(mu) ** 2, method='Golden', bracket=bracket)

        else:
            if log: print("Finding maximal S_z value.")
            if bracket is None: bracket = (1.0, 2.0)

            result = scipy.optimize.minimize_scalar(lambda mu: -f(mu), method='Golden', bracket=bracket)

        if log: print(result)

        if result.success == False: raise RuntimeError("The scipy minimizer (golden) did not converge.")

        mu_optimized = result.x

    else:
        if bracket is None: raise ValueError("When not targeting an extremal spin expectation value, a bracket should be given.")

        result = scipy.optimize.root_scalar(f, method=method, bracket=bracket, xtol=tolerance)

        if log: print(result)

        if result.converged == False: raise RuntimeError("The scipy root finder {} did not converge.".format(method))

        mu_optimized = result.root

    # Do the modified GHF calculation again to check if the optimized value actually delivers a valid state.
    if log: print("Using optimized mu {}".format(mu_optimized))
    mu_optimized_vector = mu_optimized * np.array([0, 0, 1])
    (E_constrained, GHF_parameters) = do_magnetic_GHF_calculation(molecule, basis_set, B, threshold=threshold, max_iter=max_iter, mu=mu_optimized_vector, log=log, stability=stability)

    spinor_basis = gqcpy.LondonGSpinorBasis(molecule, basis_set, B)
    S_AO = spinor_basis.overlap()

    S_calculated_vector = GHF_parameters.calculateExpectationValueOf(gqcpy.ElectronicSpinOperator(), S_AO)
    S_z_calculated = S_calculated_vector[2].real

    if np.abs(S_z_calculated - S_z_target) > np.sqrt(tolerance):
        raise RuntimeError("The constraint could not be met. Target: {}. Actual: {}.".format(S_z_target, S_z_calculated))

    return E_constrained, GHF_parameters, mu_optimized_vector


def find_FCI_spin_transition_bracket(molecule, basis_set, B, bracket, threshold=1.0e-08, log=False):
    """
        Attempt to find a FCI state with the requested spin expectation value, using a modified (mu) algorithm.

        @param bracket                      A bracket (tuple) in which the target multiplier (mu) is to be found.

        @param threshold                    The threshold for bracket convergence.
        @param log                          A boolean indicating if logging output should be given.

        @return:
            - The constrained FCI energy.
            - The constrained FCI state.
            - The MO expansion of the basis in which the FCI is done.
            - The bracket in-between which the spin transition occurs.
    """

    def f(mu):
        """
            The function that shows Heaviside-behavior.
                f(mu) = <S_z>(mu)
        """

        if log: print("\tEvaluating f(mu) with mu = {}".format(mu))

        # Do a FCI calculation with the modification dictated by the Lagrange multiplier.
        mu_vector = mu * np.array([0, 0, 1])
        FCI_energy, FCI_parameters, C = do_magnetic_FCI_calculation(molecule, basis_set, B, mu=mu_vector, n=1, basis="Lowdin")

        # Calculate the S_z expectation value for the modified FCI parameters.
        spinor_basis = gqcpy.LondonGSpinorBasis(molecule, basis_set, B)
        spinor_basis.transform(C)
        Spin_op = spinor_basis.quantize(gqcpy.ElectronicSpinOperator())

        D = FCI_parameters.calculate1DM()

        S_calculated_vector = np.array(Spin_op.calculateExpectationValue(D)).real  # Omit the imaginary component, since it's zero.
        S_z_calculated = S_calculated_vector[2].real

        if log: print("\tFound S_z expectation value of {}".format(S_z_calculated))

        return S_z_calculated


    # Determine the bracket in-between which the FCI spin transition occurs, using a bisection algorithm.
    a, b = bisection(f, bracket, threshold=threshold, log=log)
    if log: print("\tFound bracket: ({}, {})".format(a,b))



    # Do the modified FCI calculations again to check if the bracket actually contains a spin transition.
    S_z_a, S_z_b = f(a), f(b)
    if abs(S_z_a - S_z_b) < 1.0e-08: raise RuntimeError("The algorithm found a bracket that doesn't contain a spin transition.")


    return a, b
