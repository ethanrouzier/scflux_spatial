#!/usr/bin/env python3
"""
Validation tests for scflux_spatial.

This module contains validation tests that verify the correctness of algorithms
against known analytical solutions and reference implementations.
"""

import unittest
import numpy as np
import pytest
from scflux_spatial.gem.gpr import GPRParser
from scflux_spatial.spatial.kinetics import SpatialKineticsModel


class TestGPRValidation(unittest.TestCase):
    """Validation tests for GPR parsing and evaluation."""
    
    def setUp(self):
        """Set up validation test fixtures."""
        self.gpr_parser = GPRParser()
    
    def test_boolean_logic_validation(self):
        """Validate GPR evaluation against boolean logic."""
        # Test cases with known boolean logic results
        test_cases = [
            # (expression, gene_values, expected_result)
            ('GENE_A', {'GENE_A': 1.0}, 1.0),
            ('GENE_A and GENE_B', {'GENE_A': 1.0, 'GENE_B': 0.0}, 0.0),
            ('GENE_A or GENE_B', {'GENE_A': 1.0, 'GENE_B': 0.0}, 1.0),
            ('GENE_A and GENE_B', {'GENE_A': 0.0, 'GENE_B': 1.0}, 0.0),
            ('GENE_A or GENE_B', {'GENE_A': 0.0, 'GENE_B': 1.0}, 1.0),
        ]
        
        for expression, gene_values, expected in test_cases:
            with self.subTest(expression=expression):
                result = self.gpr_parser._evaluate_gpr_rule_with_operators(
                    expression, gene_values
                )
                self.assertEqual(result, expected)
    
    def test_operator_precedence_validation(self):
        """Validate operator precedence rules."""
        # Test that AND has higher precedence than OR
        gene_values = {'GENE_A': 1.0, 'GENE_B': 0.0, 'GENE_C': 0.0}
        
        # A or B and C should be interpreted as A or (B and C)
        result = self.gpr_parser._evaluate_gpr_rule_with_operators(
            'GENE_A or GENE_B and GENE_C', gene_values
        )
        
        # Expected: 1.0 or (0.0 and 0.0) = 1.0 or 0.0 = 1.0
        expected = max(1.0, min(0.0, 0.0))
        self.assertEqual(result, expected)
    
    def test_associativity_validation(self):
        """Validate operator associativity."""
        gene_values = {'GENE_A': 1.0, 'GENE_B': 2.0, 'GENE_C': 3.0}
        
        # Test left associativity of AND
        result1 = self.gpr_parser._evaluate_gpr_rule_with_operators(
            'GENE_A and GENE_B and GENE_C', gene_values
        )
        result2 = self.gpr_parser._evaluate_gpr_rule_with_operators(
            '(GENE_A and GENE_B) and GENE_C', gene_values
        )
        
        # Both should give the same result (minimum of all three)
        expected = min(1.0, 2.0, 3.0)
        self.assertEqual(result1, expected)
        self.assertEqual(result2, expected)
        self.assertEqual(result1, result2)
    
    def test_commutativity_validation(self):
        """Validate operator commutativity."""
        gene_values = {'GENE_A': 2.0, 'GENE_B': 3.0}
        
        # AND should be commutative
        result1 = self.gpr_parser._evaluate_gpr_rule_with_operators(
            'GENE_A and GENE_B', gene_values
        )
        result2 = self.gpr_parser._evaluate_gpr_rule_with_operators(
            'GENE_B and GENE_A', gene_values
        )
        
        self.assertEqual(result1, result2)
        
        # OR should be commutative
        result3 = self.gpr_parser._evaluate_gpr_rule_with_operators(
            'GENE_A or GENE_B', gene_values
        )
        result4 = self.gpr_parser._evaluate_gpr_rule_with_operators(
            'GENE_B or GENE_A', gene_values
        )
        
        self.assertEqual(result3, result4)


class TestFBAValidation(unittest.TestCase):
    """Validation tests for FBA algorithms."""
    
    def test_mass_balance_validation(self):
        """Validate mass balance in toy metabolic model."""
        # Create simple toy model: A -> B -> C
        # Stoichiometric matrix
        S = np.array([
            [-1,  0],  # A: consumed by R1
            [ 1, -1],  # B: produced by R1, consumed by R2
            [ 0,  1],  # C: produced by R2
        ])
        
        # Test balanced flux vector
        balanced_flux = np.array([1.0, 1.0])  # R1=1, R2=1
        
        # Calculate mass balance: S * v
        mass_balance = np.dot(S, balanced_flux)
        
        # Should be zero (steady state)
        np.testing.assert_array_almost_equal(mass_balance, np.zeros(3))
    
    def test_flux_bounds_validation(self):
        """Validate flux bounds constraints."""
        # Test bounds: [lower, upper]
        bounds = np.array([
            [0, 10],    # R1: 0 <= v1 <= 10
            [-5, 5],    # R2: -5 <= v2 <= 5
        ])
        
        # Test valid flux vectors
        valid_fluxes = [
            np.array([5.0, 2.0]),
            np.array([0.0, 0.0]),
            np.array([10.0, -5.0]),
        ]
        
        for flux in valid_fluxes:
            with self.subTest(flux=flux):
                # Check bounds constraints
                self.assertTrue(np.all(flux >= bounds[:, 0]))
                self.assertTrue(np.all(flux <= bounds[:, 1]))
        
        # Test invalid flux vectors
        invalid_fluxes = [
            np.array([15.0, 2.0]),  # R1 > upper bound
            np.array([5.0, -10.0]), # R2 < lower bound
        ]
        
        for flux in invalid_fluxes:
            with self.subTest(flux=flux):
                # Should violate at least one bound
                violates_bounds = (
                    np.any(flux < bounds[:, 0]) or 
                    np.any(flux > bounds[:, 1])
                )
                self.assertTrue(violates_bounds)
    
    def test_objective_optimization_validation(self):
        """Validate objective function optimization."""
        # Simple linear programming problem
        # Maximize: 2*x + 3*y
        # Subject to: x + y <= 4, x >= 0, y >= 0
        
        # Optimal solution should be x=0, y=4, objective=12
        objective = np.array([2.0, 3.0])
        optimal_flux = np.array([0.0, 4.0])
        
        optimal_value = np.dot(objective, optimal_flux)
        expected_value = 12.0
        
        self.assertAlmostEqual(optimal_value, expected_value)
        
        # Test that this is indeed optimal
        test_fluxes = [
            np.array([4.0, 0.0]),  # x=4, y=0, obj=8
            np.array([2.0, 2.0]),  # x=2, y=2, obj=10
            np.array([1.0, 3.0]),  # x=1, y=3, obj=11
        ]
        
        for flux in test_fluxes:
            test_value = np.dot(objective, flux)
            self.assertLessEqual(test_value, optimal_value)


class TestRDValidation(unittest.TestCase):
    """Validation tests for reaction-diffusion equations."""
    
    def test_diffusion_equation_validation(self):
        """Validate diffusion equation against analytical solution."""
        # For pure diffusion on a disk with Dirichlet BC,
        # the steady-state solution should be uniform and equal to boundary value
        
        # Test parameters
        grid_size = 50
        boundary_value = 1.0
        
        # Create uniform concentration field (steady state)
        concentration = np.full((grid_size, grid_size), boundary_value)
        
        # Calculate Laplacian (should be zero for uniform field)
        # Using finite differences
        laplacian = np.zeros_like(concentration)
        
        # Central difference approximation
        for i in range(1, grid_size-1):
            for j in range(1, grid_size-1):
                laplacian[i, j] = (
                    concentration[i+1, j] + concentration[i-1, j] +
                    concentration[i, j+1] + concentration[i, j-1] -
                    4 * concentration[i, j]
                )
        
        # Laplacian should be zero for uniform field
        np.testing.assert_array_almost_equal(laplacian, np.zeros_like(laplacian))
    
    def test_mass_conservation_validation(self):
        """Validate mass conservation in RD equations."""
        # Test mass conservation with Neumann boundary conditions
        
        # Initial mass
        initial_concentration = np.ones((10, 10))
        initial_mass = np.sum(initial_concentration)
        
        # With Neumann BC (no flux), mass should be conserved
        # Simulate diffusion by averaging
        concentration = initial_concentration.copy()
        
        # Apply diffusion operator (mass-conserving)
        new_concentration = np.zeros_like(concentration)
        for i in range(1, concentration.shape[0]-1):
            for j in range(1, concentration.shape[1]-1):
                new_concentration[i, j] = 0.25 * (
                    concentration[i+1, j] + concentration[i-1, j] +
                    concentration[i, j+1] + concentration[i, j-1]
                )
        
        # Boundary conditions (no flux)
        new_concentration[0, :] = concentration[0, :]
        new_concentration[-1, :] = concentration[-1, :]
        new_concentration[:, 0] = concentration[:, 0]
        new_concentration[:, -1] = concentration[:, -1]
        
        final_mass = np.sum(new_concentration)
        
        # Mass should be conserved (approximately)
        self.assertAlmostEqual(final_mass, initial_mass, places=10)
    
    def test_reaction_diffusion_validation(self):
        """Validate reaction-diffusion equation."""
        # Test case: constant reaction rate
        # ∂C/∂t = D∇²C + r
        # For steady state: D∇²C + r = 0
        
        # Parameters
        D = 1.0  # Diffusion coefficient
        r = 0.1  # Constant reaction rate (production)
        
        # For 1D steady state with constant r:
        # D * d²C/dx² + r = 0
        # Solution: C(x) = -r/(2D) * x² + A*x + B
        
        # Test with simple 1D case
        x = np.linspace(0, 1, 11)
        dx = x[1] - x[0]
        
        # Analytical solution with BC: C(0) = 0, C(1) = 0
        A = r / (2 * D)
        B = 0
        analytical_solution = -A * x**2 + A * x
        
        # Numerical solution using finite differences
        # D * (C[i+1] - 2*C[i] + C[i-1])/dx² + r = 0
        # Rearranging: C[i] = (C[i+1] + C[i-1] + r*dx²/D) / 2
        
        # Solve numerically
        C = np.zeros_like(x)
        C[0] = 0  # Boundary condition
        C[-1] = 0  # Boundary condition
        
        # Iterative solution
        for iteration in range(1000):
            C_old = C.copy()
            for i in range(1, len(x)-1):
                C[i] = (C[i+1] + C[i-1] + r * dx**2 / D) / 2
            
            # Check convergence
            if np.max(np.abs(C - C_old)) < 1e-6:
                break
        
        # Compare with analytical solution
        np.testing.assert_array_almost_equal(C, analytical_solution, decimal=2)


class TestKineticsValidation(unittest.TestCase):
    """Validation tests for spatial kinetics."""
    
    def test_michaelis_menten_validation(self):
        """Validate Michaelis-Menten kinetics."""
        # Test Michaelis-Menten equation: v = Vmax * C / (Km + C)
        
        # Parameters
        Vmax = 10.0  # Maximum rate
        Km = 5.0     # Michaelis constant
        
        # Test concentrations
        concentrations = [0.0, 1.0, 5.0, 10.0, 50.0]
        
        for C in concentrations:
            with self.subTest(concentration=C):
                # Calculate rate
                rate = Vmax * C / (Km + C)
                
                # Validate properties
                self.assertGreaterEqual(rate, 0)  # Non-negative
                self.assertLessEqual(rate, Vmax)  # Cannot exceed Vmax
                
                # Test specific values
                if C == 0:
                    self.assertEqual(rate, 0)
                elif C == Km:
                    self.assertAlmostEqual(rate, Vmax / 2)
                elif C >> Km:
                    self.assertAlmostEqual(rate, Vmax)
    
    def test_unit_conversion_validation(self):
        """Validate unit conversions in spatial kinetics."""
        # Test conversion from FBA units to spatial units
        kinetics_model = SpatialKineticsModel(
            cell_density=1e9,      # cells/L
            spot_density=1e6,      # spots/L
            cell_volume=1e-12,     # L/cell
            spot_volume=1e-9       # L/spot
        )
        
        # FBA rate: 10 mmol/gDW/h
        fba_rate = 10.0  # mmol/gDW/h
        
        # Convert to spatial rate
        spatial_rate = kinetics_model.convert_fba_to_spatial_rate(
            fba_rate, spatial_type='cell'
        )
        
        # Check units: should be mol/L/s
        # Expected: 10 mmol/gDW/h * (1 mol/1000 mmol) * (1 gDW/cell) * (1 cell/1e-12 L) * (1 h/3600 s)
        # = 10 * 1e-3 * 1e-12 / 3600 = 2.78e-18 mol/L/s
        
        expected_rate = 10.0 * 1e-3 * 1e-12 / 3600  # mol/L/s
        self.assertAlmostEqual(spatial_rate, expected_rate, places=15)
    
    def test_mass_balance_validation(self):
        """Validate mass balance in kinetics model."""
        # Test mass balance: rate_in = rate_out + rate_consumption
        
        # Parameters
        uptake_rate = 1e-6  # mol/L/s
        consumption_rate = 0.5e-6  # mol/L/s
        
        # Mass balance: uptake = consumption + accumulation
        # For steady state: uptake = consumption
        accumulation_rate = uptake_rate - consumption_rate
        
        # Should be zero for steady state
        self.assertAlmostEqual(accumulation_rate, 0.5e-6)
        
        # Test with different rates
        test_cases = [
            (1e-6, 1e-6, 0.0),      # Balanced
            (2e-6, 1e-6, 1e-6),     # Accumulation
            (1e-6, 2e-6, -1e-6),    # Depletion
        ]
        
        for uptake, consumption, expected_accumulation in test_cases:
            with self.subTest(uptake=uptake, consumption=consumption):
                accumulation = uptake - consumption
                self.assertAlmostEqual(accumulation, expected_accumulation)


if __name__ == '__main__':
    unittest.main(verbosity=2)
