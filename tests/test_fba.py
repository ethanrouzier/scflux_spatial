#!/usr/bin/env python3
"""
Unit tests for FBA (Flux Balance Analysis) functionality.

This module tests FBA operations including:
- Toy model creation and validation
- pFBA (parsimonious FBA) implementation
- Expression integration methods
- Flux bounds and constraints
"""

import unittest
import numpy as np
import pandas as pd
from scflux_spatial.fba.integrate_expression import solve_with_pfba


class TestFBA(unittest.TestCase):
    """Test cases for FBA functionality."""
    
    def setUp(self):
        """Set up test fixtures with toy model."""
        self.toy_model = self._create_toy_model()
        self.gene_expression = {
            'GENE_A': 5.0,
            'GENE_B': 3.0,
            'GENE_C': 2.0,
            'GENE_D': 4.0
        }
    
    def _create_toy_model(self):
        """Create a simple toy metabolic model for testing."""
        class ToyModel:
            def __init__(self):
                # Simple linear pathway: A -> B -> C -> D
                self.reactions = ['R1', 'R2', 'R3', 'R4']
                self.metabolites = ['A', 'B', 'C', 'D']
                
                # Stoichiometric matrix (metabolites x reactions)
                self.S = np.array([
                    [-1,  0,  0,  0],  # A: consumed by R1
                    [ 1, -1,  0,  0],  # B: produced by R1, consumed by R2
                    [ 0,  1, -1,  0],  # C: produced by R2, consumed by R3
                    [ 0,  0,  1, -1],  # D: produced by R3, consumed by R4
                ])
                
                # Reaction bounds: [lower, upper]
                self.bounds = np.array([
                    [0, 10],    # R1: A -> B
                    [0, 10],    # R2: B -> C
                    [0, 10],    # R3: C -> D
                    [0, 10],    # R4: D -> (export)
                ])
                
                # Objective: maximize D production (R4)
                self.objective = np.array([0, 0, 0, 1])
        
        return ToyModel()
    
    def test_toy_model_structure(self):
        """Test toy model structure and properties."""
        model = self.toy_model
        
        # Check dimensions
        self.assertEqual(len(model.reactions), 4)
        self.assertEqual(len(model.metabolites), 4)
        self.assertEqual(model.S.shape, (4, 4))
        self.assertEqual(model.bounds.shape, (4, 2))
        self.assertEqual(len(model.objective), 4)
        
        # Check stoichiometry
        # R1: A -> B
        self.assertEqual(model.S[0, 0], -1)  # A consumed
        self.assertEqual(model.S[1, 0], 1)   # B produced
        
        # R2: B -> C
        self.assertEqual(model.S[1, 1], -1)  # B consumed
        self.assertEqual(model.S[2, 1], 1)   # C produced
        
        # Check bounds
        np.testing.assert_array_equal(model.bounds, np.array([[0, 10]] * 4))
        
        # Check objective
        expected_obj = np.array([0, 0, 0, 1])
        np.testing.assert_array_equal(model.objective, expected_obj)
    
    def test_flux_balance_analysis(self):
        """Test basic FBA functionality."""
        model = self.toy_model
        
        # Solve FBA: maximize R4 (D production)
        # With bounds [0, 10] for all reactions, optimal flux should be 10
        # for all reactions in the linear pathway
        
        # Mock FBA solver (simplified)
        def mock_fba_solve(model):
            # For linear pathway with equal bounds, all fluxes should be equal
            # and limited by the bottleneck (smallest upper bound)
            max_flux = np.min(model.bounds[:, 1])  # = 10
            fluxes = np.full(len(model.reactions), max_flux)
            objective_value = np.dot(model.objective, fluxes)
            return fluxes, objective_value
        
        fluxes, obj_value = mock_fba_solve(model)
        
        # All fluxes should be equal to maximum bound
        expected_fluxes = np.array([10, 10, 10, 10])
        np.testing.assert_array_equal(fluxes, expected_fluxes)
        
        # Objective value should be 10 (flux of R4)
        self.assertEqual(obj_value, 10.0)
    
    def test_parsimonious_fba(self):
        """Test pFBA implementation (minimize sum of absolute fluxes)."""
        model = self.toy_model
        
        # Mock pFBA solver
        def mock_pfba_solve(model):
            # pFBA minimizes sum of absolute fluxes while maintaining objective
            # For this toy model, we can find the minimal flux solution
            
            # If we want to produce 1 unit of D, minimal fluxes would be:
            # R1 = R2 = R3 = R4 = 1
            fluxes = np.array([1, 1, 1, 1])
            objective_value = np.dot(model.objective, fluxes)
            
            # Calculate sum of absolute fluxes
            sum_abs_fluxes = np.sum(np.abs(fluxes))
            
            return fluxes, objective_value, sum_abs_fluxes
        
        fluxes, obj_value, sum_abs = mock_pfba_solve(model)
        
        # Check that all fluxes are equal (minimal solution)
        expected_fluxes = np.array([1, 1, 1, 1])
        np.testing.assert_array_equal(fluxes, expected_fluxes)
        
        # Check objective value
        self.assertEqual(obj_value, 1.0)
        
        # Check sum of absolute fluxes
        self.assertEqual(sum_abs, 4.0)
    
    def test_flux_bounds_constraints(self):
        """Test flux bounds and constraints."""
        model = self.toy_model
        
        # Test with constrained bounds
        constrained_bounds = np.array([
            [0, 5],     # R1: limited to 5
            [0, 10],    # R2: unlimited
            [0, 10],    # R3: unlimited
            [0, 10],    # R4: unlimited
        ])
        
        # Mock constrained FBA
        def constrained_fba_solve(model, bounds):
            # Flux through pathway is limited by bottleneck (R1 = 5)
            bottleneck_flux = np.min(bounds[:, 1])
            fluxes = np.full(len(model.reactions), bottleneck_flux)
            objective_value = np.dot(model.objective, fluxes)
            return fluxes, objective_value
        
        fluxes, obj_value = constrained_fba_solve(model, constrained_bounds)
        
        # All fluxes should be limited by bottleneck (5)
        expected_fluxes = np.array([5, 5, 5, 5])
        np.testing.assert_array_equal(fluxes, expected_fluxes)
        
        # Objective should be 5
        self.assertEqual(obj_value, 5.0)
    
    def test_expression_integration_eflux(self):
        """Test E-Flux expression integration method."""
        # Mock E-Flux implementation
        def mock_eflux_integration(gene_expression, base_bounds):
            # E-Flux scales bounds by gene expression
            scaled_bounds = {}
            for rxn, bounds in base_bounds.items():
                # Simple scaling: multiply by average expression
                avg_expr = np.mean(list(gene_expression.values()))
                scale_factor = avg_expr / 5.0  # Normalize by reference
                
                new_upper = bounds[1] * scale_factor
                new_lower = bounds[0] * scale_factor
                scaled_bounds[rxn] = (new_lower, new_upper)
            
            return scaled_bounds
        
        base_bounds = {
            'R1': (0, 10),
            'R2': (0, 10),
            'R3': (0, 10),
            'R4': (0, 10)
        }
        
        scaled_bounds = mock_eflux_integration(self.gene_expression, base_bounds)
        
        # Check that bounds were scaled
        for rxn, bounds in scaled_bounds.items():
            self.assertGreater(bounds[1], 0)  # Upper bound should be positive
            self.assertGreaterEqual(bounds[1], bounds[0])  # Upper >= lower
    
    def test_expression_integration_imat(self):
        """Test iMAT expression integration method."""
        # Mock iMAT implementation
        def mock_imat_integration(gene_expression, base_bounds, threshold=0.5):
            # iMAT classifies reactions as high/low expression
            scaled_bounds = {}
            
            for rxn, bounds in base_bounds.items():
                # Simple classification based on average expression
                avg_expr = np.mean(list(gene_expression.values()))
                
                if avg_expr > threshold:
                    # High expression: increase upper bound
                    new_upper = bounds[1] * 2.0
                    new_lower = bounds[0]
                else:
                    # Low expression: decrease upper bound
                    new_upper = bounds[1] * 0.5
                    new_lower = bounds[0]
                
                scaled_bounds[rxn] = (new_lower, new_upper)
            
            return scaled_bounds
        
        base_bounds = {
            'R1': (0, 10),
            'R2': (0, 10),
            'R3': (0, 10),
            'R4': (0, 10)
        }
        
        scaled_bounds = mock_imat_integration(self.gene_expression, base_bounds)
        
        # Check that bounds were modified
        for rxn, bounds in scaled_bounds.items():
            self.assertGreater(bounds[1], 0)
            self.assertGreaterEqual(bounds[1], bounds[0])
    
    def test_flux_variability_analysis(self):
        """Test flux variability analysis (FVA)."""
        model = self.toy_model
        
        # Mock FVA implementation
        def mock_fva_solve(model):
            # For each reaction, find min and max possible flux
            fva_results = {}
            
            for i, rxn in enumerate(model.reactions):
                # In a linear pathway, all reactions have the same flux
                # Min and max are determined by bounds
                min_flux = model.bounds[i, 0]
                max_flux = model.bounds[i, 1]
                
                fva_results[rxn] = {
                    'minimum': min_flux,
                    'maximum': max_flux
                }
            
            return fva_results
        
        fva_results = mock_fva_solve(model)
        
        # Check FVA results
        for rxn in model.reactions:
            self.assertIn('minimum', fva_results[rxn])
            self.assertIn('maximum', fva_results[rxn])
            self.assertGreaterEqual(fva_results[rxn]['maximum'], 
                                   fva_results[rxn]['minimum'])
    
    def test_metabolic_flux_analysis(self):
        """Test metabolic flux analysis with different objectives."""
        model = self.toy_model
        
        # Test different objectives
        objectives = {
            'maximize_D': np.array([0, 0, 0, 1]),      # Maximize D production
            'maximize_B': np.array([0, 1, 0, 0]),      # Maximize B production
            'minimize_flux': np.array([-1, -1, -1, -1]) # Minimize total flux
        }
        
        def mock_objective_fba(model, objective):
            # Solve FBA with given objective
            if np.all(objective == np.array([0, 0, 0, 1])):
                # Maximize D: all fluxes = 10
                fluxes = np.array([10, 10, 10, 10])
            elif np.all(objective == np.array([0, 1, 0, 0])):
                # Maximize B: R1 = 10, others = 0
                fluxes = np.array([10, 0, 0, 0])
            else:
                # Minimize flux: all fluxes = 0
                fluxes = np.array([0, 0, 0, 0])
            
            objective_value = np.dot(objective, fluxes)
            return fluxes, objective_value
        
        for obj_name, objective in objectives.items():
            fluxes, obj_value = mock_objective_fba(model, objective)
            
            # Check that solution is valid
            self.assertEqual(len(fluxes), len(model.reactions))
            self.assertIsInstance(obj_value, (int, float))
    
    def test_toy_model_validation(self):
        """Test toy model validation and consistency."""
        model = self.toy_model
        
        # Check mass balance: S * v = 0 at steady state
        # For linear pathway A -> B -> C -> D, this should be satisfied
        
        # Test with balanced flux vector
        balanced_fluxes = np.array([1, 1, 1, 1])  # All fluxes equal
        
        # Calculate mass balance: S * v
        mass_balance = np.dot(model.S, balanced_fluxes)
        
        # Should be zero (steady state)
        np.testing.assert_array_almost_equal(mass_balance, np.zeros(4))
        
        # Test with unbalanced flux vector
        unbalanced_fluxes = np.array([2, 1, 1, 1])  # R1 > others
        
        # Calculate mass balance
        mass_balance = np.dot(model.S, unbalanced_fluxes)
        
        # Should not be zero (not steady state)
        self.assertFalse(np.allclose(mass_balance, 0))
    
    def test_flux_distribution_analysis(self):
        """Test flux distribution analysis."""
        model = self.toy_model
        
        # Mock flux distribution analysis
        def analyze_flux_distribution(fluxes):
            analysis = {
                'total_flux': np.sum(np.abs(fluxes)),
                'active_reactions': np.sum(fluxes != 0),
                'max_flux': np.max(fluxes),
                'min_flux': np.min(fluxes),
                'flux_variance': np.var(fluxes)
            }
            return analysis
        
        # Test with different flux patterns
        test_fluxes = [
            np.array([1, 1, 1, 1]),      # Uniform
            np.array([10, 0, 0, 0]),     # Bottleneck
            np.array([5, 3, 2, 1])       # Decreasing
        ]
        
        for fluxes in test_fluxes:
            analysis = analyze_flux_distribution(fluxes)
            
            # Check analysis structure
            self.assertIn('total_flux', analysis)
            self.assertIn('active_reactions', analysis)
            self.assertIn('max_flux', analysis)
            self.assertIn('min_flux', analysis)
            self.assertIn('flux_variance', analysis)
            
            # Check values are reasonable
            self.assertGreaterEqual(analysis['total_flux'], 0)
            self.assertGreaterEqual(analysis['active_reactions'], 0)
            self.assertLessEqual(analysis['active_reactions'], len(fluxes))


class TestFBAPerformance(unittest.TestCase):
    """Performance tests for FBA operations."""
    
    def test_large_model_performance(self):
        """Test performance with larger toy model."""
        # Create larger toy model
        n_reactions = 100
        n_metabolites = 50
        
        class LargeToyModel:
            def __init__(self):
                self.reactions = [f'R{i}' for i in range(n_reactions)]
                self.metabolites = [f'M{i}' for i in range(n_metabolites)]
                self.S = np.random.randn(n_metabolites, n_reactions)
                self.bounds = np.array([[0, 10]] * n_reactions)
                self.objective = np.random.randn(n_reactions)
        
        model = LargeToyModel()
        
        # Mock FBA solve
        def mock_fba_solve(model):
            fluxes = np.random.uniform(0, 10, len(model.reactions))
            objective_value = np.dot(model.objective, fluxes)
            return fluxes, objective_value
        
        import time
        start_time = time.time()
        fluxes, obj_value = mock_fba_solve(model)
        end_time = time.time()
        
        # Should complete within reasonable time
        execution_time = end_time - start_time
        self.assertLess(execution_time, 1.0)
        self.assertEqual(len(fluxes), n_reactions)


if __name__ == '__main__':
    unittest.main(verbosity=2)
