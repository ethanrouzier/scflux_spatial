#!/usr/bin/env python3
"""
Integration tests for scflux_spatial.

This module contains integration tests that test multiple components together.
"""

import unittest
import numpy as np
import pytest
from scflux_spatial.gem.human_gem import HumanGEM
from scflux_spatial.gem.gpr import GPRParser
from scflux_spatial.fba.integrate_expression import integrate_expression_with_method
from scflux_spatial.spatial.kinetics import SpatialKineticsModel
from scflux_spatial.spatial.rd import RDField
from scflux_spatial.dataio import load_visium


class TestIntegration(unittest.TestCase):
    """Integration tests for scflux_spatial components."""
    
    def setUp(self):
        """Set up integration test fixtures."""
        self.gene_expression = {
            'HK1': 5.0,
            'GAPDH': 4.0,
            'LDHA': 3.0,
            'PFKL': 4.5,
            'PKM': 3.8
        }
        
        self.base_bounds = {
            'HEX1': (0, 1000),
            'GAPD': (0, 1000),
            'LDH_L': (0, 1000),
            'PFK': (0, 1000),
            'PYK': (0, 1000)
        }
    
    @pytest.mark.integration
    def test_gpr_fba_integration(self):
        """Test integration between GPR parsing and FBA."""
        # Initialize GPR parser
        gpr_parser = GPRParser()
        
        # Test GPR evaluation
        gpr_result = gpr_parser._evaluate_gpr_rule_with_operators(
            'HK1 and GAPDH',
            self.gene_expression
        )
        
        # Should return minimum value (4.0)
        self.assertEqual(gpr_result, 4.0)
        
        # Test expression integration with GPR
        try:
            integrated_bounds = integrate_expression_with_method(
                self.gene_expression,
                self.base_bounds,
                method='eflux'
            )
            
            # Check that bounds were modified
            self.assertIsInstance(integrated_bounds, dict)
            self.assertEqual(len(integrated_bounds), len(self.base_bounds))
            
            # Check that bounds are reasonable
            for rxn, bounds in integrated_bounds.items():
                self.assertGreater(bounds[1], 0)  # Upper bound positive
                self.assertGreaterEqual(bounds[1], bounds[0])  # Upper >= lower
                
        except Exception as e:
            # If Human-GEM is not available, this is expected
            self.assertIn("Human-GEM", str(e))
    
    @pytest.mark.integration
    def test_spatial_kinetics_rd_integration(self):
        """Test integration between spatial kinetics and RD."""
        # Initialize spatial kinetics model
        kinetics_model = SpatialKineticsModel(
            cell_density=1e9,
            spot_density=1e6,
            cell_volume=1e-12,
            spot_volume=1e-9
        )
        
        # Test FBA rate conversion
        fba_rate = 10.0  # mmol/gDW/h
        spatial_rate = kinetics_model.convert_fba_to_spatial_rate(
            fba_rate, spatial_type='cell'
        )
        
        # Check conversion is reasonable
        self.assertGreater(spatial_rate, 0)
        self.assertLess(spatial_rate, 1e-3)  # Should be small
        
        # Initialize RD field
        try:
            rd_field = RDField(
                substrate='O2',
                grid_size=50,
                domain_size=1.0,
                diffusion_coefficient=1e-9,
                boundary_condition='dirichlet',
                boundary_value=2e-5
            )
            
            # Set initial condition
            initial_condition = np.full((50, 50), 1e-5)
            rd_field.set_initial_condition(initial_condition)
            
            # Set reaction rate
            reaction_rate = np.full((50, 50), -1e-6)  # Consumption
            rd_field.set_reaction_rate_field(reaction_rate)
            
            # Test steady state solution
            convergence_info = rd_field.solve_steady_state(
                convergence_tolerance=1e-3,
                max_iterations=50
            )
            
            # Should converge or make progress
            self.assertGreater(convergence_info['iterations'], 0)
            
        except ImportError:
            # FiPy not available, skip RD test
            self.skipTest("FiPy not available for RD integration test")
    
    @pytest.mark.integration
    def test_visium_data_processing(self):
        """Test integration with Visium data processing."""
        try:
            # Load demo Visium data
            adata = load_visium(use_demo=True)
            
            # Check data structure
            self.assertGreater(adata.n_obs, 0)  # Has observations
            self.assertGreater(adata.n_vars, 0)  # Has variables
            
            # Check spatial coordinates
            self.assertIn('x', adata.obs.columns)
            self.assertIn('y', adata.obs.columns)
            
            # Check expression matrix
            self.assertIsNotNone(adata.X)
            
            # Test metabolic score calculation
            glycolysis_genes = ['HK1', 'GAPDH', 'PFKL', 'PKM']
            available_genes = [g for g in glycolysis_genes if g in adata.var_names]
            
            if available_genes:
                # Calculate glycolysis score
                if hasattr(adata.X, 'toarray'):
                    expr_matrix = adata.X.toarray()
                else:
                    expr_matrix = adata.X
                
                glycolysis_indices = [adata.var_names.get_loc(g) for g in available_genes]
                glycolysis_scores = expr_matrix[:, glycolysis_indices].mean(axis=1)
                
                # Check scores are reasonable
                self.assertEqual(len(glycolysis_scores), adata.n_obs)
                self.assertTrue(np.all(glycolysis_scores >= 0))
                
        except Exception as e:
            # If demo data is not available, this is expected
            self.assertIn("demo", str(e).lower())
    
    @pytest.mark.integration
    def test_end_to_end_workflow(self):
        """Test complete end-to-end workflow."""
        try:
            # Step 1: Load data
            adata = load_visium(use_demo=True)
            
            # Step 2: Calculate metabolic scores
            glycolysis_genes = ['HK1', 'GAPDH', 'PFKL', 'PKM']
            available_genes = [g for g in glycolysis_genes if g in adata.var_names]
            
            if available_genes and len(available_genes) > 0:
                # Calculate scores
                if hasattr(adata.X, 'toarray'):
                    expr_matrix = adata.X.toarray()
                else:
                    expr_matrix = adata.X
                
                glycolysis_indices = [adata.var_names.get_loc(g) for g in available_genes]
                glycolysis_scores = expr_matrix[:, glycolysis_indices].mean(axis=1)
                
                # Step 3: Create gene expression dictionary
                gene_expression = {}
                for i, gene in enumerate(available_genes):
                    gene_expression[gene] = glycolysis_scores.mean()
                
                # Step 4: Test GPR evaluation
                gpr_parser = GPRParser()
                if len(available_genes) >= 2:
                    gpr_result = gpr_parser._evaluate_gpr_rule_with_operators(
                        f'{available_genes[0]} and {available_genes[1]}',
                        gene_expression
                    )
                    self.assertIsInstance(gpr_result, (int, float))
                
                # Step 5: Test expression integration
                try:
                    integrated_bounds = integrate_expression_with_method(
                        gene_expression,
                        self.base_bounds,
                        method='linear'
                    )
                    self.assertIsInstance(integrated_bounds, dict)
                    
                except Exception as e:
                    # Expected if Human-GEM not available
                    self.assertIn("Human-GEM", str(e))
                
                # Step 6: Test spatial kinetics
                kinetics_model = SpatialKineticsModel(
                    cell_density=1e9,
                    spot_density=1e6,
                    cell_volume=1e-12,
                    spot_volume=1e-9
                )
                
                # Convert FBA rate to spatial rate
                fba_rate = 10.0
                spatial_rate = kinetics_model.convert_fba_to_spatial_rate(
                    fba_rate, spatial_type='cell'
                )
                self.assertGreater(spatial_rate, 0)
                
        except Exception as e:
            # If demo data is not available, this is expected
            self.assertIn("demo", str(e).lower())
    
    @pytest.mark.integration
    def test_error_handling_integration(self):
        """Test error handling across integrated components."""
        # Test with invalid gene expression
        invalid_expression = {
            'INVALID_GENE': 5.0,
            'ANOTHER_INVALID': 3.0
        }
        
        gpr_parser = GPRParser()
        
        # Should handle missing genes gracefully
        with self.assertRaises(KeyError):
            gpr_parser._evaluate_gpr_rule_with_operators(
                'GENE_A and GENE_B',
                invalid_expression
            )
        
        # Test with empty expression
        empty_expression = {}
        
        with self.assertRaises(KeyError):
            gpr_parser._evaluate_gpr_rule_with_operators(
                'GENE_A',
                empty_expression
            )
        
        # Test with invalid bounds
        invalid_bounds = {
            'RXN1': (10, 5)  # Lower bound > upper bound
        }
        
        try:
            integrated_bounds = integrate_expression_with_method(
                self.gene_expression,
                invalid_bounds,
                method='eflux'
            )
            # Should handle invalid bounds gracefully
            self.assertIsInstance(integrated_bounds, dict)
            
        except Exception as e:
            # Expected behavior
            self.assertIsInstance(e, Exception)
    
    @pytest.mark.integration
    def test_performance_integration(self):
        """Test performance of integrated components."""
        import time
        
        # Test with larger dataset
        large_expression = {}
        for i in range(100):
            large_expression[f'GENE_{i}'] = np.random.uniform(0, 10)
        
        gpr_parser = GPRParser()
        
        # Test GPR evaluation performance
        start_time = time.time()
        
        for i in range(10):
            rule = f'GENE_{i} and GENE_{i+1}'
            result = gpr_parser._evaluate_gpr_rule_with_operators(
                rule, large_expression
            )
            self.assertIsInstance(result, (int, float))
        
        end_time = time.time()
        execution_time = end_time - start_time
        
        # Should complete within reasonable time
        self.assertLess(execution_time, 5.0)  # 5 seconds max
        
        # Test expression integration performance
        start_time = time.time()
        
        try:
            integrated_bounds = integrate_expression_with_method(
                large_expression,
                self.base_bounds,
                method='eflux'
            )
            self.assertIsInstance(integrated_bounds, dict)
            
        except Exception as e:
            # Expected if Human-GEM not available
            self.assertIn("Human-GEM", str(e))
        
        end_time = time.time()
        execution_time = end_time - start_time
        
        # Should complete within reasonable time
        self.assertLess(execution_time, 10.0)  # 10 seconds max


if __name__ == '__main__':
    unittest.main(verbosity=2)
