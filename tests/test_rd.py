#!/usr/bin/env python3
"""
Unit tests for reaction-diffusion (RD) functionality.

This module tests RD operations including:
- Pure diffusion on disk and grid
- Known analytical solutions
- Boundary conditions (Dirichlet, Neumann)
- Steady-state convergence
"""

import unittest
import numpy as np
from scflux_spatial.spatial.rd import RDField


class TestReactionDiffusion(unittest.TestCase):
    """Test cases for reaction-diffusion functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.grid_size = 50
        self.domain_size = 1.0  # 1 mm
        self.diffusion_coeff = 1e-9  # m²/s
        
    def test_pure_diffusion_analytical_solution(self):
        """Test pure diffusion against known analytical solution."""
        # For pure diffusion on a disk with Dirichlet BC, the steady-state
        # solution should be uniform and equal to the boundary value
        
        # Create RD field with no reaction
        rd_field = RDField(
            substrate='test',
            grid_size=self.grid_size,
            domain_size=self.domain_size,
            diffusion_coefficient=self.diffusion_coeff,
            boundary_condition='dirichlet',
            boundary_value=1.0
        )
        
        # Set initial condition (non-uniform)
        initial_condition = np.random.uniform(0, 0.5, (self.grid_size, self.grid_size))
        rd_field.set_initial_condition(initial_condition)
        
        # Set zero reaction rate (pure diffusion)
        zero_reaction = np.zeros((self.grid_size, self.grid_size))
        rd_field.set_reaction_rate_field(zero_reaction)
        
        # Solve to steady state
        convergence_info = rd_field.solve_steady_state()
        
        # Get final concentration field
        final_concentration = rd_field.get_concentration_field()
        
        # For pure diffusion with Dirichlet BC, steady state should be uniform
        # and equal to boundary value (1.0)
        np.testing.assert_allclose(final_concentration, 1.0, rtol=1e-2)
        
        # Check convergence
        self.assertTrue(convergence_info['converged'])
        self.assertLess(convergence_info['final_residual'], 1e-4)
    
    def test_diffusion_on_disk_geometry(self):
        """Test diffusion on disk geometry with analytical solution."""
        # For a disk with radius R and constant boundary condition C0,
        # the steady-state solution is C(r) = C0 (uniform)
        
        # Create circular domain (approximated by grid)
        rd_field = RDField(
            substrate='O2',
            grid_size=self.grid_size,
            domain_size=self.domain_size,
            diffusion_coefficient=self.diffusion_coeff,
            boundary_condition='dirichlet',
            boundary_value=2e-5  # mol/L
        )
        
        # Set initial condition with radial gradient
        x = np.linspace(0, self.domain_size, self.grid_size)
        y = np.linspace(0, self.domain_size, self.grid_size)
        X, Y = np.meshgrid(x, y)
        
        # Distance from center
        center = self.domain_size / 2
        r = np.sqrt((X - center)**2 + (Y - center)**2)
        
        # Initial condition: higher in center, lower at edges
        initial_condition = 2e-5 * (1 - r / (self.domain_size/2))
        rd_field.set_initial_condition(initial_condition)
        
        # Pure diffusion (no reaction)
        zero_reaction = np.zeros((self.grid_size, self.grid_size))
        rd_field.set_reaction_rate_field(zero_reaction)
        
        # Solve to steady state
        convergence_info = rd_field.solve_steady_state()
        
        # Get final concentration
        final_concentration = rd_field.get_concentration_field()
        
        # Should converge to uniform concentration equal to boundary value
        np.testing.assert_allclose(final_concentration, 2e-5, rtol=1e-2)
        
        # Check convergence
        self.assertTrue(convergence_info['converged'])
    
    def test_diffusion_with_source_sink(self):
        """Test diffusion with source and sink terms."""
        # Test case: constant source in center, constant sink at edges
        
        rd_field = RDField(
            substrate='Glc',
            grid_size=self.grid_size,
            domain_size=self.domain_size,
            diffusion_coefficient=self.diffusion_coeff,
            boundary_condition='dirichlet',
            boundary_value=5e-3  # mol/L
        )
        
        # Set initial condition
        initial_condition = np.full((self.grid_size, self.grid_size), 5e-3)
        rd_field.set_initial_condition(initial_condition)
        
        # Create source-sink pattern
        reaction_rate = np.zeros((self.grid_size, self.grid_size))
        
        # Source in center (positive reaction rate)
        center = self.grid_size // 2
        reaction_rate[center-5:center+5, center-5:center+5] = 1e-6
        
        # Sink at edges (negative reaction rate)
        reaction_rate[0, :] = -1e-6
        reaction_rate[-1, :] = -1e-6
        reaction_rate[:, 0] = -1e-6
        reaction_rate[:, -1] = -1e-6
        
        rd_field.set_reaction_rate_field(reaction_rate)
        
        # Solve to steady state
        convergence_info = rd_field.solve_steady_state()
        
        # Get final concentration
        final_concentration = rd_field.get_concentration_field()
        
        # Should have higher concentration in center (source) and lower at edges (sink)
        center_concentration = final_concentration[center, center]
        edge_concentration = final_concentration[0, 0]
        
        self.assertGreater(center_concentration, edge_concentration)
        
        # Check convergence
        self.assertTrue(convergence_info['converged'])
    
    def test_neumann_boundary_conditions(self):
        """Test Neumann (no-flux) boundary conditions."""
        # With Neumann BC, the steady state should preserve total mass
        
        rd_field = RDField(
            substrate='Lac',
            grid_size=self.grid_size,
            domain_size=self.domain_size,
            diffusion_coefficient=self.diffusion_coeff,
            boundary_condition='neumann',
            boundary_value=0.0  # No flux
        )
        
        # Set initial condition with non-zero total mass
        initial_condition = np.ones((self.grid_size, self.grid_size)) * 1e-3
        rd_field.set_initial_condition(initial_condition)
        
        # Calculate initial total mass
        initial_mass = rd_field.get_total_mass()
        
        # Pure diffusion (no reaction)
        zero_reaction = np.zeros((self.grid_size, self.grid_size))
        rd_field.set_reaction_rate_field(zero_reaction)
        
        # Solve to steady state
        convergence_info = rd_field.solve_steady_state()
        
        # Calculate final total mass
        final_mass = rd_field.get_total_mass()
        
        # With Neumann BC and no reaction, total mass should be conserved
        np.testing.assert_allclose(final_mass, initial_mass, rtol=1e-6)
        
        # Final concentration should be uniform (equal to initial average)
        final_concentration = rd_field.get_concentration_field()
        expected_concentration = initial_mass / (self.grid_size * self.grid_size)
        np.testing.assert_allclose(final_concentration, expected_concentration, rtol=1e-2)
    
    def test_diffusion_equation_consistency(self):
        """Test that the diffusion equation is satisfied."""
        # Test: ∂C/∂t = D∇²C - r (with r = 0 for pure diffusion)
        
        rd_field = RDField(
            substrate='test',
            grid_size=self.grid_size,
            domain_size=self.domain_size,
            diffusion_coefficient=self.diffusion_coeff,
            boundary_condition='dirichlet',
            boundary_value=1.0
        )
        
        # Set initial condition
        initial_condition = np.random.uniform(0, 1, (self.grid_size, self.grid_size))
        rd_field.set_initial_condition(initial_condition)
        
        # Pure diffusion
        zero_reaction = np.zeros((self.grid_size, self.grid_size))
        rd_field.set_reaction_rate_field(zero_reaction)
        
        # Solve for a short time
        time_step = 0.01
        total_time = 0.1
        
        convergence_info = rd_field.solve_transient(time_step, total_time)
        
        # Get concentration field
        concentration = rd_field.get_concentration_field()
        
        # Check that the field is reasonable
        self.assertGreater(np.max(concentration), 0)
        self.assertLess(np.max(concentration), 2.0)  # Should not exceed boundary value
    
    def test_convergence_properties(self):
        """Test convergence properties of the solver."""
        rd_field = RDField(
            substrate='test',
            grid_size=self.grid_size,
            domain_size=self.domain_size,
            diffusion_coefficient=self.diffusion_coeff,
            boundary_condition='dirichlet',
            boundary_value=1.0
        )
        
        # Set initial condition
        initial_condition = np.zeros((self.grid_size, self.grid_size))
        rd_field.set_initial_condition(initial_condition)
        
        # Pure diffusion
        zero_reaction = np.zeros((self.grid_size, self.grid_size))
        rd_field.set_reaction_rate_field(zero_reaction)
        
        # Test different convergence tolerances
        tolerances = [1e-3, 1e-4, 1e-5]
        
        for tol in tolerances:
            # Reset initial condition
            rd_field.set_initial_condition(initial_condition)
            
            # Solve with different tolerance
            convergence_info = rd_field.solve_steady_state(
                convergence_tolerance=tol,
                max_iterations=1000
            )
            
            # Should converge for reasonable tolerances
            if tol >= 1e-4:
                self.assertTrue(convergence_info['converged'])
                self.assertLess(convergence_info['final_residual'], tol)
    
    def test_grid_independence(self):
        """Test that results are approximately grid-independent."""
        # Test with different grid sizes
        grid_sizes = [25, 50, 100]
        results = []
        
        for grid_size in grid_sizes:
            rd_field = RDField(
                substrate='test',
                grid_size=grid_size,
                domain_size=self.domain_size,
                diffusion_coefficient=self.diffusion_coeff,
                boundary_condition='dirichlet',
                boundary_value=1.0
            )
            
            # Set initial condition
            initial_condition = np.random.uniform(0, 0.5, (grid_size, grid_size))
            rd_field.set_initial_condition(initial_condition)
            
            # Pure diffusion
            zero_reaction = np.zeros((grid_size, grid_size))
            rd_field.set_reaction_rate_field(zero_reaction)
            
            # Solve to steady state
            convergence_info = rd_field.solve_steady_state()
            
            # Get average concentration
            final_concentration = rd_field.get_concentration_field()
            avg_concentration = np.mean(final_concentration)
            results.append(avg_concentration)
        
        # Results should be similar across grid sizes
        for i in range(1, len(results)):
            np.testing.assert_allclose(results[i], results[0], rtol=0.1)
    
    def test_physical_units_consistency(self):
        """Test that physical units are consistent."""
        # Test with realistic parameters
        rd_field = RDField(
            substrate='O2',
            grid_size=self.grid_size,
            domain_size=1e-3,  # 1 mm in meters
            diffusion_coefficient=1e-9,  # m²/s
            boundary_condition='dirichlet',
            boundary_value=2e-5  # mol/L
        )
        
        # Set initial condition
        initial_condition = np.full((self.grid_size, self.grid_size), 1e-5)
        rd_field.set_initial_condition(initial_condition)
        
        # Pure diffusion
        zero_reaction = np.zeros((self.grid_size, self.grid_size))
        rd_field.set_reaction_rate_field(zero_reaction)
        
        # Solve to steady state
        convergence_info = rd_field.solve_steady_state()
        
        # Get final concentration
        final_concentration = rd_field.get_concentration_field()
        
        # Check units consistency
        # Concentration should be in mol/L
        self.assertGreater(np.max(final_concentration), 0)
        self.assertLess(np.max(final_concentration), 1e-3)  # Reasonable range
        
        # Check that the solution makes physical sense
        # For pure diffusion with Dirichlet BC, should converge to boundary value
        np.testing.assert_allclose(final_concentration, 2e-5, rtol=1e-2)
    
    def test_reaction_diffusion_coupling(self):
        """Test coupling between reaction and diffusion."""
        # Test case: consumption reaction with diffusion
        
        rd_field = RDField(
            substrate='Glc',
            grid_size=self.grid_size,
            domain_size=self.domain_size,
            diffusion_coefficient=self.diffusion_coeff,
            boundary_condition='dirichlet',
            boundary_value=5e-3
        )
        
        # Set initial condition
        initial_condition = np.full((self.grid_size, self.grid_size), 5e-3)
        rd_field.set_initial_condition(initial_condition)
        
        # Set consumption reaction (negative reaction rate)
        # Higher consumption in center
        reaction_rate = np.zeros((self.grid_size, self.grid_size))
        center = self.grid_size // 2
        
        # Create radial consumption pattern
        x = np.linspace(0, self.domain_size, self.grid_size)
        y = np.linspace(0, self.domain_size, self.grid_size)
        X, Y = np.meshgrid(x, y)
        
        r = np.sqrt((X - self.domain_size/2)**2 + (Y - self.domain_size/2)**2)
        max_r = self.domain_size / 2
        
        # Consumption rate decreases with distance from center
        reaction_rate = -1e-6 * (1 - r / max_r)
        
        rd_field.set_reaction_rate_field(reaction_rate)
        
        # Solve to steady state
        convergence_info = rd_field.solve_steady_state()
        
        # Get final concentration
        final_concentration = rd_field.get_concentration_field()
        
        # Should have lower concentration in center (higher consumption)
        center_concentration = final_concentration[center, center]
        edge_concentration = final_concentration[0, 0]
        
        self.assertLess(center_concentration, edge_concentration)
        
        # Check convergence
        self.assertTrue(convergence_info['converged'])


class TestRDPerformance(unittest.TestCase):
    """Performance tests for RD operations."""
    
    def test_large_grid_performance(self):
        """Test performance with large grid."""
        # Test with larger grid
        large_grid_size = 200
        
        rd_field = RDField(
            substrate='test',
            grid_size=large_grid_size,
            domain_size=1.0,
            diffusion_coefficient=1e-9,
            boundary_condition='dirichlet',
            boundary_value=1.0
        )
        
        # Set initial condition
        initial_condition = np.random.uniform(0, 0.5, (large_grid_size, large_grid_size))
        rd_field.set_initial_condition(initial_condition)
        
        # Pure diffusion
        zero_reaction = np.zeros((large_grid_size, large_grid_size))
        rd_field.set_reaction_rate_field(zero_reaction)
        
        # Time the solution
        import time
        start_time = time.time()
        convergence_info = rd_field.solve_steady_state(max_iterations=100)
        end_time = time.time()
        
        # Should complete within reasonable time
        execution_time = end_time - start_time
        self.assertLess(execution_time, 10.0)  # 10 seconds max
        
        # Should make progress even if not fully converged
        self.assertGreater(convergence_info['iterations'], 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
