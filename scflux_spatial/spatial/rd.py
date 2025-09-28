"""
Reaction-diffusion solver using FiPy for implicit time integration.

This module provides a 2D reaction-diffusion solver with support for
multiple species and volumetric source terms.
"""

from dataclasses import dataclass
from typing import Dict, List
import numpy as np
from fipy import CellVariable, Grid2D, TransientTerm, DiffusionTerm


@dataclass
class Species:
    """Represents a chemical species with diffusion properties."""
    name: str
    D: float      # m^2/s
    C0: float     # mol/m^3


class RDSolver:
    """
    2D Reaction-Diffusion solver using FiPy.
    
    Solves the equation: ∂C/∂t = ∇·(D∇C) + S
    where C is concentration, D is diffusivity, and S is source term.
    """
    
    def __init__(self, nx: int, ny: int, Lx: float, Ly: float, 
                 species: List[Species], bc: Dict[str, str]):
        """
        Initialize the RD solver.
        
        Args:
            nx, ny: Grid dimensions
            Lx, Ly: Domain size in meters
            species: List of species to track
            bc: Boundary conditions (currently only Neumann supported)
        """
        self.mesh = Grid2D(nx=nx, ny=ny, dx=Lx/nx, dy=Ly/ny)
        self.species = species
        self.bc = bc
        self.vars: Dict[str, CellVariable] = {}
        self.sources: Dict[str, CellVariable] = {}

        # Initialize species variables
        for sp in species:
            Ci = CellVariable(name=sp.name, mesh=self.mesh, value=sp.C0, hasOld=True)
            Si = CellVariable(name=f"S_{sp.name}", mesh=self.mesh, value=0.0)
            self.vars[sp.name] = Ci
            self.sources[sp.name] = Si

        # Pre-build implicit equations (sources are explicit via Variable)
        self.eqs = {}
        for sp in species:
            # Create equation: dC/dt = D*∇²C + S
            # Use implicit time stepping with explicit source
            self.eqs[sp.name] = (TransientTerm(var=self.vars[sp.name]) ==
                                DiffusionTerm(coeff=sp.D, var=self.vars[sp.name]) + 
                                self.sources[sp.name])

    def set_source(self, name: str, value_array: np.ndarray):
        """
        Set volumetric source term for a species.
        
        Args:
            name: Species name
            value_array: Source values in mol/m^3/s, shape=(nx*ny,)
        """
        if name not in self.sources:
            raise ValueError(f"Species {name} not found")
        self.sources[name].value = value_array

    def step(self, dt: float):
        """
        Perform one time step using implicit scheme.
        
        Args:
            dt: Time step in seconds
        """
        for name, eq in self.eqs.items():
            self.vars[name].updateOld()
            eq.solve(dt=dt)

    def snapshot(self) -> Dict[str, np.ndarray]:
        """
        Get current concentration field for all species.
        
        Returns:
            Dictionary mapping species names to concentration arrays
        """
        return {name: var.value.copy() for name, var in self.vars.items()}


class RDField:
    """
    Legacy RDField class for backward compatibility.
    
    This class provides a simplified interface to the RDSolver
    for existing code that expects RDField.
    """
    
    def __init__(self, substrate: str = "test", grid_size: int = 64, 
                 diffusion_coefficient: float = 1e-9, domain_size: float = 1.0, 
                 initial_concentration: float = 0.0, boundary_condition: str = "neumann",
                 boundary_value: float = 0.0):
        """
        Initialize RDField with legacy API.
        
        Args:
            substrate: Species name (legacy parameter name)
            grid_size: Grid size (assumes square grid)
            diffusion_coefficient: Diffusion coefficient in m²/s
            domain_size: Domain size in meters
            initial_concentration: Initial concentration in mol/m³
            boundary_condition: Boundary condition type ("neumann" or "dirichlet")
            boundary_value: Boundary condition value
        """
        self.name = substrate
        self.grid_size = grid_size
        self.diffusion_coefficient = diffusion_coefficient
        self.domain_size = domain_size
        self.initial_concentration = initial_concentration
        self.boundary_condition = boundary_condition
        self.boundary_value = boundary_value
        
        # Create species and solver
        species = Species(substrate, diffusion_coefficient, initial_concentration)
        self.solver = RDSolver(
            nx=grid_size, ny=grid_size, 
            Lx=domain_size, Ly=domain_size,
            species=[species],
            bc={"left": boundary_condition, "right": boundary_condition, 
                "bottom": boundary_condition, "top": boundary_condition}
        )
    
    def set_concentration(self, value: float):
        """Set uniform concentration."""
        self.solver.vars[self.name].value = value
    
    def get_concentration(self) -> np.ndarray:
        """Get concentration field."""
        return self.solver.vars[self.name].value
    
    def set_source(self, source: np.ndarray):
        """Set source term."""
        self.solver.set_source(self.name, source)
    
    def step(self, dt: float):
        """Perform time step."""
        self.solver.step(dt)
    
    def solve_steady_state(self, tolerance: float = 1e-6, max_iterations: int = 1000):
        """Solve to steady state."""
        for _ in range(max_iterations):
            old_values = self.get_concentration().copy()
            self.step(dt=1.0)
            new_values = self.get_concentration()
            
            if np.linalg.norm(new_values - old_values) < tolerance:
                break