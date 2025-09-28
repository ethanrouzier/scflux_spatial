"""
Unit tests for reaction-diffusion solver.
"""

import numpy as np
from scflux_spatial.spatial.rd import RDSolver, Species


def test_rd_uniform_bc_no_source_stable():
    """Test that RD solver is stable with uniform BC and no sources."""
    rd = RDSolver(
        nx=16, ny=16, Lx=1.0, Ly=1.0,
        species=[Species("O2", 2e-9, 0.2)], 
        bc={"left": "neumann", "right": "neumann", "bottom": "neumann", "top": "neumann"}
    )
    base = rd.snapshot()["O2"].copy()
    
    # Run several steps with no source
    for _ in range(3):
        rd.set_source("O2", 0.0 * base)
        rd.step(dt=1.0)
    
    # Should remain stable (no change)
    diff = np.linalg.norm(rd.snapshot()["O2"] - base)
    assert diff < 1e-9


def test_rd_source_consumption():
    """Test that sources are properly applied."""
    rd = RDSolver(
        nx=8, ny=8, Lx=1.0, Ly=1.0,
        species=[Species("O2", 2e-9, 0.2)], 
        bc={"left": "neumann", "right": "neumann", "bottom": "neumann", "top": "neumann"}
    )
    
    initial = rd.snapshot()["O2"].copy()
    
    # Apply uniform consumption
    n_cells = 8 * 8
    consumption = -0.1 * np.ones(n_cells)  # mol/m^3/s
    rd.set_source("O2", consumption)
    rd.step(dt=1.0)
    
    final = rd.snapshot()["O2"]
    
    # Should decrease due to consumption
    assert np.mean(final) < np.mean(initial)


def test_rd_multiple_species():
    """Test RD solver with multiple species."""
    rd = RDSolver(
        nx=4, ny=4, Lx=1.0, Ly=1.0,
        species=[
            Species("O2", 2e-9, 0.2),
            Species("Glc", 6e-10, 5.0)
        ], 
        bc={"left": "neumann", "right": "neumann", "bottom": "neumann", "top": "neumann"}
    )
    
    # Test that both species are tracked
    snap = rd.snapshot()
    assert "O2" in snap
    assert "Glc" in snap
    assert snap["O2"].shape == (16,)  # 4x4 grid
    assert snap["Glc"].shape == (16,)
    
    # Test source setting for both species
    n_cells = 16
    rd.set_source("O2", -0.1 * np.ones(n_cells))
    rd.set_source("Glc", 0.05 * np.ones(n_cells))
    rd.step(dt=1.0)
    
    # Should have changed
    new_snap = rd.snapshot()
    assert not np.allclose(snap["O2"], new_snap["O2"])
    assert not np.allclose(snap["Glc"], new_snap["Glc"])



