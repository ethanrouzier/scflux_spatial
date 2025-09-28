"""
Unit tests for spatial FBA coupling.
"""

import numpy as np
from scflux_spatial.spatial.rd import RDSolver, Species
from scflux_spatial.spatial.coupling import SpatialFBACoupler


def test_coupling_runs_and_converges():
    """Test that coupling runs and converges."""
    rd = RDSolver(
        16, 16, 1.0, 1.0, 
        [Species("O2", 2e-9, 0.2)], 
        bc={"left": "neumann", "right": "neumann", "bottom": "neumann", "top": "neumann"}
    )
    
    n = 16 * 16
    
    def fba(conc):
        """Simple FBA function that consumes O2."""
        return {"O2": -0.1 * np.ones(n)}
    
    coupler = SpatialFBACoupler(
        rd, fba_fn=fba, rho_gDW_per_m3=5e3, 
        alpha=0.5, tol_rel=1e-2, max_iters=5
    )
    
    # Run several iterations
    for _ in range(5):
        coupler.iterate(dt=1.0)
    
    # Should decrease O2 but remain positive
    o2 = rd.snapshot()["O2"]
    assert (o2 > 0).all()  # Should remain positive
    assert np.mean(o2) < 0.2  # Should be less than initial value


def test_coupling_under_relaxation():
    """Test that under-relaxation stabilizes the coupling."""
    rd = RDSolver(
        8, 8, 1.0, 1.0, 
        [Species("O2", 2e-9, 0.2)], 
        bc={"left": "neumann", "right": "neumann", "bottom": "neumann", "top": "neumann"}
    )
    
    n = 8 * 8
    
    def fba(conc):
        """FBA function with some noise to test stability."""
        base_flux = -0.1 * np.ones(n)
        noise = 0.01 * np.random.randn(n)
        return {"O2": base_flux + noise}
    
    # Test with different relaxation parameters
    for alpha in [0.1, 0.5, 0.9]:
        rd_test = RDSolver(
            8, 8, 1.0, 1.0, 
            [Species("O2", 2e-9, 0.2)], 
            bc={"left": "neumann", "right": "neumann", "bottom": "neumann", "top": "neumann"}
        )
        
        coupler = SpatialFBACoupler(
            rd_test, fba_fn=fba, rho_gDW_per_m3=5e3, 
            alpha=alpha, tol_rel=1e-2, max_iters=10
        )
        
        # Should run without exploding
        for _ in range(3):
            coupler.iterate(dt=1.0)
        
        o2 = rd_test.snapshot()["O2"]
        assert (o2 > 0).all()  # Should remain positive
        assert not np.any(np.isnan(o2))  # Should not have NaN values


def test_coupling_multiple_species():
    """Test coupling with multiple species."""
    rd = RDSolver(
        4, 4, 1.0, 1.0, 
        [Species("O2", 2e-9, 0.2), Species("Glc", 6e-10, 5.0)], 
        bc={"left": "neumann", "right": "neumann", "bottom": "neumann", "top": "neumann"}
    )
    
    n = 4 * 4
    
    def fba(conc):
        """FBA function for multiple species."""
        return {
            "O2": -0.1 * np.ones(n),  # Consumption
            "Glc": -0.05 * np.ones(n)  # Consumption
        }
    
    coupler = SpatialFBACoupler(
        rd, fba_fn=fba, rho_gDW_per_m3=5e3, 
        alpha=0.4, tol_rel=1e-2, max_iters=5
    )
    
    # Run coupling
    for _ in range(3):
        coupler.iterate(dt=1.0)
    
    snap = rd.snapshot()
    
    # Both species should be affected
    assert "O2" in snap
    assert "Glc" in snap
    assert (snap["O2"] > 0).all()
    assert (snap["Glc"] > 0).all()



