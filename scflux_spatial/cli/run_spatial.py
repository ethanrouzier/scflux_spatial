"""
Smoke test CLI for spatial flux analysis.

This module provides a command-line interface to test the basic
functionality of the spatial flux analysis system.
"""

import argparse
import numpy as np
from scflux_spatial.spatial.rd import RDSolver, Species
from scflux_spatial.spatial.coupling import SpatialFBACoupler


def dummy_fba(conc):  # conc: Dict[str, np.ndarray] mol/m^3
    """
    Dummy FBA function for testing.
    
    Args:
        conc: Dictionary of concentrations in mol/m^3
        
    Returns:
        Dictionary of fluxes in mmol/gDW/h
    """
    # Uniform O2 consumption (mmol/gDW/h), zero production for Glc
    n = next(iter(conc.values())).size
    return {"O2": -0.5*np.ones(n), "Glc": 0.0*np.ones(n)}  # simple demo


def main():
    """Main CLI function."""
    ap = argparse.ArgumentParser(description="Spatial flux analysis smoke test")
    ap.add_argument("--grid", type=int, default=64, help="Grid size (nx=ny)")
    ap.add_argument("--iters", type=int, default=5, help="Number of time steps")
    ap.add_argument("--domain", type=float, default=1.0, help="Domain size in meters")
    args = ap.parse_args()

    print(f"Running spatial flux analysis smoke test...")
    print(f"  Grid: {args.grid}x{args.grid}")
    print(f"  Domain: {args.domain}m x {args.domain}m")
    print(f"  Time steps: {args.iters}")

    # Create RD solver
    rd = RDSolver(
        nx=args.grid, ny=args.grid, Lx=args.domain, Ly=args.domain,
        species=[Species("O2", 2e-9, 0.2), Species("Glc", 6e-10, 5.0)],
        bc={"left": "neumann", "right": "neumann", "bottom": "neumann", "top": "neumann"}
    )
    
    # Create coupler
    coupler = SpatialFBACoupler(
        rd, fba_fn=dummy_fba, rho_gDW_per_m3=5e3, 
        alpha=0.4, tol_rel=1e-3, max_iters=20
    )

    # Calculate safe time step (not critical for implicit scheme)
    dt = (args.domain/args.grid)**2 / (4*2e-9)
    print(f"  Time step: {dt:.2e} s")

    # Run simulation
    for i in range(args.iters):
        coupler.iterate(dt=dt)
        if i % max(1, args.iters//5) == 0:
            print(f"  Step {i+1}/{args.iters}")

    # Get final results
    snap = rd.snapshot()
    print(f"\nResults:")
    print(f"  O2 mean: {float(np.mean(snap['O2'])):.3e} mol/m^3")
    print(f"  Glc mean: {float(np.mean(snap['Glc'])):.3e} mol/m^3")
    print(f"  O2 min/max: {float(np.min(snap['O2'])):.3e} / {float(np.max(snap['O2'])):.3e} mol/m^3")
    print(f"  Glc min/max: {float(np.min(snap['Glc'])):.3e} / {float(np.max(snap['Glc'])):.3e} mol/m^3")
    
    print("\n✅ Smoke test completed successfully!")


if __name__ == "__main__":
    main()