#!/usr/bin/env python3
"""
Example usage of CLI tools for scflux_spatial.

This script demonstrates how to use the command-line interfaces
for running flux balance analysis and spatial simulations.
"""

import subprocess
import sys
from pathlib import Path


def run_command(cmd, description):
    """Run a command and display the result."""
    print(f"\n🔧 {description}")
    print(f"Command: {' '.join(cmd)}")
    print("-" * 60)
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        
        if result.stdout:
            print("STDOUT:")
            print(result.stdout)
        
        if result.stderr:
            print("STDERR:")
            print(result.stderr)
        
        print(f"Return code: {result.returncode}")
        
        if result.returncode == 0:
            print("✅ Command executed successfully")
        else:
            print("❌ Command failed")
            
    except subprocess.TimeoutExpired:
        print("⏰ Command timed out")
    except Exception as e:
        print(f"❌ Error running command: {e}")


def main():
    """Main function demonstrating CLI usage."""
    print("🧬 scflux_spatial CLI Usage Examples")
    print("=" * 60)
    
    # Example 1: Run FBA with E-Flux method
    print("\n📊 Example 1: FBA with E-Flux method")
    run_command([
        "python", "-m", "scflux_spatial.cli.run_flux",
        "--demo",
        "--method", "eflux",
        "--objective", "ATP",
        "--output", "fba_eflux_results"
    ], "Run FBA analysis with E-Flux integration method")
    
    # Example 2: Run FBA with iMAT method
    print("\n📊 Example 2: FBA with iMAT method")
    run_command([
        "python", "-m", "scflux_spatial.cli.run_flux",
        "--demo",
        "--method", "imat",
        "--objective", "biomass",
        "--pfba",
        "--output", "fba_imat_results"
    ], "Run FBA analysis with iMAT integration method and pFBA")
    
    # Example 3: Run spatial simulation with O2 and glucose
    print("\n🌊 Example 3: Spatial simulation with O2 and glucose")
    run_command([
        "python", "-m", "scflux_spatial.cli.run_spatial",
        "--substrates", "O2,Glc",
        "--iters", "5",
        "--tol", "1e-3",
        "--grid", "50",
        "--demo",
        "--output", "spatial_o2_glc_results"
    ], "Run spatial RD simulation for O2 and glucose")
    
    # Example 4: Run spatial simulation with multiple substrates
    print("\n🌊 Example 4: Spatial simulation with multiple substrates")
    run_command([
        "python", "-m", "scflux_spatial.cli.run_spatial",
        "--substrates", "O2,Glc,Lac",
        "--iters", "10",
        "--tol", "1e-4",
        "--grid", "100",
        "--domain", "2.0",
        "--D_O2", "1e-9",
        "--D_Glc", "1e-9",
        "--D_Lac", "1e-9",
        "--demo",
        "--output", "spatial_multi_substrate_results"
    ], "Run spatial RD simulation for multiple substrates")
    
    # Example 5: Help for FBA CLI
    print("\n❓ Example 5: FBA CLI help")
    run_command([
        "python", "-m", "scflux_spatial.cli.run_flux",
        "--help"
    ], "Display help for FBA CLI")
    
    # Example 6: Help for spatial CLI
    print("\n❓ Example 6: Spatial CLI help")
    run_command([
        "python", "-m", "scflux_spatial.cli.run_spatial",
        "--help"
    ], "Display help for spatial CLI")
    
    print(f"\n🎉 CLI usage examples completed!")
    print(f"\nKey CLI commands:")
    print(f"  FBA Analysis:")
    print(f"    python -m scflux_spatial.cli.run_flux --demo --method eflux --objective ATP")
    print(f"    python -m scflux_spatial.cli.run_flux --method imat --objective biomass --pfba")
    print(f"  Spatial Simulation:")
    print(f"    python -m scflux_spatial.cli.run_spatial --substrates O2,Glc --iters 10 --tol 1e-4")
    print(f"    python -m scflux_spatial.cli.run_spatial --substrates O2,Glc,Lac --grid 100 --demo")


if __name__ == "__main__":
    main()
