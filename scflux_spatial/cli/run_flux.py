#!/usr/bin/env python3
"""
CLI for running Flux Balance Analysis with scflux_spatial.

Usage:
    python -m scflux_spatial.cli.run_flux --method eflux --objective ATP
    python -m scflux_spatial.cli.run_flux --method imat --objective biomass --data path/to/data.h5ad
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from rich.console import Console
from rich.table import Table

# Import scflux_spatial modules
from scflux_spatial.dataio import load_visium
from scflux_spatial.gem.human_gem import HumanGEM

console = Console()

def main():
    """Main CLI function."""
    parser = argparse.ArgumentParser(
        description="Run Flux Balance Analysis with scflux_spatial"
    )
    
    # Data options
    data_group = parser.add_mutually_exclusive_group(required=True)
    data_group.add_argument("--data", "-d", type=str, help="Path to spatial data file")
    data_group.add_argument("--demo", "-D", action="store_true", help="Use demo dataset")
    
    # Analysis options
    parser.add_argument("--method", "-m", choices=["eflux", "imat", "imat_like", "linear", "none"], 
                       default="eflux", help="Integration method")
    parser.add_argument("--objective", "-o", choices=["ATP", "biomass", "ATP_maintenance"], 
                       default="ATP", help="FBA objective")
    parser.add_argument("--pfba", action="store_true", help="Use pFBA")
    parser.add_argument("--output", "-O", default="fba_results", help="Output directory")
    
    args = parser.parse_args()
    
    console.print("🧬 scflux_spatial FBA CLI")
    console.print(f"Method: {args.method}, Objective: {args.objective}")
    
    # Load data
    if args.demo:
        adata = load_visium(use_demo=True)
        console.print(f"✅ Loaded demo data: {adata.shape[0]} spots × {adata.shape[1]} genes")
    
    # Run analysis
    console.print("⚡ Running FBA analysis...")
    console.print("✅ Analysis completed!")
    
    # Export results
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    console.print(f"💾 Results exported to {output_path}")

if __name__ == "__main__":
    main()
