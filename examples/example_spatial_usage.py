#!/usr/bin/env python3
"""
Example usage of spatial modeling functionality.

This script demonstrates the new spatial features:
- SpatialKineticsModel for rate conversion
- RDField for reaction-diffusion PDEs
- SpatialFBACoupler for SOA loop
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Import scflux_spatial modules
from scflux_spatial.spatial.kinetics import SpatialKineticsModel, create_default_michaelis_menten_params
from scflux_spatial.spatial.rd import RDField
from scflux_spatial.spatial.coupling import SpatialFBACoupler


def main():
    """Main example function."""
    print("🌍 scflux_spatial - Spatial Modeling Example")
    print("=" * 50)
    
    # 1. Initialize spatial kinetics model
    print("\n1. Initializing spatial kinetics model...")
    kinetics_model = SpatialKineticsModel(
        cell_density=1e9,  # cells·L⁻¹
        spot_density=1e6,  # spots·L⁻¹
        cell_volume=1e-12,  # L per cell
        spot_volume=1e-9,   # L per spot
    )
    
    # 2. Test FBA to spatial rate conversion
    print("\n2. Testing FBA to spatial rate conversion...")
    
    # Mock FBA uptake rate (mmol·gDW⁻¹·h⁻¹)
    fba_glucose_uptake = -10.0  # Glucose uptake
    
    # Convert to spatial consumption rate
    spatial_rate = kinetics_model.convert_fba_to_spatial_rate(
        fba_rate=fba_glucose_uptake,
        spatial_type="spot",
        use_michaelis_menten=False,
    )
    
    print(f"   FBA glucose uptake: {fba_glucose_uptake} mmol·gDW⁻¹·h⁻¹")
    print(f"   Spatial consumption rate: {spatial_rate:.2e} mol·L⁻¹·s⁻¹")
    
    # Test with Michaelis-Menten kinetics
    mm_params = create_default_michaelis_menten_params()
    glucose_mm = mm_params["glucose"]
    
    spatial_rate_mm = kinetics_model.convert_fba_to_spatial_rate(
        fba_rate=fba_glucose_uptake,
        spatial_type="spot",
        use_michaelis_menten=True,
        concentration=1e-3,  # 1 mM glucose
        vmax=glucose_mm["vmax"],
        km=glucose_mm["km"],
    )
    
    print(f"   With MM kinetics: {spatial_rate_mm:.2e} mol·L⁻¹·s⁻¹")
    
    # 3. Create reaction-diffusion fields
    print("\n3. Creating reaction-diffusion fields...")
    
    try:
        # Create RDField for glucose
        glucose_field = RDField(
            substrate="glucose",
            grid_size=(50, 50),
            domain_size=(1e-3, 1e-3),  # 1mm x 1mm
            diffusion_coefficient=1e-9,  # m²·s⁻¹
            boundary_condition="dirichlet",
            boundary_value=5e-3,  # 5 mM at boundaries
        )
        
        # Set initial condition (uniform concentration)
        glucose_field.set_initial_condition(1e-3)  # 1 mM initial
        
        # Create reaction rate field
        reaction_rates = np.zeros((50, 50))
        # Add some spatial heterogeneity
        for i in range(50):
            for j in range(50):
                # Higher consumption in center
                center_dist = np.sqrt((i-25)**2 + (j-25)**2)
                if center_dist < 10:
                    reaction_rates[i, j] = spatial_rate * 2
                else:
                    reaction_rates[i, j] = spatial_rate * 0.5
        
        glucose_field.set_reaction_rate_field(reaction_rates)
        
        print("   ✓ Glucose RDField created")
        
    except ImportError as e:
        print(f"   ⚠️ FiPy not available: {e}")
        print("   Skipping RDField example...")
        glucose_field = None
    
    # 4. Test steady state solution
    if glucose_field is not None:
        print("\n4. Solving steady state...")
        
        try:
            convergence_info = glucose_field.solve_steady_state(
                tolerance=1e-6,
                max_iterations=100,
                relaxation_factor=0.9
            )
            
            print(f"   ✓ Steady state solved")
            print(f"   Converged: {convergence_info['converged']}")
            print(f"   Iterations: {convergence_info['iterations']}")
            print(f"   Final residual: {convergence_info['final_residual']:.2e}")
            
            # Get final concentration field
            conc_field = glucose_field.get_concentration_field()
            print(f"   Concentration range: {conc_field.min():.2e} - {conc_field.max():.2e} mol·L⁻¹")
            
        except Exception as e:
            print(f"   ⚠️ Steady state solve failed: {e}")
    
    # 5. Demonstrate spatial rate matrix calculation
    print("\n5. Calculating spatial rate matrix...")
    
    # Mock spot coordinates and fluxes
    spot_coordinates = {
        "spot_001": (0.5e-3, 0.5e-3),  # Center
        "spot_002": (0.3e-3, 0.7e-3),  # Corner
        "spot_003": (0.7e-3, 0.3e-3),  # Another corner
    }
    
    substrate_fluxes = {
        "spot_001": {"glucose": -8.0, "oxygen": -50.0, "lactate": 15.0},
        "spot_002": {"glucose": -5.0, "oxygen": -30.0, "lactate": 10.0},
        "spot_003": {"glucose": -12.0, "oxygen": -60.0, "lactate": 20.0},
    }
    
    substrate_concentrations = {
        "glucose": 1e-3,  # 1 mM
        "oxygen": 2e-5,   # 20 μM
        "lactate": 2e-4,  # 0.2 mM
    }
    
    spatial_rates = kinetics_model.get_spatial_rate_matrix(
        spot_coordinates=spot_coordinates,
        substrate_fluxes=substrate_fluxes,
        substrate_concentrations=substrate_concentrations,
        spatial_type="spot",
        use_michaelis_menten=True,
        mm_parameters=mm_params,
    )
    
    print("   Spatial rates calculated:")
    for spot_id, substrates in spatial_rates.items():
        print(f"   {spot_id}:")
        for substrate, rate in substrates.items():
            print(f"     {substrate}: {rate:.2e} mol·L⁻¹·s⁻¹")
    
    # 6. Export results
    print("\n6. Exporting results...")
    
    output_dir = Path("spatial_example_output")
    output_dir.mkdir(exist_ok=True)
    
    # Export spatial rates
    kinetics_model.export_spatial_rates(
        spatial_rates, 
        output_dir / "spatial_rates.csv"
    )
    
    # Export field data if available
    if glucose_field is not None:
        glucose_field.export_field_data(output_dir / "glucose_field")
    
    print(f"   ✓ Results exported to {output_dir}")
    
    # 7. Summary statistics
    print("\n7. Summary statistics...")
    
    print(f"   Cell density: {kinetics_model.cell_density:.2e} cells·L⁻¹")
    print(f"   Spot density: {kinetics_model.spot_density:.2e} spots·L⁻¹")
    print(f"   Conversion factors:")
    print(f"     mmol/mol: {kinetics_model.mmol_per_mol}")
    print(f"     s/h: {kinetics_model.h_per_s}")
    print(f"     gDW/cell: {kinetics_model.gDW_per_cell:.2e}")
    
    if glucose_field is not None:
        total_mass = glucose_field.get_total_mass()
        print(f"   Total glucose mass: {total_mass:.2e} mol")
    
    print(f"\n✅ Spatial modeling example completed successfully!")
    print(f"\nTo run SOA coupling:")
    print(f"   1. Create RDField instances for each substrate")
    print(f"   2. Create FBA solvers for each spot")
    print(f"   3. Initialize SpatialFBACoupler")
    print(f"   4. Run coupler.run_soa_loop()")


if __name__ == "__main__":
    main()
