#!/usr/bin/env python3
"""
Example usage of visualization functionality.

This script demonstrates the new visualization features:
- heatmap_concentration with spot overlays
- field_quiver for gradient visualization
- spot_metric_map for FBA results
- Escher cell type flux maps
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Import scflux_spatial modules
from scflux_spatial.viz.maps import SpatialMapper
from scflux_spatial.viz.escher_view import EscherViewer
from scflux_spatial.dataio import load_visium


def main():
    """Main example function."""
    print("🎨 scflux_spatial - Visualization Example")
    print("=" * 50)
    
    # 1. Load demo data
    print("\n1. Loading demo Visium data...")
    adata = load_visium(use_demo=True)
    print(f"   ✓ Loaded: {adata.shape[0]} spots × {adata.shape[1]} genes")
    
    # 2. Initialize spatial mapper
    print("\n2. Initializing spatial mapper...")
    mapper = SpatialMapper(
        colormap="viridis",
        figure_size=(10, 8),
        dpi=300
    )
    
    # 3. Create mock concentration fields
    print("\n3. Creating mock concentration fields...")
    np.random.seed(42)
    
    # Oxygen concentration field (2D)
    O2_field = np.random.uniform(1e-5, 2e-5, (100, 100))
    O2_field[0, :] = O2_field[-1, :] = O2_field[:, 0] = O2_field[:, -1] = 2e-5  # Boundaries
    
    # Glucose concentration field (2D)
    Glc_field = np.random.uniform(1e-3, 5e-3, (100, 100))
    Glc_field[0, :] = Glc_field[-1, :] = Glc_field[:, 0] = Glc_field[:, -1] = 5e-3  # Boundaries
    
    print(f"   ✓ Created O2 field: {O2_field.shape}, range: {O2_field.min():.2e} - {O2_field.max():.2e} mol/L")
    print(f"   ✓ Created Glc field: {Glc_field.shape}, range: {Glc_field.min():.2e} - {Glc_field.max():.2e} mol/L")
    
    # 4. Create concentration heatmaps with spot overlays
    print("\n4. Creating concentration heatmaps with spot overlays...")
    
    # Oxygen heatmap
    fig_o2 = mapper.heatmap_concentration(
        adata=adata,
        C_field=O2_field,
        field_name="oxygen",
        overlay_spots=True,
        spot_size=6,
        spot_opacity=0.8,
        title="Oxygen Concentration with Visium Spots"
    )
    
    # Glucose heatmap
    fig_glc = mapper.heatmap_concentration(
        adata=adata,
        C_field=Glc_field,
        field_name="glucose",
        overlay_spots=True,
        spot_size=6,
        spot_opacity=0.8,
        title="Glucose Concentration with Visium Spots"
    )
    
    print("   ✓ Created concentration heatmaps with spot overlays")
    
    # 5. Create gradient quiver plots
    print("\n5. Creating gradient quiver plots...")
    
    # Oxygen gradients
    fig_o2_quiver = mapper.field_quiver(
        C_field=O2_field,
        field_name="oxygen",
        skip=8,
        scale=0.5,
        title="Oxygen Concentration Gradients"
    )
    
    # Glucose gradients
    fig_glc_quiver = mapper.field_quiver(
        C_field=Glc_field,
        field_name="glucose",
        skip=8,
        scale=0.5,
        title="Glucose Concentration Gradients"
    )
    
    print("   ✓ Created gradient quiver plots")
    
    # 6. Create spot metric maps
    print("\n6. Creating spot metric maps...")
    
    # Mock ATP flux data
    np.random.seed(123)
    atp_flux = np.random.normal(10, 2, adata.n_obs)
    
    # Mock lactate export data
    lactate_export = np.random.normal(5, 1.5, adata.n_obs)
    
    # ATP flux map
    fig_atp = mapper.spot_metric_map(
        adata=adata,
        metric_name="ATP Flux",
        metric_values=atp_flux,
        colormap="RdYlBu_r",
        title="ATP Production Flux Distribution",
        show_clusters=True
    )
    
    # Lactate export map
    fig_lactate = mapper.spot_metric_map(
        adata=adata,
        metric_name="Lactate Export",
        metric_values=lactate_export,
        colormap="Reds",
        title="Lactate Export Distribution",
        show_clusters=False
    )
    
    print("   ✓ Created spot metric maps")
    
    # 7. Create multi-field comparison
    print("\n7. Creating multi-field comparison...")
    
    field_dict = {
        "Oxygen": O2_field,
        "Glucose": Glc_field,
    }
    
    fig_comparison = mapper.create_multi_field_comparison(
        adata=adata,
        field_dict=field_dict,
        overlay_spots=True,
        title="Multi-Field Concentration Comparison"
    )
    
    print("   ✓ Created multi-field comparison")
    
    # 8. Escher visualization
    print("\n8. Creating Escher cell type flux maps...")
    
    try:
        # Create mock cell type fluxes
        cell_types = adata.obs['leiden'].unique() if 'leiden' in adata.obs else ['Cluster_0', 'Cluster_1', 'Cluster_2']
        
        cell_type_fluxes = {}
        for i, ct in enumerate(cell_types):
            # Mock central carbon metabolism fluxes
            cell_type_fluxes[ct] = {
                'HEX1': np.random.normal(2.5, 0.5),    # Hexokinase
                'PFK': np.random.normal(1.8, 0.3),     # Phosphofructokinase
                'GAPDH': np.random.normal(3.2, 0.7),   # Glyceraldehyde-3-phosphate dehydrogenase
                'LDHA': np.random.normal(2.8, 0.6),    # Lactate dehydrogenase A
                'CS': np.random.normal(1.5, 0.4),      # Citrate synthase
                'MDH': np.random.normal(2.0, 0.5),     # Malate dehydrogenase
                'PDH': np.random.normal(1.2, 0.3),     # Pyruvate dehydrogenase
                'ACONT': np.random.normal(1.0, 0.2),   # Aconitase
                'ICDHyr': np.random.normal(0.8, 0.2),  # Isocitrate dehydrogenase
                'AKGDH': np.random.normal(0.9, 0.2),   # α-ketoglutarate dehydrogenase
            }
        
        # Create Escher viewer
        escher_viewer = EscherViewer()
        
        # Create cell type flux maps
        cell_type_maps = escher_viewer.create_cell_type_flux_map(
            cell_type_fluxes=cell_type_fluxes,
            map_type="central_carbon",
            title_prefix="Cluster",
            width=800,
            height=600
        )
        
        print(f"   ✓ Created Escher maps for {len(cell_type_maps)} cell types")
        
        # Create flux comparison map
        comparison_map = escher_viewer.create_flux_comparison_map(
            cell_type_fluxes=cell_type_fluxes,
            map_type="central_carbon",
            comparison_method="difference",
            reference_type=cell_types[0],
            title="Flux Comparison"
        )
        
        print("   ✓ Created flux comparison map")
        
        # Export flux statistics
        escher_viewer.create_flux_statistics_summary(
            cell_type_fluxes=cell_type_fluxes,
            output_path="escher_flux_summary"
        )
        
        print("   ✓ Exported flux statistics")
        
    except ImportError:
        print("   ⚠️ Escher not available, skipping Escher visualizations")
    except Exception as e:
        print(f"   ⚠️ Error creating Escher maps: {e}")
    
    # 9. Export visualizations
    print("\n9. Exporting visualizations...")
    
    output_dir = Path("visualization_example_output")
    output_dir.mkdir(exist_ok=True)
    
    # Export concentration heatmaps
    mapper.save_figure(fig_o2, output_dir / "oxygen_concentration.html", format="html")
    mapper.save_figure(fig_glc, output_dir / "glucose_concentration.html", format="html")
    
    # Export quiver plots
    mapper.save_figure(fig_o2_quiver, output_dir / "oxygen_gradients.html", format="html")
    mapper.save_figure(fig_glc_quiver, output_dir / "glucose_gradients.html", format="html")
    
    # Export metric maps
    mapper.save_figure(fig_atp, output_dir / "atp_flux_map.html", format="html")
    mapper.save_figure(fig_lactate, output_dir / "lactate_export_map.html", format="html")
    
    # Export comparison
    mapper.save_figure(fig_comparison, output_dir / "multi_field_comparison.html", format="html")
    
    print(f"   ✓ Visualizations exported to {output_dir}")
    
    # 10. Summary statistics
    print("\n10. Summary statistics...")
    
    print(f"   Data dimensions: {adata.shape[0]} spots × {adata.shape[1]} genes")
    print(f"   Concentration fields: {len(field_dict)} fields")
    print(f"   O2 concentration range: {O2_field.min():.2e} - {O2_field.max():.2e} mol/L")
    print(f"   Glucose concentration range: {Glc_field.min():.2e} - {Glc_field.max():.2e} mol/L")
    print(f"   ATP flux range: {atp_flux.min():.2f} - {atp_flux.max():.2f}")
    print(f"   Lactate export range: {lactate_export.min():.2f} - {lactate_export.max():.2f}")
    
    if 'cell_type_fluxes' in locals():
        print(f"   Cell types analyzed: {len(cell_type_fluxes)}")
        print(f"   Reactions per cell type: {len(list(cell_type_fluxes.values())[0])}")
    
    print(f"\n✅ Visualization example completed successfully!")
    print(f"\nGenerated files:")
    print(f"   - Concentration heatmaps with spot overlays")
    print(f"   - Gradient quiver plots")
    print(f"   - Spot metric maps (ATP flux, lactate export)")
    print(f"   - Multi-field comparison")
    print(f"   - Escher cell type flux maps (if available)")
    print(f"   - Flux statistics summary")


if __name__ == "__main__":
    main()
