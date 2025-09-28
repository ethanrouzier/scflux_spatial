#!/usr/bin/env python3
"""
Example usage of the scflux_spatial load_visium function.

This script demonstrates how to use the new load_visium function
to load and analyze Visium data with metabolic pathway scoring.
"""

from scflux_spatial.dataio import load_visium
import matplotlib.pyplot as plt
import numpy as np


def main():
    """Main example function."""
    print("🧬 scflux_spatial - Visium Data Loading Example")
    print("=" * 50)
    
    # Load demo dataset
    print("\n1. Loading demo mouse brain Visium dataset...")
    adata = load_visium(use_demo=True)
    
    print(f"   ✓ Dataset loaded successfully!")
    print(f"   ✓ Shape: {adata.shape[0]} spots × {adata.shape[1]} genes")
    print(f"   ✓ Spatial coordinates: {adata.obsm['spatial_coordinates'].shape}")
    
    # Display available columns
    print(f"\n2. Available observation columns:")
    for col in adata.obs.columns:
        print(f"   - {col}")
    
    # Display metabolic scores
    print(f"\n3. Metabolic pathway scores:")
    print(f"   - Glycolysis score: {adata.obs['glycolysis_score'].mean():.3f} ± {adata.obs['glycolysis_score'].std():.3f}")
    print(f"   - OXPHOS score: {adata.obs['oxphos_score'].mean():.3f} ± {adata.obs['oxphos_score'].std():.3f}")
    print(f"   - Metabolic activity: {adata.obs['metabolic_activity'].mean():.3f} ± {adata.obs['metabolic_activity'].std():.3f}")
    
    # Create a simple visualization
    print(f"\n4. Creating visualization...")
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Plot glycolysis scores
    scatter1 = axes[0].scatter(adata.obs['x'], adata.obs['y'], 
                              c=adata.obs['glycolysis_score'], 
                              cmap='viridis', s=10, alpha=0.7)
    axes[0].set_title('Glycolysis Score')
    axes[0].set_xlabel('X Position')
    axes[0].set_ylabel('Y Position')
    plt.colorbar(scatter1, ax=axes[0])
    
    # Plot OXPHOS scores
    scatter2 = axes[1].scatter(adata.obs['x'], adata.obs['y'], 
                              c=adata.obs['oxphos_score'], 
                              cmap='plasma', s=10, alpha=0.7)
    axes[1].set_title('OXPHOS Score')
    axes[1].set_xlabel('X Position')
    axes[1].set_ylabel('Y Position')
    plt.colorbar(scatter2, ax=axes[1])
    
    # Plot combined metabolic activity
    scatter3 = axes[2].scatter(adata.obs['x'], adata.obs['y'], 
                              c=adata.obs['metabolic_activity'], 
                              cmap='coolwarm', s=10, alpha=0.7)
    axes[2].set_title('Combined Metabolic Activity')
    axes[2].set_xlabel('X Position')
    axes[2].set_ylabel('Y Position')
    plt.colorbar(scatter3, ax=axes[2])
    
    plt.tight_layout()
    plt.savefig('visium_metabolic_scores.png', dpi=150, bbox_inches='tight')
    print(f"   ✓ Visualization saved as 'visium_metabolic_scores.png'")
    
    # Show correlation between pathways
    print(f"\n5. Pathway correlations:")
    glycolysis_oxphos_corr = np.corrcoef(adata.obs['glycolysis_score'], adata.obs['oxphos_score'])[0, 1]
    print(f"   - Glycolysis vs OXPHOS correlation: {glycolysis_oxphos_corr:.3f}")
    
    # Find top metabolic spots
    top_metabolic_spots = adata.obs.nlargest(5, 'metabolic_activity')
    print(f"\n6. Top 5 most metabolically active spots:")
    for i, (spot_id, row) in enumerate(top_metabolic_spots.iterrows(), 1):
        print(f"   {i}. {spot_id}: activity={row['metabolic_activity']:.3f}, "
              f"glycolysis={row['glycolysis_score']:.3f}, oxphos={row['oxphos_score']:.3f}")
    
    print(f"\n✅ Example completed successfully!")
    print(f"\nTo load your own Visium data, use:")
    print(f"   adata = load_visium(adata_path='/path/to/your/visium/data', use_demo=False)")


if __name__ == "__main__":
    main()
