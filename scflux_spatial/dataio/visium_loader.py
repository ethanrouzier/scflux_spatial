"""
Visium data loader for spatial transcriptomics data.

This module provides functionality to load and preprocess 10X Visium data
using squidpy and scanpy.

Dataset Sources:
- Demo dataset: Mouse brain Visium data from 10X Genomics
  https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-1-sagittal-posterior-1-standard-1-1-0
- Tutorial: Squidpy H&E tutorial
  https://squidpy.readthedocs.io/en/stable/auto_tutorials/tutorial_spatialgeneexpression.html
- 10X Visium documentation:
  https://support.10xgenomics.com/spatial-gene-expression/datasets
"""

import os
from pathlib import Path
from typing import Optional, Tuple, Union

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from rich.console import Console

console = Console()


class VisiumLoader:
    """Loader for 10X Visium spatial transcriptomics data."""

    def __init__(
        self,
        data_path: Union[str, Path],
        library_id: str = "spaceranger_output",
        spatial_key: str = "spatial",
    ):
        """
        Initialize Visium data loader.

        Args:
            data_path: Path to the Visium data directory
            library_id: Library ID for the dataset
            spatial_key: Key for spatial coordinates in AnnData
        """
        self.data_path = Path(data_path)
        self.library_id = library_id
        self.spatial_key = spatial_key
        self.adata: Optional[ad.AnnData] = None

    def load_data(
        self,
        filter_missing_genes: bool = True,
        min_cells_per_gene: int = 10,
        min_genes_per_cell: int = 200,
    ) -> ad.AnnData:
        """
        Load Visium data using squidpy.

        Args:
            filter_missing_genes: Whether to filter genes with low expression
            min_cells_per_gene: Minimum number of cells expressing a gene
            min_genes_per_cell: Minimum number of genes per cell

        Returns:
            AnnData object with spatial information
        """
        console.print(f"Loading Visium data from {self.data_path}")

        # Check if data path exists
        if not self.data_path.exists():
            raise FileNotFoundError(f"Data path {self.data_path} does not exist")

        # Load data using squidpy
        self.adata = sq.read.visium(self.data_path, library_id=self.library_id)

        # Store spatial key
        self.adata.uns[self.spatial_key] = {
            "library_id": self.library_id,
            "use_quality": "hires",
        }

        # Basic preprocessing
        if filter_missing_genes:
            self._filter_genes_and_cells(min_cells_per_gene, min_genes_per_cell)

        console.print(f"Loaded data with shape: {self.adata.shape}")
        return self.adata

    def _filter_genes_and_cells(
        self, min_cells_per_gene: int, min_genes_per_cell: int
    ) -> None:
        """Filter genes and cells based on expression thresholds."""
        console.print("Filtering genes and cells...")

        # Filter genes
        sc.pp.filter_genes(self.adata, min_cells=min_cells_per_gene)
        console.print(f"After gene filtering: {self.adata.shape}")

        # Filter cells
        sc.pp.filter_cells(self.adata, min_genes=min_genes_per_cell)
        console.print(f"After cell filtering: {self.adata.shape}")

    def get_spatial_coordinates(self) -> np.ndarray:
        """
        Extract spatial coordinates from the loaded data.

        Returns:
            Array of spatial coordinates (n_cells, 2)
        """
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        # Get spatial coordinates
        spatial_coords = self.adata.obsm[f"{self.spatial_key}_coordinates"]
        return spatial_coords

    def get_gene_expression_matrix(self) -> np.ndarray:
        """
        Get the gene expression matrix.

        Returns:
            Gene expression matrix (n_cells, n_genes)
        """
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        return self.adata.X.toarray() if hasattr(self.adata.X, "toarray") else self.adata.X

    def get_gene_names(self) -> np.ndarray:
        """
        Get gene names.

        Returns:
            Array of gene names
        """
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        return self.adata.var_names.values

    def get_cell_barcodes(self) -> np.ndarray:
        """
        Get cell barcodes.

        Returns:
            Array of cell barcodes
        """
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        return self.adata.obs_names.values

    def add_metadata(self, metadata: pd.DataFrame) -> None:
        """
        Add metadata to the AnnData object.

        Args:
            metadata: DataFrame with cell metadata
        """
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        # Align metadata with cell barcodes
        common_barcodes = self.adata.obs_names.intersection(metadata.index)
        metadata_aligned = metadata.loc[common_barcodes]

        # Add to obs
        for col in metadata_aligned.columns:
            self.adata.obs[col] = metadata_aligned[col]

        console.print(f"Added metadata columns: {list(metadata_aligned.columns)}")

    def save_data(self, output_path: Union[str, Path]) -> None:
        """
        Save the loaded data.

        Args:
            output_path: Path to save the data
        """
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        self.adata.write(output_path)
        console.print(f"Saved data to {output_path}")

    def get_tissue_image(self) -> Optional[np.ndarray]:
        """
        Get the tissue image if available.

        Returns:
            Tissue image array or None if not available
        """
        if self.adata is None:
            raise ValueError("Data not loaded. Call load_data() first.")

        # Check if image is available in uns
        if self.spatial_key in self.adata.uns:
            spatial_info = self.adata.uns[self.spatial_key]
            if "images" in spatial_info:
                return spatial_info["images"]["hires"]

        return None

    def load_visium(
        self, 
        adata_path: Optional[str] = None, 
        use_demo: bool = True
    ) -> ad.AnnData:
        """
        Load Visium data with preprocessing and metabolic pathway scoring.

        Args:
            adata_path: Path to Visium data directory (if use_demo=False)
            use_demo: If True, load demo mouse brain dataset from squidpy

        Returns:
            AnnData object with preprocessed Visium data and metabolic scores

        Notes:
            Demo dataset: Mouse brain Visium data from 10X Genomics
            https://www.10xgenomics.com/resources/datasets/mouse-brain-serial-section-1-sagittal-posterior-1-standard-1-1-0
            
            Tutorial reference: Squidpy H&E tutorial
            https://squidpy.readthedocs.io/en/stable/auto_tutorials/tutorial_spatialgeneexpression.html
        """
        if use_demo:
            console.print("Loading demo mouse brain Visium dataset...")
            adata = self._load_demo_dataset()
        else:
            if adata_path is None:
                raise ValueError("adata_path must be provided when use_demo=False")
            console.print(f"Loading Visium data from {adata_path}")
            adata = self._load_visium_directory(adata_path)

        # Preprocess the data
        adata = self._preprocess_visium_data(adata)
        
        # Add metabolic pathway scores
        adata = self._add_metabolic_scores(adata)
        
        console.print(f"Loaded and preprocessed data with shape: {adata.shape}")
        return adata

    def _load_demo_dataset(self) -> ad.AnnData:
        """Load demo mouse brain Visium dataset from squidpy."""
        try:
            # Load the demo dataset from squidpy
            adata = sq.datasets.visium_fluo_adata_crop()
            console.print("Successfully loaded demo mouse brain dataset")
            return adata
        except Exception as e:
            console.print(f"Error loading demo dataset: {e}")
            console.print("Creating mock dataset for demonstration...")
            return self._create_mock_dataset()

    def _load_visium_directory(self, adata_path: str) -> ad.AnnData:
        """Load Visium data from standard directory structure."""
        adata_path = Path(adata_path)
        
        if not adata_path.exists():
            raise FileNotFoundError(f"Data path {adata_path} does not exist")
        
        try:
            # Use scanpy to read Visium data
            adata = sc.read_visium(adata_path)
            console.print("Successfully loaded Visium data from directory")
            return adata
        except Exception as e:
            console.print(f"Error loading Visium directory: {e}")
            raise

    def _preprocess_visium_data(self, adata: ad.AnnData) -> ad.AnnData:
        """Preprocess Visium data with filtering and normalization."""
        console.print("Preprocessing Visium data...")
        
        # Calculate QC metrics
        adata.var['mt'] = adata.var_names.str.startswith('MT-')  # Mitochondrial genes
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        
        # Filter cells and genes
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=10)
        
        # Filter mitochondrial genes
        adata = adata[:, ~adata.var['mt']].copy()
        
        # Normalize and log-transform
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Find highly variable genes
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata.raw = adata
        
        # Scale data
        sc.pp.scale(adata, max_value=10)
        
        console.print("Data preprocessing completed")
        return adata

    def _add_metabolic_scores(self, adata: ad.AnnData) -> ad.AnnData:
        """Add metabolic pathway scores (glycolysis and OXPHOS)."""
        console.print("Calculating metabolic pathway scores...")
        
        # Define gene sets for metabolic pathways
        glycolysis_genes = [
            'HK1', 'HK2', 'HK3', 'GPI', 'PFKL', 'PFKM', 'PFKP', 'ALDOA', 'ALDOB', 'ALDOC',
            'TPI1', 'GAPDH', 'PGK1', 'PGK2', 'PGAM1', 'PGAM2', 'ENO1', 'ENO2', 'ENO3',
            'PKM', 'PKLR', 'LDHA', 'LDHB', 'LDHC', 'LDHD'
        ]
        
        oxphos_genes = [
            'NDUFA1', 'NDUFA2', 'NDUFA3', 'NDUFA4', 'NDUFA5', 'NDUFA6', 'NDUFA7', 'NDUFA8',
            'NDUFA9', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFAB1', 'NDUFAF1',
            'NDUFAF2', 'NDUFAF3', 'NDUFAF4', 'NDUFAF5', 'NDUFAF6', 'NDUFAF7', 'NDUFAF8',
            'NDUFB1', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFB8',
            'NDUFB9', 'NDUFB10', 'NDUFB11', 'NDUFV1', 'NDUFV2', 'NDUFV3', 'SDHA', 'SDHB',
            'SDHC', 'SDHD', 'UQCRC1', 'UQCRC2', 'UQCRFS1', 'UQCRQ', 'UQCRB', 'UQCRH',
            'CYC1', 'CYCS', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C',
            'COX7A1', 'COX7A2', 'COX7B', 'COX7C', 'COX8A', 'ATP5F1A', 'ATP5F1B', 'ATP5F1C',
            'ATP5F1D', 'ATP5F1E', 'ATP5PB', 'ATP5PD', 'ATP5PF', 'ATP5PO', 'ATP5MG'
        ]
        
        # Filter genes to those present in the dataset
        available_glycolysis = [gene for gene in glycolysis_genes if gene in adata.var_names]
        available_oxphos = [gene for gene in oxphos_genes if gene in adata.var_names]
        
        console.print(f"Found {len(available_glycolysis)} glycolysis genes and {len(available_oxphos)} OXPHOS genes")
        
        # Calculate pathway scores
        if available_glycolysis:
            sc.tl.score_genes(adata, available_glycolysis, score_name='glycolysis_score')
        else:
            adata.obs['glycolysis_score'] = 0.0
            
        if available_oxphos:
            sc.tl.score_genes(adata, available_oxphos, score_name='oxphos_score')
        else:
            adata.obs['oxphos_score'] = 0.0
        
        # Add combined metabolic activity score
        adata.obs['metabolic_activity'] = adata.obs['glycolysis_score'] + adata.obs['oxphos_score']
        
        # Add spatial coordinates if not present
        if 'spatial_coordinates' not in adata.obsm:
            if 'spatial' in adata.obsm:
                adata.obsm['spatial_coordinates'] = adata.obsm['spatial']
            else:
                # Create mock spatial coordinates
                n_spots = adata.n_obs
                x_coords = np.random.uniform(0, 10, n_spots)
                y_coords = np.random.uniform(0, 10, n_spots)
                adata.obsm['spatial_coordinates'] = np.column_stack([x_coords, y_coords])
        
        # Add x, y coordinates to obs for easy access
        spatial_coords = adata.obsm['spatial_coordinates']
        adata.obs['x'] = spatial_coords[:, 0]
        adata.obs['y'] = spatial_coords[:, 1]
        
        console.print("Metabolic pathway scores calculated")
        return adata

    def _create_mock_dataset(self) -> ad.AnnData:
        """Create a mock Visium dataset for demonstration purposes."""
        console.print("Creating mock Visium dataset...")
        
        # Create mock expression matrix
        n_spots = 1000
        n_genes = 2000
        
        # Create realistic expression distribution
        X = np.random.negative_binomial(5, 0.3, size=(n_spots, n_genes)).astype(float)
        
        # Add some spatial structure
        x_coords = np.random.uniform(0, 10, n_spots)
        y_coords = np.random.uniform(0, 10, n_spots)
        
        # Create gene names
        gene_names = [f"GENE_{i:04d}" for i in range(n_genes)]
        
        # Create spot barcodes
        spot_barcodes = [f"SPOT_{i:06d}" for i in range(n_spots)]
        
        # Create AnnData object
        adata = ad.AnnData(X=X)
        adata.var_names = gene_names
        adata.obs_names = spot_barcodes
        
        # Add spatial coordinates
        adata.obsm['spatial_coordinates'] = np.column_stack([x_coords, y_coords])
        adata.obsm['spatial'] = adata.obsm['spatial_coordinates']
        
        console.print("Mock dataset created successfully")
        return adata


# Convenience function for direct usage
def load_visium(
    adata_path: Optional[str] = None, 
    use_demo: bool = True
) -> ad.AnnData:
    """
    Convenience function to load Visium data with preprocessing.
    
    Args:
        adata_path: Path to Visium data directory (if use_demo=False)
        use_demo: If True, load demo mouse brain dataset from squidpy
        
    Returns:
        AnnData object with preprocessed Visium data and metabolic scores
        
    Examples:
        >>> # Load demo dataset
        >>> adata = load_visium(use_demo=True)
        
        >>> # Load custom dataset
        >>> adata = load_visium(adata_path="/path/to/visium/data", use_demo=False)
    """
    # For demo data, we don't need a real path
    demo_path = "demo" if use_demo else adata_path
    loader = VisiumLoader(demo_path)
    return loader.load_visium(adata_path, use_demo)
