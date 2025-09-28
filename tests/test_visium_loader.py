"""
Tests for Visium data loader module.
"""

import pytest
import numpy as np
import anndata as ad
from pathlib import Path

from scflux_spatial.dataio import VisiumLoader, load_visium


class TestVisiumLoader:
    """Test cases for VisiumLoader class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.loader = VisiumLoader()

    def test_load_visium_demo(self):
        """Test loading demo dataset."""
        adata = self.loader.load_visium(use_demo=True)
        
        assert isinstance(adata, ad.AnnData)
        assert adata.n_obs > 0
        assert adata.n_vars > 0
        
        # Check that spatial coordinates are present
        assert 'spatial_coordinates' in adata.obsm
        assert 'x' in adata.obs
        assert 'y' in adata.obs
        
        # Check that metabolic scores are calculated
        assert 'glycolysis_score' in adata.obs
        assert 'oxphos_score' in adata.obs
        assert 'metabolic_activity' in adata.obs

    def test_load_visium_mock_fallback(self):
        """Test mock dataset creation when demo fails."""
        # This will test the mock dataset creation
        adata = self.loader._create_mock_dataset()
        
        assert isinstance(adata, ad.AnnData)
        assert adata.n_obs == 1000
        assert adata.n_vars == 2000
        assert 'spatial_coordinates' in adata.obsm
        assert 'x' in adata.obs
        assert 'y' in adata.obs

    def test_preprocess_visium_data(self):
        """Test Visium data preprocessing."""
        # Create mock data
        adata = self.loader._create_mock_dataset()
        
        # Preprocess
        processed_adata = self.loader._preprocess_visium_data(adata)
        
        assert isinstance(processed_adata, ad.AnnData)
        assert processed_adata.n_obs <= adata.n_obs  # Some cells might be filtered
        assert processed_adata.n_vars <= adata.n_vars  # Some genes might be filtered
        assert 'mt' in processed_adata.var.columns

    def test_add_metabolic_scores(self):
        """Test metabolic pathway score calculation."""
        # Create mock data with some metabolic genes
        n_spots = 100
        n_genes = 50
        
        # Create expression matrix
        X = np.random.negative_binomial(5, 0.3, size=(n_spots, n_genes)).astype(float)
        
        # Create gene names including some metabolic genes
        gene_names = ['HK1', 'GAPDH', 'LDHA', 'NDUFA1', 'SDHA', 'COX4I1'] + [f"GENE_{i}" for i in range(44)]
        
        # Create AnnData
        adata = ad.AnnData(X=X)
        adata.var_names = gene_names
        adata.obsm['spatial_coordinates'] = np.random.uniform(0, 10, size=(n_spots, 2))
        
        # Add metabolic scores
        scored_adata = self.loader._add_metabolic_scores(adata)
        
        assert 'glycolysis_score' in scored_adata.obs
        assert 'oxphos_score' in scored_adata.obs
        assert 'metabolic_activity' in scored_adata.obs
        assert 'x' in scored_adata.obs
        assert 'y' in scored_adata.obs
        
        # Check that scores are calculated
        assert not scored_adata.obs['glycolysis_score'].isna().all()
        assert not scored_adata.obs['oxphos_score'].isna().all()

    def test_load_visium_invalid_path(self):
        """Test loading with invalid path."""
        with pytest.raises(ValueError):
            self.loader.load_visium(adata_path=None, use_demo=False)

    def test_spatial_coordinates_handling(self):
        """Test spatial coordinates handling."""
        # Create mock data without spatial coordinates
        n_spots = 50
        n_genes = 100
        X = np.random.negative_binomial(5, 0.3, size=(n_spots, n_genes)).astype(float)
        
        adata = ad.AnnData(X=X)
        adata.var_names = [f"GENE_{i}" for i in range(n_genes)]
        
        # Add metabolic scores (should create spatial coordinates)
        scored_adata = self.loader._add_metabolic_scores(adata)
        
        assert 'spatial_coordinates' in scored_adata.obsm
        assert 'x' in scored_adata.obs
        assert 'y' in scored_adata.obs
        assert scored_adata.obs['x'].shape[0] == n_spots
        assert scored_adata.obs['y'].shape[0] == n_spots


class TestLoadVisiumFunction:
    """Test cases for the convenience function load_visium."""

    def test_load_visium_function_demo(self):
        """Test the convenience function with demo data."""
        adata = load_visium(use_demo=True)
        
        assert isinstance(adata, ad.AnnData)
        assert adata.n_obs > 0
        assert adata.n_vars > 0
        assert 'glycolysis_score' in adata.obs
        assert 'oxphos_score' in adata.obs

    def test_load_visium_function_invalid_path(self):
        """Test the convenience function with invalid path."""
        with pytest.raises(ValueError):
            load_visium(adata_path=None, use_demo=False)

    def test_load_visium_function_parameters(self):
        """Test the convenience function parameter handling."""
        # Test with use_demo=True (should work)
        adata1 = load_visium(use_demo=True)
        assert isinstance(adata1, ad.AnnData)
        
        # Test with use_demo=False but no path (should raise error)
        with pytest.raises(ValueError):
            load_visium(use_demo=False)


def test_metabolic_gene_sets():
    """Test that metabolic gene sets are properly defined."""
    loader = VisiumLoader()
    
    # Create mock data with all metabolic genes
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
    
    all_genes = glycolysis_genes + oxphos_genes + [f"OTHER_{i}" for i in range(100)]
    
    n_spots = 50
    n_genes = len(all_genes)
    X = np.random.negative_binomial(5, 0.3, size=(n_spots, n_genes)).astype(float)
    
    adata = ad.AnnData(X=X)
    adata.var_names = all_genes
    adata.obsm['spatial_coordinates'] = np.random.uniform(0, 10, size=(n_spots, 2))
    
    # Add metabolic scores
    scored_adata = loader._add_metabolic_scores(adata)
    
    # Check that scores are calculated for all spots
    assert scored_adata.obs['glycolysis_score'].shape[0] == n_spots
    assert scored_adata.obs['oxphos_score'].shape[0] == n_spots
    assert scored_adata.obs['metabolic_activity'].shape[0] == n_spots
    
    # Check that scores are not all zero (at least some metabolic genes present)
    assert scored_adata.obs['glycolysis_score'].abs().sum() > 0
    assert scored_adata.obs['oxphos_score'].abs().sum() > 0
