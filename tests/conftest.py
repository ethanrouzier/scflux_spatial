"""
Configuration and fixtures for pytest tests.

This module provides common fixtures and configuration for all tests.
"""

import pytest
import numpy as np
import tempfile
import os
from pathlib import Path


@pytest.fixture
def random_seed():
    """Set random seed for reproducible tests."""
    np.random.seed(42)
    return 42


@pytest.fixture
def temp_dir():
    """Create temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture
def mock_gene_expression():
    """Mock gene expression data for testing."""
    return {
        'GENE_A': 5.0,
        'GENE_B': 3.0,
        'GENE_C': 2.0,
        'GENE_D': 4.0,
        'GENE_E': 1.0,
        'GENE_F': 6.0,
        'HK1': 4.5,
        'GAPDH': 3.8,
        'LDHA': 2.1,
        'PFKL': 3.2
    }


@pytest.fixture
def mock_gpr_rules():
    """Mock GPR rules for testing."""
    return {
        'simple_and': 'GENE_A and GENE_B',
        'simple_or': 'GENE_C or GENE_D',
        'complex_nested': '(GENE_A and GENE_B) or (GENE_C and GENE_D)',
        'triple_and': 'GENE_A and GENE_B and GENE_C',
        'triple_or': 'GENE_A or GENE_B or GENE_C',
        'mixed_operators': 'GENE_A and (GENE_B or GENE_C) and GENE_D',
        'single_gene': 'GENE_A',
        'glycolysis_rule': 'HK1 and GAPDH and PFKL'
    }


@pytest.fixture
def mock_flux_bounds():
    """Mock flux bounds for testing."""
    return {
        'R1': (0, 10),
        'R2': (0, 10),
        'R3': (0, 10),
        'R4': (0, 10),
        'HEX1': (0, 1000),
        'GAPD': (0, 1000),
        'LDH_L': (0, 1000),
        'PFK': (0, 1000),
        'PYK': (0, 1000)
    }


@pytest.fixture
def mock_toy_model():
    """Mock toy metabolic model for testing."""
    class MockToyModel:
        def __init__(self):
            self.reactions = ['R1', 'R2', 'R3', 'R4']
            self.metabolites = ['A', 'B', 'C', 'D']
            
            # Simple linear pathway: A -> B -> C -> D
            self.S = np.array([
                [-1,  0,  0,  0],  # A: consumed by R1
                [ 1, -1,  0,  0],  # B: produced by R1, consumed by R2
                [ 0,  1, -1,  0],  # C: produced by R2, consumed by R3
                [ 0,  0,  1, -1],  # D: produced by R3, consumed by R4
            ])
            
            self.bounds = np.array([
                [0, 10],    # R1: A -> B
                [0, 10],    # R2: B -> C
                [0, 10],    # R3: C -> D
                [0, 10],    # R4: D -> (export)
            ])
            
            # Objective: maximize D production (R4)
            self.objective = np.array([0, 0, 0, 1])
    
    return MockToyModel()


@pytest.fixture
def mock_spatial_coordinates():
    """Mock spatial coordinates for testing."""
    return {
        'spot_001': (0.1, 0.2),
        'spot_002': (0.3, 0.4),
        'spot_003': (0.5, 0.6),
        'spot_004': (0.7, 0.8),
        'spot_005': (0.9, 1.0)
    }


@pytest.fixture
def mock_concentration_field():
    """Mock concentration field for testing."""
    grid_size = 50
    # Create a simple concentration field with gradient
    x = np.linspace(0, 1, grid_size)
    y = np.linspace(0, 1, grid_size)
    X, Y = np.meshgrid(x, y)
    
    # Linear gradient from 0 to 1
    concentration = X + Y
    return concentration


@pytest.fixture
def mock_reaction_rate_field():
    """Mock reaction rate field for testing."""
    grid_size = 50
    # Create a simple reaction rate field
    reaction_rate = np.ones((grid_size, grid_size)) * 1e-6
    return reaction_rate


@pytest.fixture
def mock_rd_parameters():
    """Mock RD parameters for testing."""
    return {
        'grid_size': 50,
        'domain_size': 1.0,  # mm
        'diffusion_coefficient': 1e-9,  # m²/s
        'boundary_condition': 'dirichlet',
        'boundary_value': 1e-3,  # mol/L
        'convergence_tolerance': 1e-4,
        'max_iterations': 100
    }


@pytest.fixture
def mock_kinetics_parameters():
    """Mock kinetics parameters for testing."""
    return {
        'cell_density': 1e9,      # cells/L
        'spot_density': 1e6,      # spots/L
        'cell_volume': 1e-12,     # L/cell
        'spot_volume': 1e-9,      # L/spot
        'mmol_per_mol': 1e-3,
        'h_per_s': 3600,
        'gDW_per_cell': 1e-12
    }


@pytest.fixture
def mock_michaelis_menten_params():
    """Mock Michaelis-Menten parameters for testing."""
    return {
        'O2': {'Vmax': 1e-6, 'Km': 1e-5},
        'Glc': {'Vmax': 1e-5, 'Km': 1e-4},
        'Lac': {'Vmax': 1e-6, 'Km': 1e-5}
    }


@pytest.fixture
def mock_soa_parameters():
    """Mock SOA loop parameters for testing."""
    return {
        'convergence_tolerance': 1e-4,
        'max_iterations': 10,
        'use_pfba': True,
        'relaxation_factor': 0.1
    }


@pytest.fixture
def mock_visium_data():
    """Mock Visium data structure for testing."""
    import pandas as pd
    import numpy as np
    
    # Create mock AnnData-like structure
    n_spots = 100
    n_genes = 1000
    
    # Mock expression matrix
    expression_matrix = np.random.poisson(5, (n_spots, n_genes))
    
    # Mock spatial coordinates
    spatial_coords = np.random.uniform(0, 10, (n_spots, 2))
    
    # Mock gene names
    gene_names = [f'GENE_{i:04d}' for i in range(n_genes)]
    
    # Mock spot names
    spot_names = [f'spot_{i:03d}' for i in range(n_spots)]
    
    # Create mock metadata
    metadata = pd.DataFrame({
        'x': spatial_coords[:, 0],
        'y': spatial_coords[:, 1],
        'total_counts': expression_matrix.sum(axis=1),
        'n_genes_by_counts': (expression_matrix > 0).sum(axis=1)
    }, index=spot_names)
    
    return {
        'expression_matrix': expression_matrix,
        'spatial_coords': spatial_coords,
        'gene_names': gene_names,
        'spot_names': spot_names,
        'metadata': metadata
    }


@pytest.fixture(scope="session")
def test_data_dir():
    """Create test data directory for session-scoped fixtures."""
    test_dir = Path(__file__).parent / "test_data"
    test_dir.mkdir(exist_ok=True)
    return test_dir


# Test configuration
def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
    config.addinivalue_line(
        "markers", "unit: marks tests as unit tests"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test collection to add markers automatically."""
    for item in items:
        # Add unit marker to all tests by default
        if not any(marker in item.keywords for marker in ["integration", "slow"]):
            item.add_marker(pytest.mark.unit)
        
        # Add slow marker to tests that take longer
        if "performance" in item.name or "large" in item.name:
            item.add_marker(pytest.mark.slow)
