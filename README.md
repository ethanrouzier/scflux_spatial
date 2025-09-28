# scflux_spatial

**Spatial Flux Balance Analysis for Single-Cell Data**

A Python package for integrating spatial transcriptomics data with metabolic modeling through reaction-diffusion coupling and flux balance analysis.

## Features

- **Spatial Reaction-Diffusion Modeling**: 2D PDE solving with FiPy for oxygen and glucose transport
- **Flux Balance Analysis Integration**: COBRApy-based metabolic modeling with Human-GEM
- **Expression Integration**: E-Flux, iMAT-like, and linear methods for gene expression integration
- **Self-Organizing Adaptive (SOA) Loop**: Iterative coupling between spatial and metabolic processes
- **Visium Data Support**: Direct integration with 10x Genomics Visium spatial transcriptomics
- **Interactive Visualization**: UMAP embeddings colored by predicted metabolic fluxes

## Installation

### Prerequisites

- Python 3.12+
- Poetry (recommended) or pip

### Install with Poetry (Recommended)

```bash
# Clone the repository
git clone https://github.com/ethanrouzier/scflux_spatial.git
cd scflux_spatial

# Install dependencies
poetry install

# Activate the environment
poetry shell
```

### Install with pip

```bash
# Clone the repository
git clone https://github.com/ethanrouzier/scflux_spatial.git
cd scflux_spatial

# Install dependencies
pip install -e .
```

## Quickstart

### 1. Basic Usage

```python
import scflux_spatial
from scflux_spatial.dataio import load_visium
from scflux_spatial.spatial.rd import RDSolver, Species
from scflux_spatial.spatial.coupling import SpatialFBACoupler

# Load spatial transcriptomics data
adata = load_visium(use_demo=True)

# Set up reaction-diffusion solver
rd_solver = RDSolver(
    species=[
        Species(name="O2", D=1e-9, initial_conc=2e-5),
        Species(name="Glc", D=1e-9, initial_conc=5e-3)
    ],
    domain_size=(2.0, 2.0),  # mm
    grid_size=(64, 64)
)

# Set up spatial FBA coupling
coupler = SpatialFBACoupler(
    rd=rd_solver,
    fba_fn=your_fba_function,
    rho_gDW_per_m3=1e6,  # Biomass density
    alpha=0.4  # Under-relaxation parameter
)

# Run coupling simulation
for iteration in range(20):
    concentrations = coupler.iterate(dt=0.1)
    print(f"Iteration {iteration}: Max O2 = {concentrations['O2'].max():.2e} mol/L")
```

### 2. Jupyter Notebooks

Explore the comprehensive examples:

```bash
# Launch Jupyter
poetry run jupyter notebook notebooks/

# Or with JupyterLab
poetry run jupyter lab notebooks/
```

**Available Notebooks:**
- `01_quickstart_visium.ipynb`: Introduction to Visium data loading and preprocessing
- `02_flux_validation.ipynb`: Flux balance analysis validation with hypoxia correlations
- `03_spatial_coupling_demo.ipynb`: Complete spatial coupling simulation with SOA loop

### 3. Command Line Interface

```bash
# Run spatial simulation
poetry run python -m scflux_spatial.cli.run_spatial --grid-size 32 --iterations 10

# Run flux analysis
poetry run python -m scflux_spatial.cli.run_flux --method eflux
```

## Scientific Background

This package implements a novel approach to spatial metabolic modeling by coupling:

1. **Reaction-Diffusion PDEs**: Transport of oxygen and glucose through tissue
2. **Flux Balance Analysis**: Metabolic flux predictions based on gene expression
3. **Self-Organizing Adaptive Loop**: Iterative coupling for convergence

### Key Components

- **RDSolver**: 2D reaction-diffusion solver using FiPy finite volume methods
- **SpatialFBACoupler**: Orchestrates the SOA loop between spatial and metabolic processes
- **ExpressionIntegrator**: Integrates gene expression data into metabolic models
- **Unit Conversion**: Converts between FBA fluxes (mmol·gDW⁻¹·h⁻¹) and volumetric sources (mol·m⁻³·s⁻¹)

## Data Integration

### Supported Data Types

- **10x Genomics Visium**: Spatial transcriptomics with spot coordinates
- **Human-GEM**: Genome-scale metabolic model
- **Gene Expression**: Integration via E-Flux, iMAT-like, or linear methods

## Testing

Run the test suite:

```bash
# Run all tests
poetry run pytest

# Run specific test categories
poetry run pytest tests/test_units.py
poetry run pytest tests/test_rd_diffusion.py
poetry run pytest tests/test_coupling_dummy.py

# Run notebooks as tests
poetry run pytest --nbmake notebooks/
```

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Install development dependencies
poetry install --with dev

# Install pre-commit hooks
poetry run pre-commit install

# Run linting
poetry run black .
poetry run ruff check .
poetry run isort .
```

## Architecture

```
scflux_spatial/
├── spatial/          # Reaction-diffusion modeling
│   ├── rd.py         # RDSolver, Species, RDField
│   ├── coupling.py   # SpatialFBACoupler
│   ├── kinetics.py   # KineticsModel
│   └── units.py      # Unit conversion utilities
├── fba/              # Flux balance analysis
│   ├── core.py       # FBA core functions
│   ├── integrate_expression.py  # Expression integration
│   └── objectives.py # Objective functions
├── gem/              # Genome-scale models
│   ├── human_gem.py  # Human-GEM integration
│   └── gpr.py        # Gene-protein-reaction rules
├── dataio/           # Data loading
│   └── visium_loader.py  # Visium data support
├── viz/              # Visualization
│   ├── maps.py       # Spatial mapping
│   └── escher_view.py # Escher pathway visualization
└── cli/              # Command line interface
    ├── run_spatial.py # Spatial simulation CLI
    └── run_flux.py   # Flux analysis CLI
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this software in your research, please cite:

```bibtex
@software{scflux_spatial,
  title={scflux_spatial: Spatial Flux Balance Analysis for Single-Cell Data},
  author={Ethan Rouzier},
  year={2024},
  url={https://github.com/ethanrouzier/scflux_spatial}
}
```

## Support

- **Issues**: [GitHub Issues](https://github.com/ethanrouzier/scflux_spatial/issues)
- **Discussions**: [GitHub Discussions](https://github.com/ethanrouzier/scflux_spatial/discussions)