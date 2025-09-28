# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial release of scflux_spatial
- Spatial reaction-diffusion modeling with FiPy
- Flux balance analysis integration with COBRApy
- Self-Organizing Adaptive (SOA) loop for RD↔FBA coupling
- Visium spatial transcriptomics data support
- Human-GEM metabolic model integration
- Expression integration methods (E-Flux, iMAT-like, linear)
- Unit conversion utilities for FBA and spatial modeling
- Interactive Jupyter notebooks with examples
- Command-line interface for spatial simulations
- Comprehensive test suite
- Pre-commit hooks for code quality
- GitHub Actions CI/CD pipeline

### Changed
- N/A

### Deprecated
- N/A

### Removed
- N/A

### Fixed
- N/A

### Security
- N/A

## [0.1.0] - 2024-09-27

### Added
- Core spatial modeling functionality
- Reaction-diffusion solver with implicit time integration
- Spatial FBA coupler with under-relaxation
- Expression integration for metabolic modeling
- Visium data loader with demo dataset
- Unit conversion between FBA and spatial scales
- Comprehensive documentation and examples
- Test suite with unit and integration tests
- Poetry-based dependency management
- Pre-commit hooks for code formatting and linting
- GitHub Actions workflow for continuous integration

### Technical Details
- **RDSolver**: 2D reaction-diffusion solver using FiPy
- **SpatialFBACoupler**: SOA loop implementation with convergence tracking
- **ExpressionIntegrator**: Gene expression integration methods
- **Unit Conversion**: FBA fluxes ↔ volumetric sources
- **Cache System**: Efficient data caching for large models
- **Timeout Protection**: Prevents infinite loading of large models

### Dependencies
- Python 3.12+
- FiPy 3.5.0+ (PDE solving)
- COBRApy 0.28.0+ (metabolic modeling)
- NumPy 2.0.0+ (numerical computing)
- SciPy (scientific computing)
- Matplotlib (plotting)
- Jupyter (notebooks)
- Poetry (dependency management)

### Performance
- Optimized for 64×64 spatial grids
- Efficient sparse matrix operations
- Adaptive time stepping for stability
- Under-relaxation for convergence
- Memory-efficient caching system

### Documentation
- Comprehensive README with quickstart guide
- API documentation with examples
- Interactive Jupyter notebooks
- Contributing guidelines
- Code style guidelines
- Test documentation

### Testing
- Unit tests for core functionality
- Integration tests for component interactions
- Notebook tests for reproducibility
- Performance benchmarks
- Memory usage monitoring

### CI/CD
- GitHub Actions workflow
- Automated testing on Python 3.12
- Code formatting with black
- Linting with ruff
- Import sorting with isort
- Notebook execution testing

---

## Release Notes

### v0.1.0 - Initial Release

This is the initial release of scflux_spatial, a comprehensive package for spatial flux balance analysis. The package provides:

- **Spatial Modeling**: 2D reaction-diffusion PDEs for oxygen and glucose transport
- **Metabolic Integration**: COBRApy-based flux balance analysis with Human-GEM
- **Expression Integration**: Multiple methods for incorporating gene expression data
- **Coupling Framework**: Self-Organizing Adaptive loop for spatial-metabolic coupling
- **Data Support**: Direct integration with 10x Genomics Visium data
- **Visualization**: Interactive plots and UMAP embeddings
- **Performance**: Optimized for large-scale spatial simulations

### Key Features

1. **RDSolver**: Implicit time integration for reaction-diffusion equations
2. **SpatialFBACoupler**: Iterative coupling with convergence tracking
3. **ExpressionIntegrator**: E-Flux, iMAT-like, and linear integration methods
4. **Unit Conversion**: Seamless conversion between FBA and spatial scales
5. **Cache System**: Efficient handling of large metabolic models
6. **CLI Interface**: Command-line tools for batch processing
7. **Jupyter Integration**: Interactive notebooks with examples

### Scientific Applications

- **Tumor Metabolism**: Spatial analysis of cancer cell metabolism
- **Tissue Engineering**: Metabolic modeling in engineered tissues
- **Drug Discovery**: Spatial drug distribution and metabolism
- **Systems Biology**: Multi-scale modeling of biological systems

### Performance Benchmarks

- **Grid Size**: 64×64 spatial grids (adjustable)
- **Time Step**: Adaptive stepping for numerical stability
- **Memory**: Efficient sparse matrix operations
- **Convergence**: Under-relaxation for stable coupling
- **Caching**: Fast model loading and data persistence

### Getting Started

```bash
# Install with Poetry
poetry install

# Run examples
poetry run jupyter notebook notebooks/

# Run CLI
poetry run python -m scflux_spatial.cli.run_spatial
```

### Documentation

- **README.md**: Comprehensive overview and quickstart
- **CONTRIBUTING.md**: Guidelines for contributors
- **API Docs**: Detailed function and class documentation
- **Notebooks**: Interactive examples and tutorials
- **Tests**: Comprehensive test suite with examples

### Community

- **GitHub**: Issues, discussions, and pull requests
- **Documentation**: Comprehensive guides and examples
- **Support**: Active community support
- **Contributing**: Open to contributions and improvements

---

*For more information, see the [README.md](README.md) and [CONTRIBUTING.md](CONTRIBUTING.md) files.*

