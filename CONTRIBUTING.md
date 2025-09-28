# Contributing to scflux_spatial

Thank you for your interest in contributing to scflux_spatial! This document provides guidelines for contributing to the project.

## 🚀 Getting Started

### Development Setup

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/yourusername/scflux_spatial.git
   cd scflux_spatial
   ```

3. **Install development dependencies**:
   ```bash
   poetry install --with dev
   ```

4. **Install pre-commit hooks**:
   ```bash
   poetry run pre-commit install
   ```

### Development Workflow

1. **Create a feature branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** and ensure they pass tests:
   ```bash
   poetry run pytest
   poetry run black .
   poetry run ruff check .
   poetry run isort .
   ```

3. **Commit your changes**:
   ```bash
   git add .
   git commit -m "Add your feature description"
   ```

4. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

5. **Create a Pull Request** on GitHub

## 📝 Code Style

### Python Code

- **Formatting**: Use `black` with line length 100
- **Linting**: Use `ruff` for fast linting
- **Import sorting**: Use `isort`
- **Type hints**: Use type hints for function signatures
- **Docstrings**: Use Google-style docstrings

### Example

```python
def flux_to_volumetric_source(
    v_mmol_gDW_h: float, 
    rho_gDW_per_m3: float
) -> float:
    """
    Convert FBA flux to volumetric source.
    
    Args:
        v_mmol_gDW_h: Flux in mmol·gDW⁻¹·h⁻¹
        rho_gDW_per_m3: Biomass density in gDW·m⁻³
        
    Returns:
        Volumetric source in mol·m⁻³·s⁻¹
    """
    return float(v_mmol_gDW_h) * float(rho_gDW_per_m3) * 1e-3 / 3600.0
```

## 🧪 Testing

### Running Tests

```bash
# Run all tests
poetry run pytest

# Run specific test files
poetry run pytest tests/test_units.py

# Run with coverage
poetry run pytest --cov=scflux_spatial

# Run notebooks as tests
poetry run pytest --nbmake notebooks/
```

### Writing Tests

- **Unit tests**: Test individual functions and classes
- **Integration tests**: Test component interactions
- **Notebook tests**: Ensure notebooks run without errors

### Test Structure

```python
def test_flux_conversion():
    """Test flux to volumetric source conversion."""
    result = flux_to_volumetric_source(1.0, 1e6)
    expected = 1.0 * 1e6 * 1e-3 / 3600.0
    assert abs(result - expected) < 1e-12
```

## 📚 Documentation

### Code Documentation

- **Docstrings**: All public functions and classes need docstrings
- **Type hints**: Use type hints for better IDE support
- **Examples**: Include usage examples in docstrings

### Notebook Documentation

- **Clear explanations**: Explain the scientific background
- **Reproducibility**: Use random seeds and cache data
- **Visualization**: Include plots and figures
- **Performance**: Add timing information for long operations

## 🐛 Bug Reports

When reporting bugs, please include:

1. **Environment**: Python version, OS, package versions
2. **Reproduction**: Minimal code to reproduce the issue
3. **Expected behavior**: What should happen
4. **Actual behavior**: What actually happens
5. **Error messages**: Full traceback if applicable

## 💡 Feature Requests

When requesting features, please include:

1. **Use case**: Why is this feature needed?
2. **Proposed solution**: How should it work?
3. **Alternatives**: Other approaches considered
4. **Implementation**: Any ideas for implementation

## 🔬 Scientific Contributions

We welcome contributions in:

- **Algorithms**: New methods for spatial metabolic modeling
- **Integration**: Support for new data types and formats
- **Visualization**: New plotting and analysis tools
- **Performance**: Optimizations for large-scale simulations
- **Documentation**: Tutorials and examples

## 📋 Pull Request Checklist

Before submitting a PR, ensure:

- [ ] Code follows the style guidelines
- [ ] Tests pass (`poetry run pytest`)
- [ ] Code is formatted (`poetry run black .`)
- [ ] Linting passes (`poetry run ruff check .`)
- [ ] Imports are sorted (`poetry run isort .`)
- [ ] Documentation is updated
- [ ] Notebooks run without errors
- [ ] Commit messages are clear and descriptive

## 🏷️ Release Process

1. **Version bump**: Update version in `pyproject.toml`
2. **Changelog**: Update `CHANGELOG.md`
3. **Tag**: Create a git tag for the release
4. **Build**: Test the build process
5. **Publish**: Publish to PyPI

## 📞 Getting Help

- **GitHub Issues**: For bug reports and feature requests
- **GitHub Discussions**: For questions and general discussion
- **Email**: your.email@example.com

## 🙏 Recognition

Contributors will be recognized in:
- **README.md**: Listed as contributors
- **CHANGELOG.md**: Mentioned in release notes
- **GitHub**: Listed as contributors

Thank you for contributing to scflux_spatial! 🎉

