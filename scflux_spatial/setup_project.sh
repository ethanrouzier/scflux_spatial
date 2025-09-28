#!/bin/bash

# Setup script for scflux_spatial project
echo "🚀 Setting up scflux_spatial project..."

# Check if Poetry is installed
if ! command -v poetry &> /dev/null; then
    echo "❌ Poetry is not installed. Please install Poetry first:"
    echo "   curl -sSL https://install.python-poetry.org | python3 -"
    exit 1
fi

# Install dependencies
echo "📦 Installing dependencies with Poetry..."
poetry install

# Install pre-commit hooks
echo "🔧 Installing pre-commit hooks..."
poetry run pre-commit install

# Run pre-commit on all files
echo "🎨 Running pre-commit on all files..."
poetry run pre-commit run --all-files

# Run tests
echo "🧪 Running tests..."
poetry run pytest

echo "✅ Project setup completed successfully!"
echo ""
echo "To activate the environment, run:"
echo "   poetry shell"
echo ""
echo "To run the Streamlit app:"
echo "   streamlit run app/streamlit_app.py"
echo ""
echo "To run flux analysis:"
echo "   poetry run python -m cli.run_flux"
echo ""
echo "To run spatial simulation:"
echo "   poetry run python -m cli.run_spatial"
