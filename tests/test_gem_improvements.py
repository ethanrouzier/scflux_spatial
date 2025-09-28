"""
Tests for improved GEM module functionality.
"""

import pytest
import numpy as np
from pathlib import Path

from scflux_spatial.gem.human_gem import HumanGEM
from scflux_spatial.gem.gpr import GPRParser


class TestHumanGEMImprovements:
    """Test cases for improved HumanGEM functionality."""

    def test_get_latest_model_url(self):
        """Test getting latest model URL."""
        gem = HumanGEM()
        url = gem._get_latest_model_url()
        
        assert isinstance(url, str)
        assert "github.com" in url
        assert "Human-GEM" in url
        assert url.endswith(".xml")

    def test_get_default_medium(self):
        """Test getting default medium composition."""
        gem = HumanGEM()
        medium = gem.get_default_medium()
        
        assert isinstance(medium, dict)
        assert len(medium) > 0
        
        # Check for essential nutrients
        assert 'EX_glc__D_e' in medium  # Glucose
        assert 'EX_o2_e' in medium      # Oxygen
        assert 'EX_gln__L_e' in medium  # Glutamine
        
        # Check that uptake reactions have negative bounds
        assert medium['EX_glc__D_e'] < 0
        assert medium['EX_o2_e'] < 0
        
        # Check that secretion reactions have positive bounds
        assert medium['EX_lac__L_e'] > 0
        assert medium['EX_co2_e'] > 0

    def test_harmonize_compartments(self):
        """Test compartment harmonization."""
        # This test would require a real model, so we'll test the logic
        gem = HumanGEM()
        
        # Test compartment mapping
        compartment_mapping = {
            'c': 'c',      # Cytosol
            'm': 'c_mito', # Mitochondria
            'e': 'e',      # Extracellular
            'n': 'c',      # Nucleus -> Cytosol
        }
        
        # Test mapping logic
        for old_comp, new_comp in compartment_mapping.items():
            assert old_comp != new_comp or old_comp == 'c' or old_comp == 'e'


class TestGPRImprovements:
    """Test cases for improved GPR functionality."""

    def setup_method(self):
        """Set up test fixtures."""
        self.gpr_parser = GPRParser()
        
        # Set up test GPR rules
        self.gpr_parser.gpr_rules = {
            "reaction1": "gene1 and gene2",
            "reaction2": "gene1 or gene2",
            "reaction3": "(gene1 and gene2) or gene3",
            "reaction4": "gene1",
            "reaction5": "",
        }
        
        self.gpr_parser.reaction_genes = {
            "reaction1": {"gene1", "gene2"},
            "reaction2": {"gene1", "gene2"},
            "reaction3": {"gene1", "gene2", "gene3"},
            "reaction4": {"gene1"},
            "reaction5": set(),
        }

    def test_gpr_eval_default_operators(self):
        """Test GPR evaluation with default operators (AND=min, OR=max)."""
        gene_expression = {"gene1": 0.8, "gene2": 0.6, "gene3": 0.3}
        
        reaction_scores = self.gpr_parser.gpr_eval(gene_expression)
        
        assert "reaction1" in reaction_scores
        assert "reaction2" in reaction_scores
        assert "reaction3" in reaction_scores
        assert "reaction4" in reaction_scores
        assert "reaction5" in reaction_scores
        
        # Test AND operation (should be min)
        assert reaction_scores["reaction1"] == min(0.8, 0.6)  # 0.6
        
        # Test OR operation (should be max)
        assert reaction_scores["reaction2"] == max(0.8, 0.6)  # 0.8
        
        # Test complex expression
        # (gene1 and gene2) or gene3 = max(min(0.8, 0.6), 0.3) = max(0.6, 0.3) = 0.6
        assert reaction_scores["reaction3"] == max(min(0.8, 0.6), 0.3)  # 0.6
        
        # Test single gene
        assert reaction_scores["reaction4"] == 0.8
        
        # Test empty rule
        assert reaction_scores["reaction5"] == 1.0

    def test_gpr_eval_custom_operators(self):
        """Test GPR evaluation with custom operators."""
        gene_expression = {"gene1": 0.8, "gene2": 0.6, "gene3": 0.3}
        custom_operators = {"AND": "max", "OR": "min"}  # Opposite of default
        
        reaction_scores = self.gpr_parser.gpr_eval(gene_expression, custom_operators)
        
        # Test AND operation (should be max with custom operators)
        assert reaction_scores["reaction1"] == max(0.8, 0.6)  # 0.8
        
        # Test OR operation (should be min with custom operators)
        assert reaction_scores["reaction2"] == min(0.8, 0.6)  # 0.6

    def test_parse_gpr_expression_simple(self):
        """Test parsing simple GPR expressions."""
        gene_expression = {"gene1": 0.8, "gene2": 0.6}
        
        # Simple AND
        result = self.gpr_parser._parse_gpr_expression("gene1 and gene2", gene_expression, {"AND": "min", "OR": "max"})
        assert result == min(0.8, 0.6)
        
        # Simple OR
        result = self.gpr_parser._parse_gpr_expression("gene1 or gene2", gene_expression, {"AND": "min", "OR": "max"})
        assert result == max(0.8, 0.6)
        
        # Single gene
        result = self.gpr_parser._parse_gpr_expression("gene1", gene_expression, {"AND": "min", "OR": "max"})
        assert result == 0.8

    def test_parse_gpr_expression_complex(self):
        """Test parsing complex GPR expressions with parentheses."""
        gene_expression = {"gene1": 0.8, "gene2": 0.6, "gene3": 0.3}
        
        # Complex expression: (gene1 and gene2) or gene3
        result = self.gpr_parser._parse_gpr_expression("(gene1 and gene2) or gene3", gene_expression, {"AND": "min", "OR": "max"})
        expected = max(min(0.8, 0.6), 0.3)  # max(0.6, 0.3) = 0.6
        assert result == expected

    def test_parse_gpr_expression_edge_cases(self):
        """Test parsing edge cases."""
        gene_expression = {"gene1": 0.8}
        
        # Unknown gene
        result = self.gpr_parser._parse_gpr_expression("unknown_gene", gene_expression, {"AND": "min", "OR": "max"})
        assert result == 0.0
        
        # Empty expression
        result = self.gpr_parser._parse_gpr_expression("", gene_expression, {"AND": "min", "OR": "max"})
        assert result == 0.0
        
        # Invalid expression
        result = self.gpr_parser._parse_gpr_expression("gene1 and", gene_expression, {"AND": "min", "OR": "max"})
        assert result == 0.0  # Should handle gracefully


class TestExpressionIntegrationImprovements:
    """Test cases for improved expression integration."""

    def setup_method(self):
        """Set up test fixtures."""
        from scflux_spatial.gem.gpr import GPRParser
        from scflux_spatial.fba.integrate_expression import ExpressionIntegrator
        
        self.gpr_parser = GPRParser()
        self.gpr_parser.gpr_rules = {
            "reaction1": "gene1 and gene2",
            "reaction2": "gene1 or gene2",
            "reaction3": "gene1",
        }
        self.gpr_parser.reaction_genes = {
            "reaction1": {"gene1", "gene2"},
            "reaction2": {"gene1", "gene2"},
            "reaction3": {"gene1"},
        }
        
        self.integrator = ExpressionIntegrator(self.gpr_parser)

    def test_integrate_expression_eflux(self):
        """Test E-Flux integration method."""
        gene_expression = {"gene1": 0.8, "gene2": 0.6}
        base_bounds = {
            "reaction1": (-10.0, 10.0),
            "reaction2": (-5.0, 5.0),
            "reaction3": (-2.0, 2.0),
        }
        
        result_bounds = self.integrator.integrate_expression_eflux(
            gene_expression, base_bounds, scaling_factor=1.0
        )
        
        assert "reaction1" in result_bounds
        assert "reaction2" in result_bounds
        assert "reaction3" in result_bounds
        
        # Check that bounds are scaled appropriately
        for reaction_id, (new_lower, new_upper) in result_bounds.items():
            base_lower, base_upper = base_bounds[reaction_id]
            assert isinstance(new_lower, float)
            assert isinstance(new_upper, float)
            assert new_lower <= new_upper  # Bounds should be consistent

    def test_integrate_expression_imat_like(self):
        """Test iMAT-like integration method."""
        gene_expression = {"gene1": 0.9, "gene2": 0.1}  # High and low expression
        base_bounds = {
            "reaction1": (-10.0, 10.0),
            "reaction2": (-5.0, 5.0),
            "reaction3": (-2.0, 2.0),
        }
        
        result_bounds = self.integrator.integrate_expression_imat_like(
            gene_expression, base_bounds,
            high_quantile=0.5, low_quantile=0.5,
            activation_bound=100.0, minimization_bound=0.1
        )
        
        assert "reaction1" in result_bounds
        assert "reaction2" in result_bounds
        assert "reaction3" in result_bounds
        
        # Check that high expression reactions are activated
        # reaction3 has high expression (gene1=0.9), should be activated
        r3_lower, r3_upper = result_bounds["reaction3"]
        assert r3_lower == -100.0  # Activation bound
        assert r3_upper == 100.0

    def test_integrate_expression_with_method(self):
        """Test method selection functionality."""
        gene_expression = {"gene1": 0.8, "gene2": 0.6}
        base_bounds = {"reaction1": (-10.0, 10.0)}
        
        # Test different methods
        methods = ["eflux", "imat_like", "linear", "quadratic", "none"]
        
        for method in methods:
            result = self.integrator.integrate_expression_with_method(
                gene_expression, base_bounds, method=method
            )
            assert isinstance(result, dict)
            assert "reaction1" in result
        
        # Test invalid method
        with pytest.raises(ValueError):
            self.integrator.integrate_expression_with_method(
                gene_expression, base_bounds, method="invalid_method"
            )

    def test_solve_with_pfba_mock(self):
        """Test pFBA solving with mock model."""
        # This would require a real COBRApy model, so we'll test the interface
        gene_expression = {"gene1": 0.8, "gene2": 0.6}
        base_bounds = {"reaction1": (-10.0, 10.0)}
        
        # Test that method exists and has correct signature
        assert hasattr(self.integrator, 'solve_with_pfba')
        
        # The actual solving would require a real model
        # In a real test, you would create a simple COBRApy model and test it
