#!/usr/bin/env python3
"""
Unit tests for GPR (Gene-Protein-Reaction) parsing and evaluation.
"""

import unittest
import numpy as np
from scflux_spatial.gem.gpr import GPRParser


class TestGPRParser(unittest.TestCase):
    """Test cases for GPR parser functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.gpr_parser = GPRParser()
        self.gene_expression = {
            'GENE_A': 5.0,
            'GENE_B': 3.0,
            'GENE_C': 2.0,
            'GENE_D': 4.0,
            'GENE_E': 1.0,
            'GENE_F': 6.0
        }
    
    def test_simple_and_operator(self):
        """Test simple AND operation."""
        result = self.gpr_parser._evaluate_gpr_rule_with_operators(
            'GENE_A and GENE_B', 
            {'GENE_A': 5.0, 'GENE_B': 3.0}
        )
        # AND should return minimum value
        expected = min(5.0, 3.0)  # = 3.0
        self.assertEqual(result, 3.0)
    
    def test_simple_or_operator(self):
        """Test simple OR operation."""
        result = self.gpr_parser._evaluate_gpr_rule_with_operators(
            'GENE_A or GENE_B', 
            {'GENE_A': 5.0, 'GENE_B': 3.0}
        )
        # OR should return maximum value
        expected = max(5.0, 3.0)  # = 5.0
        self.assertEqual(result, 5.0)
    
    def test_complex_nested_expressions(self):
        """Test complex nested logical expressions."""
        result = self.gpr_parser._evaluate_gpr_rule_with_operators(
            '(GENE_A and GENE_B) or (GENE_C and GENE_D)',
            self.gene_expression
        )
        # (5.0 and 3.0) or (2.0 and 4.0) = 3.0 or 2.0 = 3.0
        expected = max(min(5.0, 3.0), min(2.0, 4.0))
        self.assertEqual(result, expected)
    
    def test_edge_cases(self):
        """Test edge cases and error handling."""
        # Test with zero values
        zero_expression = {'GENE_A': 0.0, 'GENE_B': 3.0}
        result = self.gpr_parser._evaluate_gpr_rule_with_operators(
            'GENE_A and GENE_B',
            zero_expression
        )
        self.assertEqual(result, 0.0)
        
        # Test with negative values
        negative_expression = {'GENE_A': -2.0, 'GENE_B': 3.0}
        result = self.gpr_parser._evaluate_gpr_rule_with_operators(
            'GENE_A or GENE_B',
            negative_expression
        )
        self.assertEqual(result, 3.0)
    
    def test_missing_genes(self):
        """Test handling of missing genes."""
        with self.assertRaises(KeyError):
            self.gpr_parser._evaluate_gpr_rule_with_operators(
                'GENE_X and GENE_Y',
                self.gene_expression
            )
    
    def test_empty_expression(self):
        """Test empty expression handling."""
        with self.assertRaises(ValueError):
            self.gpr_parser._evaluate_gpr_rule_with_operators(
                '',
                self.gene_expression
            )


if __name__ == '__main__':
    unittest.main(verbosity=2)
