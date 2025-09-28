"""
Gene-Protein-Reaction (GPR) parsing and mapping.

This module provides functionality to parse GPR rules and map genes to reactions
for constraint-based modeling.
"""

import re
from typing import Dict, List, Set, Tuple, Union, Optional

import numpy as np
import pandas as pd
from rich.console import Console

console = Console()


class GPRParser:
    """Parser for Gene-Protein-Reaction (GPR) rules."""

    def __init__(self):
        """Initialize GPR parser."""
        self.gpr_rules: Dict[str, str] = {}
        self.gene_reactions: Dict[str, Set[str]] = {}
        self.reaction_genes: Dict[str, Set[str]] = {}

    def parse_gpr_rule(self, gpr_string: str) -> Tuple[Set[str], str]:
        """
        Parse a GPR rule string.

        Args:
            gpr_string: GPR rule string (e.g., "gene1 and gene2 or gene3")

        Returns:
            Tuple of (gene_set, parsed_expression)
        """
        if not gpr_string or gpr_string == "":
            return set(), ""

        # Clean the GPR string
        gpr_string = gpr_string.strip()
        
        # Extract gene IDs using regex
        gene_pattern = r'\b[A-Za-z0-9_]+\b'
        genes = set(re.findall(gene_pattern, gpr_string))
        
        # Remove common keywords that are not genes
        keywords = {'and', 'or', 'not', '(', ')'}
        genes = genes - keywords
        
        return genes, gpr_string

    def load_gpr_rules(self, model_reactions: List) -> None:
        """
        Load GPR rules from model reactions.

        Args:
            model_reactions: List of COBRApy reaction objects
        """
        console.print("Loading GPR rules from model reactions...")

        for reaction in model_reactions:
            reaction_id = reaction.id
            gpr_string = reaction.gene_reaction_rule
            
            # Parse GPR rule
            genes, parsed_expression = self.parse_gpr_rule(gpr_string)
            
            # Store GPR rule
            self.gpr_rules[reaction_id] = parsed_expression
            
            # Store gene-reaction mappings
            self.reaction_genes[reaction_id] = genes
            for gene in genes:
                if gene not in self.gene_reactions:
                    self.gene_reactions[gene] = set()
                self.gene_reactions[gene].add(reaction_id)

        console.print(f"Loaded GPR rules for {len(self.gpr_rules)} reactions")

    def gpr_eval(self, expression: Dict[str, float], operators: Optional[Dict[str, str]] = None) -> Dict[str, float]:
        """
        Evaluate GPR rules for all reactions using custom operators.
        
        Args:
            expression: Dictionary mapping gene IDs to expression values
            operators: Custom operators (default: AND=min, OR=max)
            
        Returns:
            Dictionary mapping reaction IDs to GPR scores
        """
        if operators is None:
            operators = {'AND': 'min', 'OR': 'max'}
        
        reaction_scores = {}
        
        for reaction_id in self.gpr_rules:
            score = self._evaluate_gpr_rule_with_operators(
                reaction_id, expression, operators
            )
            reaction_scores[reaction_id] = score
            
        return reaction_scores

    def _evaluate_gpr_rule_with_operators(
        self, 
        reaction_id: str, 
        gene_expression: Dict[str, float],
        operators: Dict[str, str]
    ) -> float:
        """
        Evaluate GPR rule using custom operators (AND=min, OR=max).
        
        Args:
            reaction_id: ID of the reaction
            gene_expression: Dictionary mapping gene IDs to expression values
            operators: Custom operators mapping
            
        Returns:
            Evaluated GPR value (0-1)
        """
        if reaction_id not in self.gpr_rules:
            return 1.0  # Default to full activity if no GPR rule

        gpr_rule = self.gpr_rules[reaction_id]
        
        if not gpr_rule:
            return 1.0

        # Parse and evaluate using custom operators
        return self._parse_gpr_expression(gpr_rule, gene_expression, operators)

    def _parse_gpr_expression(
        self, 
        gpr_rule: str, 
        gene_expression: Dict[str, float],
        operators: Dict[str, str]
    ) -> float:
        """
        Parse and evaluate GPR expression using custom operators.
        
        Args:
            gpr_rule: GPR rule string
            gene_expression: Gene expression dictionary
            operators: Custom operators mapping
            
        Returns:
            Evaluated score
        """
        try:
            # Handle parentheses by recursive evaluation
            while '(' in gpr_rule and ')' in gpr_rule:
                # Find innermost parentheses
                start = gpr_rule.rfind('(')
                end = gpr_rule.find(')', start)
                
                if start == -1 or end == -1:
                    break
                
                # Extract and evaluate inner expression
                inner_expr = gpr_rule[start+1:end]
                inner_result = self._parse_gpr_expression(inner_expr, gene_expression, operators)
                
                # Replace with result
                gpr_rule = gpr_rule[:start] + str(inner_result) + gpr_rule[end+1:]
            
            # Split by OR operators first (lower precedence)
            or_parts = gpr_rule.split(' or ')
            if len(or_parts) > 1:
                results = [self._parse_gpr_expression(part.strip(), gene_expression, operators) 
                          for part in or_parts]
                return max(results) if operators.get('OR') == 'max' else min(results)
            
            # Split by AND operators (higher precedence)
            and_parts = gpr_rule.split(' and ')
            if len(and_parts) > 1:
                results = [self._parse_gpr_expression(part.strip(), gene_expression, operators) 
                          for part in and_parts]
                return min(results) if operators.get('AND') == 'min' else max(results)
            
            # Single gene or number
            gpr_rule = gpr_rule.strip()
            
            # Check if it's a gene
            if gpr_rule in gene_expression:
                return gene_expression[gpr_rule]
            
            # Check if it's a number
            try:
                return float(gpr_rule)
            except ValueError:
                return 0.0  # Unknown gene or invalid expression
                
        except Exception as e:
            console.print(f"Error parsing GPR expression '{gpr_rule}': {e}")
            return 0.0

    def evaluate_gpr_rule(self, reaction_id: str, gene_expression: Dict[str, float]) -> float:
        """
        Evaluate GPR rule for a reaction given gene expression (legacy method).

        Args:
            reaction_id: ID of the reaction
            gene_expression: Dictionary mapping gene IDs to expression values

        Returns:
            Evaluated GPR value (0-1)
        """
        if reaction_id not in self.gpr_rules:
            return 1.0  # Default to full activity if no GPR rule

        gpr_rule = self.gpr_rules[reaction_id]
        
        if not gpr_rule:
            return 1.0

        # Replace gene IDs with expression values
        expression_rule = gpr_rule
        for gene in self.reaction_genes[reaction_id]:
            if gene in gene_expression:
                expr_value = gene_expression[gene]
            else:
                expr_value = 0.0  # Default to no expression if gene not found
            
            # Replace gene ID with expression value
            expression_rule = re.sub(r'\b' + re.escape(gene) + r'\b', str(expr_value), expression_rule)

        # Replace logical operators with Python equivalents
        expression_rule = expression_rule.replace('and', '&')
        expression_rule = expression_rule.replace('or', '|')
        expression_rule = expression_rule.replace('not', '~')

        try:
            # Evaluate the expression
            result = eval(expression_rule)
            
            # Convert boolean result to float
            if isinstance(result, bool):
                return float(result)
            elif isinstance(result, (int, float)):
                return float(result)
            else:
                return 1.0
                
        except Exception as e:
            console.print(f"Error evaluating GPR rule for {reaction_id}: {e}")
            return 1.0

    def get_reaction_bounds_from_expression(
        self, 
        reaction_id: str, 
        gene_expression: Dict[str, float],
        method: str = "imat"
    ) -> Tuple[float, float]:
        """
        Calculate reaction bounds based on gene expression.

        Args:
            reaction_id: ID of the reaction
            gene_expression: Dictionary mapping gene IDs to expression values
            method: Method for calculating bounds ('imat', 'eflux', 'linear')

        Returns:
            Tuple of (lower_bound, upper_bound)
        """
        gpr_value = self.evaluate_gpr_rule(reaction_id, gene_expression)
        
        if method == "imat":
            # iMAT-like: binary bounds based on threshold
            threshold = 0.5
            if gpr_value > threshold:
                return (-1000.0, 1000.0)  # Active
            else:
                return (0.0, 0.0)  # Inactive
                
        elif method == "eflux":
            # E-Flux: proportional bounds
            max_flux = 100.0  # Maximum flux value
            bound = gpr_value * max_flux
            return (-bound, bound)
            
        elif method == "linear":
            # Linear scaling
            max_flux = 100.0
            bound = gpr_value * max_flux
            return (-bound, bound)
            
        else:
            raise ValueError(f"Unknown method: {method}")

    def map_expression_to_reactions(
        self, 
        gene_expression: Dict[str, float],
        method: str = "imat"
    ) -> Dict[str, Tuple[float, float]]:
        """
        Map gene expression to reaction bounds for all reactions.

        Args:
            gene_expression: Dictionary mapping gene IDs to expression values
            method: Method for calculating bounds

        Returns:
            Dictionary mapping reaction IDs to (lower_bound, upper_bound) tuples
        """
        reaction_bounds = {}
        
        for reaction_id in self.reaction_genes:
            bounds = self.get_reaction_bounds_from_expression(
                reaction_id, gene_expression, method
            )
            reaction_bounds[reaction_id] = bounds
            
        return reaction_bounds

    def get_essential_genes(self, threshold: float = 0.1) -> Set[str]:
        """
        Get genes that are essential for multiple reactions.

        Args:
            threshold: Minimum number of reactions a gene must be involved in

        Returns:
            Set of essential gene IDs
        """
        essential_genes = set()
        
        for gene, reactions in self.gene_reactions.items():
            if len(reactions) >= threshold:
                essential_genes.add(gene)
                
        return essential_genes

    def get_reaction_complexity(self) -> Dict[str, int]:
        """
        Get complexity (number of genes) for each reaction.

        Returns:
            Dictionary mapping reaction IDs to number of genes
        """
        complexity = {}
        
        for reaction_id, genes in self.reaction_genes.items():
            complexity[reaction_id] = len(genes)
            
        return complexity

    def export_gpr_network(self, output_path: str) -> None:
        """
        Export GPR network to CSV files.

        Args:
            output_path: Path to save the GPR network files
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        console.print(f"Exporting GPR network to {output_path}")

        # Export reaction-gene mappings
        reaction_gene_data = []
        for reaction_id, genes in self.reaction_genes.items():
            for gene in genes:
                reaction_gene_data.append({
                    "reaction_id": reaction_id,
                    "gene_id": gene,
                    "gpr_rule": self.gpr_rules.get(reaction_id, "")
                })

        reaction_gene_df = pd.DataFrame(reaction_gene_data)
        reaction_gene_df.to_csv(output_path / "reaction_gene_mappings.csv", index=False)

        # Export gene-reaction mappings
        gene_reaction_data = []
        for gene, reactions in self.gene_reactions.items():
            for reaction in reactions:
                gene_reaction_data.append({
                    "gene_id": gene,
                    "reaction_id": reaction,
                    "reaction_count": len(self.gene_reactions[gene])
                })

        gene_reaction_df = pd.DataFrame(gene_reaction_data)
        gene_reaction_df.to_csv(output_path / "gene_reaction_mappings.csv", index=False)

        console.print("GPR network exported successfully")

    def get_gene_coverage(self, gene_expression: Dict[str, float]) -> Dict[str, float]:
        """
        Get coverage of genes in the expression data.

        Args:
            gene_expression: Dictionary mapping gene IDs to expression values

        Returns:
            Dictionary with coverage statistics
        """
        total_genes = len(self.gene_reactions)
        covered_genes = len(set(gene_expression.keys()) & set(self.gene_reactions.keys()))
        
        coverage_stats = {
            "total_genes": total_genes,
            "covered_genes": covered_genes,
            "coverage_fraction": covered_genes / total_genes if total_genes > 0 else 0.0,
            "missing_genes": total_genes - covered_genes
        }
        
        return coverage_stats
