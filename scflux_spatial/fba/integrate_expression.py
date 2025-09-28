"""
Minimal E-Flux integration for gene expression data.

This module provides basic E-Flux functionality for integrating
gene expression data into metabolic models.
"""

from typing import Dict, Iterable
import numpy as np
from cobra import Model


def eflux_apply(model: Model, rxn_to_gene_activity: Dict[str, float],
                lb_scale: float = 1.0, ub_scale: float = 1.0) -> None:
    """
    E-Flux minimal: scale reaction bounds proportionally to gene activity.
    
    Args:
        model: COBRA model
        rxn_to_gene_activity: dict reaction_id -> activity in [0,1]
        lb_scale: Scaling factor for lower bounds
        ub_scale: Scaling factor for upper bounds
    """
    for rxn in model.reactions:
        a = float(rxn_to_gene_activity.get(rxn.id, 1.0))
        if rxn.lower_bound < 0:
            rxn.lower_bound = rxn.lower_bound * (lb_scale * a)
        if rxn.upper_bound > 0:
            rxn.upper_bound = rxn.upper_bound * (ub_scale * a)


def normalize_activities(vals: Iterable[float]) -> Iterable[float]:
    """
    Normalize gene activities to [0,1] range using percentile scaling.
    
    Args:
        vals: Iterable of gene expression values
        
    Returns:
        Normalized activities in [0,1] range
    """
    v = np.asarray(list(vals), dtype=float)
    v = np.nan_to_num(v, nan=0.0, posinf=0.0, neginf=0.0)
    if v.size == 0:
        return v
    lo, hi = np.percentile(v, [5, 95])
    if hi <= lo:
        return np.clip((v - lo), 0, 1)
    return np.clip((v - lo) / (hi - lo), 0, 1)


def imat_apply(model: Model, rxn_to_gene_activity: Dict[str, float]) -> None:
    """
    iMAT integration (placeholder for future implementation).
    
    Args:
        model: COBRA model
        rxn_to_gene_activity: dict reaction_id -> activity in [0,1]
    """
    raise NotImplementedError("iMAT integration not yet implemented. Use eflux_apply for now.")


def linear_apply(model: Model, rxn_to_gene_activity: Dict[str, float]) -> None:
    """
    Linear integration method (placeholder for future implementation).
    
    Args:
        model: COBRA model
        rxn_to_gene_activity: dict reaction_id -> activity in [0,1]
    """
    raise NotImplementedError("Linear integration not yet implemented. Use eflux_apply for now.")


def solve_with_pfba(model: Model, objective: str = "biomass_reaction") -> Dict[str, float]:
    """
    Solve model using parsimonious FBA.
    
    Args:
        model: COBRA model
        objective: Objective reaction ID
        
    Returns:
        Dictionary of reaction fluxes
    """
    try:
        from cobra.flux_analysis import pfba
        solution = pfba(model)
        return {rxn.id: solution.fluxes[rxn.id] for rxn in model.reactions}
    except ImportError:
        raise NotImplementedError("COBRApy not available for FBA solving")


def integrate_expression_with_method(gene_expression: Dict[str, float], 
                                   method: str = "eflux") -> Dict[str, float]:
    """
    Integrate gene expression using specified method.
    
    Args:
        gene_expression: Gene expression values
        method: Integration method ("eflux", "imat", "linear")
        
    Returns:
        Dictionary of reaction activities
    """
    if method == "eflux":
        # Simple E-Flux implementation
        return {f"RXN_{gene}": expr for gene, expr in gene_expression.items()}
    elif method == "imat":
        raise NotImplementedError("iMAT integration not yet implemented")
    elif method == "linear":
        raise NotImplementedError("Linear integration not yet implemented")
    else:
        raise ValueError(f"Unknown method: {method}")


class ExpressionIntegrator:
    """
    Expression integrator for metabolic models.
    """
    
    def __init__(self, method: str = "eflux"):
        self.method = method
    
    def integrate(self, gene_expression: Dict[str, float]) -> Dict[str, float]:
        """Integrate gene expression data."""
        return integrate_expression_with_method(gene_expression, self.method)