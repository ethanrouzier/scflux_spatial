"""
Flux balance analysis module.

This module provides tools for integrating gene expression data
into metabolic models.
"""

from .integrate_expression import (
    eflux_apply, normalize_activities, solve_with_pfba, 
    integrate_expression_with_method, ExpressionIntegrator
)

__all__ = [
    "eflux_apply",
    "normalize_activities",
    "solve_with_pfba",
    "integrate_expression_with_method",
    "ExpressionIntegrator"
]