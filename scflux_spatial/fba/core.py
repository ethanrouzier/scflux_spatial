"""
Flux Balance Analysis (FBA) core functionality.

This module provides the core FBA solver with support for FBA, pFBA, and FVA.
"""

from typing import Dict, List, Optional, Tuple, Union

import cobra
import numpy as np
import pandas as pd
from rich.console import Console

console = Console()


class FBASolver:
    """Flux Balance Analysis solver with multiple methods."""

    def __init__(
        self,
        model: cobra.Model,
        solver: str = "glpk",
        tolerance: float = 1e-6,
        max_iterations: int = 1000,
    ):
        """
        Initialize FBA solver.

        Args:
            model: COBRApy model object
            solver: Solver to use ('glpk', 'cplex', 'gurobi')
            tolerance: Solver tolerance
            max_iterations: Maximum number of iterations
        """
        self.model = model
        self.solver = solver
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        
        # Set solver
        self.model.solver = solver
        
        # Store results
        self.last_solution = None
        self.last_objective = None

    def solve_fba(self, objective: Optional[str] = None) -> Dict[str, float]:
        """
        Solve standard FBA.

        Args:
            objective: Objective reaction ID (if None, uses model objective)

        Returns:
            Dictionary mapping reaction IDs to flux values
        """
        console.print("Solving FBA...")

        # Set objective if specified
        if objective:
            self._set_objective(objective)

        # Solve the model
        solution = self.model.optimize()
        
        if solution.status == "optimal":
            self.last_solution = solution
            self.last_objective = solution.objective_value
            
            # Convert to dictionary
            flux_dict = {rxn.id: solution.fluxes[rxn.id] for rxn in self.model.reactions}
            
            console.print(f"FBA solved successfully. Objective value: {solution.objective_value:.6f}")
            return flux_dict
        else:
            console.print(f"FBA failed with status: {solution.status}")
            return {}

    def solve_pfba(self, objective: Optional[str] = None) -> Dict[str, float]:
        """
        Solve parsimonious FBA (pFBA).

        Args:
            objective: Objective reaction ID (if None, uses model objective)

        Returns:
            Dictionary mapping reaction IDs to flux values
        """
        console.print("Solving parsimonious FBA...")

        # Set objective if specified
        if objective:
            self._set_objective(objective)

        # Solve pFBA
        solution = cobra.flux_analysis.pfba(self.model)
        
        if solution.status == "optimal":
            self.last_solution = solution
            self.last_objective = solution.objective_value
            
            # Convert to dictionary
            flux_dict = {rxn.id: solution.fluxes[rxn.id] for rxn in self.model.reactions}
            
            console.print(f"pFBA solved successfully. Objective value: {solution.objective_value:.6f}")
            return flux_dict
        else:
            console.print(f"pFBA failed with status: {solution.status}")
            return {}

    def solve_fva(
        self, 
        fraction_of_optimum: float = 0.9,
        processes: int = 1
    ) -> Tuple[Dict[str, float], Dict[str, float]]:
        """
        Solve Flux Variability Analysis (FVA).

        Args:
            fraction_of_optimum: Fraction of optimum objective value to maintain
            processes: Number of processes to use

        Returns:
            Tuple of (minimum_fluxes, maximum_fluxes) dictionaries
        """
        console.print("Solving Flux Variability Analysis...")

        # First solve FBA to get optimum
        self.solve_fba()
        
        if self.last_solution is None:
            console.print("FBA failed, cannot perform FVA")
            return {}, {}

        # Set objective constraint
        objective_value = self.last_objective * fraction_of_optimum
        
        # Perform FVA
        fva_results = cobra.flux_analysis.flux_variability_analysis(
            self.model,
            fraction_of_optimum=fraction_of_optimum,
            processes=processes
        )
        
        # Convert to dictionaries
        min_fluxes = fva_results['minimum'].to_dict()
        max_fluxes = fva_results['maximum'].to_dict()
        
        console.print(f"FVA completed for {len(min_fluxes)} reactions")
        return min_fluxes, max_fluxes

    def solve_moma(self, wild_type_fluxes: Dict[str, float]) -> Dict[str, float]:
        """
        Solve Minimization of Metabolic Adjustment (MOMA).

        Args:
            wild_type_fluxes: Wild-type flux distribution

        Returns:
            Dictionary mapping reaction IDs to flux values
        """
        console.print("Solving MOMA...")

        # Create wild-type solution object
        wt_solution = cobra.Solution(
            objective_value=0.0,
            status="optimal",
            fluxes=pd.Series(wild_type_fluxes)
        )

        # Solve MOMA
        moma_solution = cobra.flux_analysis.moma(self.model, wt_solution)
        
        if moma_solution.status == "optimal":
            self.last_solution = moma_solution
            self.last_objective = moma_solution.objective_value
            
            # Convert to dictionary
            flux_dict = {rxn.id: moma_solution.fluxes[rxn.id] for rxn in self.model.reactions}
            
            console.print(f"MOMA solved successfully. Distance: {moma_solution.objective_value:.6f}")
            return flux_dict
        else:
            console.print(f"MOMA failed with status: {moma_solution.status}")
            return {}

    def set_reaction_bounds(
        self, 
        reaction_bounds: Dict[str, Tuple[float, float]]
    ) -> None:
        """
        Set bounds for multiple reactions.

        Args:
            reaction_bounds: Dictionary mapping reaction IDs to (lower, upper) bounds
        """
        for reaction_id, (lower, upper) in reaction_bounds.items():
            if reaction_id in self.model.reactions:
                reaction = self.model.reactions.get_by_id(reaction_id)
                reaction.lower_bound = lower
                reaction.upper_bound = upper

    def get_reaction_bounds(self) -> Dict[str, Tuple[float, float]]:
        """
        Get current bounds for all reactions.

        Returns:
            Dictionary mapping reaction IDs to (lower, upper) bounds
        """
        bounds = {}
        for reaction in self.model.reactions:
            bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)
        return bounds

    def _set_objective(self, objective: str) -> None:
        """Set the objective function."""
        if objective in self.model.reactions:
            self.model.objective = objective
        else:
            console.print(f"Warning: Objective reaction {objective} not found in model")

    def get_objective_value(self) -> Optional[float]:
        """Get the last objective value."""
        return self.last_objective

    def get_solution_status(self) -> Optional[str]:
        """Get the last solution status."""
        if self.last_solution:
            return self.last_solution.status
        return None

    def get_active_reactions(
        self, 
        flux_threshold: float = 1e-6
    ) -> List[str]:
        """
        Get list of active reactions from last solution.

        Args:
            flux_threshold: Minimum absolute flux to consider active

        Returns:
            List of active reaction IDs
        """
        if self.last_solution is None:
            return []

        active_reactions = []
        for reaction in self.model.reactions:
            flux = self.last_solution.fluxes[reaction.id]
            if abs(flux) > flux_threshold:
                active_reactions.append(reaction.id)

        return active_reactions

    def get_flux_summary(self) -> Dict[str, Union[int, float]]:
        """
        Get summary statistics of the last solution.

        Returns:
            Dictionary with flux summary statistics
        """
        if self.last_solution is None:
            return {}

        fluxes = self.last_solution.fluxes
        active_fluxes = fluxes[abs(fluxes) > 1e-6]
        
        summary = {
            "total_reactions": len(self.model.reactions),
            "active_reactions": len(active_fluxes),
            "total_flux": fluxes.sum(),
            "max_flux": fluxes.max(),
            "min_flux": fluxes.min(),
            "mean_flux": active_fluxes.mean() if len(active_fluxes) > 0 else 0.0,
            "objective_value": self.last_objective
        }

        return summary

    def export_flux_results(self, output_path: str) -> None:
        """
        Export flux results to CSV.

        Args:
            output_path: Path to save the results
        """
        if self.last_solution is None:
            console.print("No solution to export")
            return

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Create results DataFrame
        results_data = []
        for reaction in self.model.reactions:
            results_data.append({
                "reaction_id": reaction.id,
                "reaction_name": reaction.name,
                "flux": self.last_solution.fluxes[reaction.id],
                "lower_bound": reaction.lower_bound,
                "upper_bound": reaction.upper_bound,
                "objective_coefficient": reaction.objective_coefficient
            })

        results_df = pd.DataFrame(results_data)
        results_df.to_csv(output_path, index=False)
        
        console.print(f"Flux results exported to {output_path}")

    def check_model_consistency(self) -> Dict[str, bool]:
        """
        Check model consistency and feasibility.

        Returns:
            Dictionary with consistency check results
        """
        checks = {}
        
        # Check if model can be solved
        try:
            solution = self.model.optimize()
            checks["solvable"] = solution.status == "optimal"
            checks["objective_value"] = solution.objective_value if solution.status == "optimal" else None
        except Exception as e:
            checks["solvable"] = False
            checks["error"] = str(e)

        # Check for blocked reactions
        try:
            blocked_reactions = cobra.flux_analysis.find_blocked_reactions(self.model)
            checks["has_blocked_reactions"] = len(blocked_reactions) > 0
            checks["num_blocked_reactions"] = len(blocked_reactions)
        except Exception as e:
            checks["blocked_reactions_check"] = False
            checks["blocked_reactions_error"] = str(e)

        # Check for unused metabolites
        unused_metabolites = []
        for metabolite in self.model.metabolites:
            if len(metabolite.reactions) == 0:
                unused_metabolites.append(metabolite.id)
        
        checks["has_unused_metabolites"] = len(unused_metabolites) > 0
        checks["num_unused_metabolites"] = len(unused_metabolites)

        return checks
