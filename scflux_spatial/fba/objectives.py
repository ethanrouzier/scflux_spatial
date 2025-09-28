"""
Objective functions for flux balance analysis.

This module provides various objective functions including ATP maintenance,
biomass proxy, and custom objectives for spatial flux analysis.
"""

from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from rich.console import Console

console = Console()


class ObjectiveFunction:
    """Objective function manager for FBA."""

    def __init__(self, model):
        """
        Initialize objective function manager.

        Args:
            model: COBRApy model object
        """
        self.model = model
        self.objective_functions = {}

    def set_biomass_objective(self, biomass_reaction: str = "biomass_human") -> None:
        """
        Set biomass as objective function.

        Args:
            biomass_reaction: ID of the biomass reaction
        """
        if biomass_reaction in self.model.reactions:
            self.model.objective = biomass_reaction
            self.objective_functions["biomass"] = biomass_reaction
            console.print(f"Set biomass objective: {biomass_reaction}")
        else:
            console.print(f"Warning: Biomass reaction {biomass_reaction} not found")

    def set_atp_maintenance_objective(self, atp_reaction: str = "ATPM") -> None:
        """
        Set ATP maintenance as objective function.

        Args:
            atp_reaction: ID of the ATP maintenance reaction
        """
        if atp_reaction in self.model.reactions:
            self.model.objective = atp_reaction
            self.objective_functions["atp_maintenance"] = atp_reaction
            console.print(f"Set ATP maintenance objective: {atp_reaction}")
        else:
            console.print(f"Warning: ATP maintenance reaction {atp_reaction} not found")

    def set_glucose_uptake_objective(self, glucose_reaction: str = "EX_glc__D_e") -> None:
        """
        Set glucose uptake as objective function.

        Args:
            glucose_reaction: ID of the glucose exchange reaction
        """
        if glucose_reaction in self.model.reactions:
            self.model.objective = glucose_reaction
            self.objective_functions["glucose_uptake"] = glucose_reaction
            console.print(f"Set glucose uptake objective: {glucose_reaction}")
        else:
            console.print(f"Warning: Glucose reaction {glucose_reaction} not found")

    def set_custom_objective(
        self, 
        objective_dict: Dict[str, float],
        direction: str = "max"
    ) -> None:
        """
        Set custom objective function.

        Args:
            objective_dict: Dictionary mapping reaction IDs to coefficients
            direction: Optimization direction ('max' or 'min')
        """
        # Clear existing objective
        self.model.objective = 0
        
        # Add reactions to objective
        for reaction_id, coefficient in objective_dict.items():
            if reaction_id in self.model.reactions:
                reaction = self.model.reactions.get_by_id(reaction_id)
                self.model.objective += coefficient * reaction
            else:
                console.print(f"Warning: Reaction {reaction_id} not found in model")

        # Store objective information
        self.objective_functions["custom"] = {
            "reactions": objective_dict,
            "direction": direction
        }
        
        console.print(f"Set custom objective with {len(objective_dict)} reactions")

    def set_multi_objective(
        self, 
        objectives: List[Dict[str, float]],
        weights: Optional[List[float]] = None
    ) -> None:
        """
        Set multi-objective function.

        Args:
            objectives: List of objective dictionaries
            weights: Weights for each objective (if None, equal weights)
        """
        if weights is None:
            weights = [1.0 / len(objectives)] * len(objectives)
        
        if len(objectives) != len(weights):
            raise ValueError("Number of objectives and weights must match")

        # Clear existing objective
        self.model.objective = 0
        
        # Combine objectives
        combined_objective = {}
        for obj_dict, weight in zip(objectives, weights):
            for reaction_id, coefficient in obj_dict.items():
                if reaction_id in combined_objective:
                    combined_objective[reaction_id] += weight * coefficient
                else:
                    combined_objective[reaction_id] = weight * coefficient

        # Set combined objective
        for reaction_id, coefficient in combined_objective.items():
            if reaction_id in self.model.reactions:
                reaction = self.model.reactions.get_by_id(reaction_id)
                self.model.objective += coefficient * reaction

        self.objective_functions["multi_objective"] = {
            "objectives": objectives,
            "weights": weights
        }
        
        console.print(f"Set multi-objective with {len(objectives)} objectives")

    def create_biomass_proxy(
        self, 
        essential_reactions: List[str],
        weights: Optional[List[float]] = None
    ) -> Dict[str, float]:
        """
        Create a biomass proxy objective from essential reactions.

        Args:
            essential_reactions: List of essential reaction IDs
            weights: Weights for each reaction (if None, equal weights)

        Returns:
            Dictionary mapping reaction IDs to coefficients
        """
        if weights is None:
            weights = [1.0 / len(essential_reactions)] * len(essential_reactions)
        
        if len(essential_reactions) != len(weights):
            raise ValueError("Number of reactions and weights must match")

        biomass_proxy = {}
        for reaction_id, weight in zip(essential_reactions, weights):
            if reaction_id in self.model.reactions:
                biomass_proxy[reaction_id] = weight

        console.print(f"Created biomass proxy with {len(biomass_proxy)} reactions")
        return biomass_proxy

    def create_atp_production_objective(
        self, 
        atp_reactions: List[str],
        weights: Optional[List[float]] = None
    ) -> Dict[str, float]:
        """
        Create ATP production objective.

        Args:
            atp_reactions: List of ATP-producing reaction IDs
            weights: Weights for each reaction (if None, equal weights)

        Returns:
            Dictionary mapping reaction IDs to coefficients
        """
        if weights is None:
            weights = [1.0 / len(atp_reactions)] * len(atp_reactions)
        
        if len(atp_reactions) != len(weights):
            raise ValueError("Number of reactions and weights must match")

        atp_objective = {}
        for reaction_id, weight in zip(atp_reactions, weights):
            if reaction_id in self.model.reactions:
                atp_objective[reaction_id] = weight

        console.print(f"Created ATP production objective with {len(atp_objective)} reactions")
        return atp_objective

    def get_objective_value(self) -> float:
        """Get current objective value."""
        return self.model.objective.expression

    def get_objective_reactions(self) -> List[str]:
        """Get reactions involved in current objective."""
        objective_reactions = []
        for reaction in self.model.reactions:
            if reaction.objective_coefficient != 0:
                objective_reactions.append(reaction.id)
        return objective_reactions

    def analyze_objective_sensitivity(
        self, 
        perturbation_range: float = 0.1
    ) -> Dict[str, float]:
        """
        Analyze sensitivity of objective to reaction perturbations.

        Args:
            perturbation_range: Range of perturbation (±fraction)

        Returns:
            Dictionary mapping reaction IDs to sensitivity scores
        """
        console.print("Analyzing objective sensitivity...")

        # Get current objective value
        current_objective = self.model.optimize().objective_value
        
        sensitivities = {}
        
        for reaction in self.model.reactions:
            if reaction.objective_coefficient != 0:
                # Store original bounds
                original_lower = reaction.lower_bound
                original_upper = reaction.upper_bound
                
                # Perturb bounds
                perturbation = abs(reaction.upper_bound - reaction.lower_bound) * perturbation_range
                reaction.lower_bound = original_lower - perturbation
                reaction.upper_bound = original_upper + perturbation
                
                # Solve and get new objective
                new_solution = self.model.optimize()
                if new_solution.status == "optimal":
                    new_objective = new_solution.objective_value
                    sensitivity = abs(new_objective - current_objective) / current_objective
                    sensitivities[reaction.id] = sensitivity
                
                # Restore original bounds
                reaction.lower_bound = original_lower
                reaction.upper_bound = original_upper

        console.print(f"Analyzed sensitivity for {len(sensitivities)} reactions")
        return sensitivities

    def export_objective_info(self, output_path: str) -> None:
        """
        Export objective function information.

        Args:
            output_path: Path to save the objective information
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Export objective reactions
        objective_data = []
        for reaction in self.model.reactions:
            if reaction.objective_coefficient != 0:
                objective_data.append({
                    "reaction_id": reaction.id,
                    "reaction_name": reaction.name,
                    "objective_coefficient": reaction.objective_coefficient,
                    "lower_bound": reaction.lower_bound,
                    "upper_bound": reaction.upper_bound
                })

        objective_df = pd.DataFrame(objective_data)
        objective_df.to_csv(output_path / "objective_reactions.csv", index=False)

        # Export objective functions summary
        functions_df = pd.DataFrame([
            {"function_type": func_type, "details": str(details)}
            for func_type, details in self.objective_functions.items()
        ])
        functions_df.to_csv(output_path / "objective_functions.csv", index=False)

        console.print(f"Objective information exported to {output_path}")

    def get_available_objectives(self) -> List[str]:
        """Get list of available objective reactions in the model."""
        available = []
        
        # Common objective reactions
        common_objectives = [
            "biomass_human",
            "ATPM",
            "EX_glc__D_e",
            "EX_o2_e",
            "EX_lac__L_e"
        ]
        
        for obj_reaction in common_objectives:
            if obj_reaction in self.model.reactions:
                available.append(obj_reaction)
        
        return available
