"""
Kinetics models for cellular metabolism.

This module provides various kinetic models including Michaelis-Menten
and consumption kinetics for spatial flux analysis.

Units:
- Concentrations: mol·L⁻¹
- Time: seconds (s)
- Rates: mol·L⁻¹·s⁻¹
- Uptake rates: mol·L⁻¹·s⁻¹ (converted from mmol·gDW⁻¹·h⁻¹)
- Cell density: cells·L⁻¹ or spots·L⁻¹
"""

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from rich.console import Console

console = Console()


class KineticsModel:
    """Kinetics model for cellular metabolism."""

    def __init__(self):
        """Initialize kinetics model."""
        self.kinetic_parameters: Dict[str, Dict[str, float]] = {}
        self.reaction_rates: Dict[str, callable] = {}

    def add_michaelis_menten(
        self,
        reaction_name: str,
        substrate: str,
        vmax: float,
        km: float,
        inhibitor: Optional[str] = None,
        ki: Optional[float] = None,
    ) -> None:
        """
        Add Michaelis-Menten kinetics for a reaction.

        Args:
            reaction_name: Name of the reaction
            substrate: Name of the substrate
            vmax: Maximum velocity
            km: Michaelis constant
            inhibitor: Name of inhibitor (optional)
            ki: Inhibition constant (optional)
        """
        console.print(f"Adding Michaelis-Menten kinetics for {reaction_name}")

        # Store parameters
        self.kinetic_parameters[reaction_name] = {
            "type": "michaelis_menten",
            "substrate": substrate,
            "vmax": vmax,
            "km": km,
        }

        if inhibitor and ki:
            self.kinetic_parameters[reaction_name]["inhibitor"] = inhibitor
            self.kinetic_parameters[reaction_name]["ki"] = ki

        # Create rate function
        def rate_function(concentrations):
            if substrate not in concentrations:
                return 0.0

            substrate_conc = concentrations[substrate]
            rate = (vmax * substrate_conc) / (km + substrate_conc)

            # Apply competitive inhibition if present
            if inhibitor and ki and inhibitor in concentrations:
                inhibitor_conc = concentrations[inhibitor]
                rate = rate / (1 + inhibitor_conc / ki)

            return rate

        self.reaction_rates[reaction_name] = rate_function

    def add_hill_kinetics(
        self,
        reaction_name: str,
        substrate: str,
        vmax: float,
        km: float,
        hill_coefficient: float = 1.0,
    ) -> None:
        """
        Add Hill kinetics for a reaction.

        Args:
            reaction_name: Name of the reaction
            substrate: Name of the substrate
            vmax: Maximum velocity
            km: Half-saturation constant
            hill_coefficient: Hill coefficient
        """
        console.print(f"Adding Hill kinetics for {reaction_name}")

        # Store parameters
        self.kinetic_parameters[reaction_name] = {
            "type": "hill",
            "substrate": substrate,
            "vmax": vmax,
            "km": km,
            "hill_coefficient": hill_coefficient,
        }

        # Create rate function
        def rate_function(concentrations):
            if substrate not in concentrations:
                return 0.0

            substrate_conc = concentrations[substrate]
            rate = (vmax * substrate_conc ** hill_coefficient) / (
                km ** hill_coefficient + substrate_conc ** hill_coefficient
            )

            return rate

        self.reaction_rates[reaction_name] = rate_function

    def add_first_order(
        self,
        reaction_name: str,
        substrate: str,
        rate_constant: float,
    ) -> None:
        """
        Add first-order kinetics for a reaction.

        Args:
            reaction_name: Name of the reaction
            substrate: Name of the substrate
            rate_constant: First-order rate constant
        """
        console.print(f"Adding first-order kinetics for {reaction_name}")

        # Store parameters
        self.kinetic_parameters[reaction_name] = {
            "type": "first_order",
            "substrate": substrate,
            "rate_constant": rate_constant,
        }

        # Create rate function
        def rate_function(concentrations):
            if substrate not in concentrations:
                return 0.0

            substrate_conc = concentrations[substrate]
            return rate_constant * substrate_conc

        self.reaction_rates[reaction_name] = rate_function

    def add_zero_order(
        self,
        reaction_name: str,
        rate_constant: float,
    ) -> None:
        """
        Add zero-order kinetics for a reaction.

        Args:
            reaction_name: Name of the reaction
            rate_constant: Zero-order rate constant
        """
        console.print(f"Adding zero-order kinetics for {reaction_name}")

        # Store parameters
        self.kinetic_parameters[reaction_name] = {
            "type": "zero_order",
            "rate_constant": rate_constant,
        }

        # Create rate function
        def rate_function(concentrations):
            return rate_constant

        self.reaction_rates[reaction_name] = rate_function

    def calculate_reaction_rate(
        self, 
        reaction_name: str, 
        concentrations: Dict[str, float]
    ) -> float:
        """
        Calculate reaction rate for a given reaction.

        Args:
            reaction_name: Name of the reaction
            concentrations: Dictionary of metabolite concentrations

        Returns:
            Reaction rate
        """
        if reaction_name not in self.reaction_rates:
            console.print(f"Warning: Reaction {reaction_name} not found")
            return 0.0

        return self.reaction_rates[reaction_name](concentrations)

    def calculate_all_rates(
        self, 
        concentrations: Dict[str, float]
    ) -> Dict[str, float]:
        """
        Calculate reaction rates for all reactions.

        Args:
            concentrations: Dictionary of metabolite concentrations

        Returns:
            Dictionary mapping reaction names to rates
        """
        rates = {}
        for reaction_name in self.reaction_rates:
            rates[reaction_name] = self.calculate_reaction_rate(reaction_name, concentrations)
        return rates

    def get_consumption_rate(
        self,
        metabolite: str,
        concentrations: Dict[str, float],
    ) -> float:
        """
        Calculate consumption rate of a metabolite.

        Args:
            metabolite: Name of the metabolite
            concentrations: Dictionary of metabolite concentrations

        Returns:
            Total consumption rate
        """
        total_consumption = 0.0

        for reaction_name, params in self.kinetic_parameters.items():
            if params.get("substrate") == metabolite:
                rate = self.calculate_reaction_rate(reaction_name, concentrations)
                total_consumption += rate

        return total_consumption

    def get_production_rate(
        self,
        metabolite: str,
        concentrations: Dict[str, float],
        stoichiometry: Dict[str, Dict[str, float]],
    ) -> float:
        """
        Calculate production rate of a metabolite.

        Args:
            metabolite: Name of the metabolite
            concentrations: Dictionary of metabolite concentrations
            stoichiometry: Stoichiometry matrix

        Returns:
            Total production rate
        """
        total_production = 0.0

        for reaction_name in self.reaction_rates:
            if reaction_name in stoichiometry and metabolite in stoichiometry[reaction_name]:
                rate = self.calculate_reaction_rate(reaction_name, concentrations)
                stoichiometric_coefficient = stoichiometry[reaction_name][metabolite]
                
                if stoichiometric_coefficient > 0:  # Product
                    total_production += rate * stoichiometric_coefficient

        return total_production

    def simulate_time_course(
        self,
        initial_concentrations: Dict[str, float],
        time_points: np.ndarray,
        stoichiometry: Dict[str, Dict[str, float]],
    ) -> Dict[str, np.ndarray]:
        """
        Simulate time course of metabolite concentrations.

        Args:
            initial_concentrations: Initial concentrations
            time_points: Array of time points
            stoichiometry: Stoichiometry matrix

        Returns:
            Dictionary mapping metabolite names to concentration time series
        """
        console.print("Simulating time course...")

        # Initialize concentrations
        concentrations = initial_concentrations.copy()
        results = {metabolite: [] for metabolite in concentrations.keys()}

        # Save initial concentrations
        for metabolite, conc in concentrations.items():
            results[metabolite].append(conc)

        # Time integration (Euler method)
        dt = time_points[1] - time_points[0] if len(time_points) > 1 else 0.1

        for t in time_points[1:]:
            # Calculate derivatives
            derivatives = {}
            
            for metabolite in concentrations:
                production = self.get_production_rate(metabolite, concentrations, stoichiometry)
                consumption = self.get_consumption_rate(metabolite, concentrations)
                derivatives[metabolite] = production - consumption

            # Update concentrations
            for metabolite in concentrations:
                concentrations[metabolite] += derivatives[metabolite] * dt
                concentrations[metabolite] = max(0, concentrations[metabolite])  # Non-negative

            # Save concentrations
            for metabolite, conc in concentrations.items():
                results[metabolite].append(conc)

        # Convert to numpy arrays
        for metabolite in results:
            results[metabolite] = np.array(results[metabolite])

        console.print("Time course simulation completed")
        return results

    def get_steady_state(
        self,
        initial_concentrations: Dict[str, float],
        stoichiometry: Dict[str, Dict[str, float]],
        tolerance: float = 1e-6,
        max_iterations: int = 1000,
    ) -> Dict[str, float]:
        """
        Find steady state concentrations.

        Args:
            initial_concentrations: Initial concentrations
            stoichiometry: Stoichiometry matrix
            tolerance: Convergence tolerance
            max_iterations: Maximum number of iterations

        Returns:
            Dictionary of steady state concentrations
        """
        console.print("Finding steady state...")

        concentrations = initial_concentrations.copy()

        for iteration in range(max_iterations):
            # Calculate derivatives
            derivatives = {}
            
            for metabolite in concentrations:
                production = self.get_production_rate(metabolite, concentrations, stoichiometry)
                consumption = self.get_consumption_rate(metabolite, concentrations)
                derivatives[metabolite] = production - consumption

            # Check convergence
            max_derivative = max(abs(d) for d in derivatives.values())
            if max_derivative < tolerance:
                console.print(f"Steady state found after {iteration + 1} iterations")
                break

            # Update concentrations (simple Euler step)
            dt = 0.01  # Small time step
            for metabolite in concentrations:
                concentrations[metabolite] += derivatives[metabolite] * dt
                concentrations[metabolite] = max(0, concentrations[metabolite])

        else:
            console.print("Warning: Steady state not reached within max iterations")

        return concentrations

    def export_kinetics_data(self, output_path: str) -> None:
        """
        Export kinetics parameters and functions.

        Args:
            output_path: Path to save the kinetics data
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        console.print(f"Exporting kinetics data to {output_path}")

        # Export parameters
        params_data = []
        for reaction_name, params in self.kinetic_parameters.items():
            params_data.append({
                "reaction_name": reaction_name,
                "type": params["type"],
                **{k: v for k, v in params.items() if k != "type"}
            })

        params_df = pd.DataFrame(params_data)
        params_df.to_csv(output_path / "kinetic_parameters.csv", index=False)

        console.print("Kinetics data exported successfully")


class SpatialKineticsModel(KineticsModel):
    """
    Spatial kinetics model for consumption/production rates per cell/spot.
    
    Handles conversion of FBA uptake rates to spatial consumption rates
    with proper unit conversions and Michaelis-Menten kinetics.
    """

    def __init__(
        self,
        cell_density: float = 1e9,  # cells·L⁻¹
        spot_density: float = 1e6,  # spots·L⁻¹
        cell_volume: float = 1e-12,  # L per cell
        spot_volume: float = 1e-9,   # L per spot
    ):
        """
        Initialize spatial kinetics model.
        
        Args:
            cell_density: Cell density in cells·L⁻¹
            spot_density: Spot density in spots·L⁻¹
            cell_volume: Volume per cell in L
            spot_volume: Volume per spot in L
        """
        super().__init__()
        self.cell_density = cell_density
        self.spot_density = spot_density
        self.cell_volume = cell_volume
        self.spot_volume = spot_volume
        
        # Conversion factors
        self.mmol_per_mol = 1e3  # mmol/mol
        self.h_per_s = 3600     # s/h
        self.gDW_per_cell = 1e-12  # gDW per cell (approximate)
        
        console.print(f"Initialized spatial kinetics with cell density: {cell_density:.2e} cells·L⁻¹")

    def convert_fba_to_spatial_rate(
        self,
        fba_rate: float,  # mmol·gDW⁻¹·h⁻¹
        spatial_type: str = "cell",  # "cell" or "spot"
        use_michaelis_menten: bool = False,
        concentration: float = 1.0,  # mol·L⁻¹
        vmax: Optional[float] = None,
        km: Optional[float] = None,
    ) -> float:
        """
        Convert FBA uptake rate to spatial consumption rate.
        
        Args:
            fba_rate: FBA uptake rate in mmol·gDW⁻¹·h⁻¹
            spatial_type: Type of spatial unit ("cell" or "spot")
            use_michaelis_menten: Whether to apply Michaelis-Menten kinetics
            concentration: Current substrate concentration in mol·L⁻¹
            vmax: Maximum velocity for Michaelis-Menten (mol·L⁻¹·s⁻¹)
            km: Michaelis constant for Michaelis-Menten (mol·L⁻¹)
            
        Returns:
            Spatial consumption rate in mol·L⁻¹·s⁻¹
        """
        # Convert FBA rate to mol·L⁻¹·s⁻¹
        # fba_rate [mmol·gDW⁻¹·h⁻¹] -> [mol·L⁻¹·s⁻¹]
        
        if spatial_type == "cell":
            density = self.cell_density
            volume = self.cell_volume
        elif spatial_type == "spot":
            density = self.spot_density
            volume = self.spot_volume
        else:
            raise ValueError(f"Unknown spatial type: {spatial_type}")
        
        # Conversion: mmol·gDW⁻¹·h⁻¹ -> mol·L⁻¹·s⁻¹
        # Step 1: mmol -> mol
        rate_mol = fba_rate / self.mmol_per_mol
        
        # Step 2: gDW⁻¹ -> L⁻¹ (using cell density and gDW per cell)
        rate_per_L = rate_mol * density * self.gDW_per_cell
        
        # Step 3: h⁻¹ -> s⁻¹
        rate_per_s = rate_per_L / self.h_per_s
        
        # Apply Michaelis-Menten kinetics if requested
        if use_michaelis_menten and vmax is not None and km is not None:
            rate_per_s = self._apply_michaelis_menten(rate_per_s, concentration, vmax, km)
        
        return rate_per_s

    def _apply_michaelis_menten(
        self,
        base_rate: float,  # mol·L⁻¹·s⁻¹
        concentration: float,  # mol·L⁻¹
        vmax: float,  # mol·L⁻¹·s⁻¹
        km: float,  # mol·L⁻¹
    ) -> float:
        """
        Apply Michaelis-Menten kinetics: v = Vmax * C / (Km + C)
        
        Args:
            base_rate: Base consumption rate
            concentration: Substrate concentration
            vmax: Maximum velocity
            km: Michaelis constant
            
        Returns:
            Michaelis-Menten rate in mol·L⁻¹·s⁻¹
        """
        if concentration <= 0:
            return 0.0
        
        # Michaelis-Menten equation
        mm_rate = vmax * concentration / (km + concentration)
        
        # Use the minimum of base rate and MM rate
        return min(base_rate, mm_rate)

    def calculate_spatial_consumption_rate(
        self,
        spot_id: str,
        substrate: str,
        fba_uptake_flux: float,  # mmol·gDW⁻¹·h⁻¹
        substrate_concentration: float,  # mol·L⁻¹
        spatial_type: str = "spot",
        use_michaelis_menten: bool = False,
        mm_params: Optional[Dict[str, float]] = None,
    ) -> float:
        """
        Calculate spatial consumption rate for a specific spot/cell.
        
        Args:
            spot_id: Identifier for the spatial location
            substrate: Name of the substrate
            fba_uptake_flux: FBA uptake flux in mmol·gDW⁻¹·h⁻¹
            substrate_concentration: Current substrate concentration in mol·L⁻¹
            spatial_type: Type of spatial unit ("cell" or "spot")
            use_michaelis_menten: Whether to apply Michaelis-Menten kinetics
            mm_params: Michaelis-Menten parameters (vmax, km)
            
        Returns:
            Consumption rate in mol·L⁻¹·s⁻¹
        """
        vmax = mm_params.get("vmax") if mm_params else None
        km = mm_params.get("km") if mm_params else None
        
        consumption_rate = self.convert_fba_to_spatial_rate(
            fba_rate=fba_uptake_flux,
            spatial_type=spatial_type,
            use_michaelis_menten=use_michaelis_menten,
            concentration=substrate_concentration,
            vmax=vmax,
            km=km,
        )
        
        console.print(f"Spot {spot_id}, {substrate}: consumption = {consumption_rate:.2e} mol·L⁻¹·s⁻¹")
        return consumption_rate

    def calculate_spatial_production_rate(
        self,
        spot_id: str,
        metabolite: str,
        fba_secretion_flux: float,  # mmol·gDW⁻¹·h⁻¹
        spatial_type: str = "spot",
    ) -> float:
        """
        Calculate spatial production rate for a specific spot/cell.
        
        Args:
            spot_id: Identifier for the spatial location
            metabolite: Name of the metabolite
            fba_secretion_flux: FBA secretion flux in mmol·gDW⁻¹·h⁻¹
            spatial_type: Type of spatial unit ("cell" or "spot")
            
        Returns:
            Production rate in mol·L⁻¹·s⁻¹
        """
        # Convert FBA secretion rate to spatial production rate
        production_rate = self.convert_fba_to_spatial_rate(
            fba_rate=abs(fba_secretion_flux),  # Ensure positive for production
            spatial_type=spatial_type,
            use_michaelis_menten=False,
        )
        
        console.print(f"Spot {spot_id}, {metabolite}: production = {production_rate:.2e} mol·L⁻¹·s⁻¹")
        return production_rate

    def get_spatial_rate_matrix(
        self,
        spot_coordinates: Dict[str, Tuple[float, float]],
        substrate_fluxes: Dict[str, Dict[str, float]],  # spot_id -> substrate -> flux
        substrate_concentrations: Dict[str, float],  # substrate -> concentration
        spatial_type: str = "spot",
        use_michaelis_menten: bool = False,
        mm_parameters: Optional[Dict[str, Dict[str, float]]] = None,  # substrate -> {vmax, km}
    ) -> Dict[str, Dict[str, float]]:
        """
        Calculate spatial rate matrix for all spots and substrates.
        
        Args:
            spot_coordinates: Dictionary mapping spot IDs to (x, y) coordinates
            substrate_fluxes: FBA fluxes for each spot and substrate
            substrate_concentrations: Current substrate concentrations
            spatial_type: Type of spatial unit ("cell" or "spot")
            use_michaelis_menten: Whether to apply Michaelis-Menten kinetics
            mm_parameters: Michaelis-Menten parameters for each substrate
            
        Returns:
            Dictionary mapping (spot_id, substrate) to spatial rate
        """
        console.print("Calculating spatial rate matrix...")
        
        spatial_rates = {}
        
        for spot_id, coordinates in spot_coordinates.items():
            spatial_rates[spot_id] = {}
            
            if spot_id in substrate_fluxes:
                for substrate, flux in substrate_fluxes[spot_id].items():
                    if substrate in substrate_concentrations:
                        concentration = substrate_concentrations[substrate]
                        mm_params = mm_parameters.get(substrate) if mm_parameters else None
                        
                        if flux < 0:  # Uptake (consumption)
                            rate = self.calculate_spatial_consumption_rate(
                                spot_id=spot_id,
                                substrate=substrate,
                                fba_uptake_flux=abs(flux),
                                substrate_concentration=concentration,
                                spatial_type=spatial_type,
                                use_michaelis_menten=use_michaelis_menten,
                                mm_params=mm_params,
                            )
                        else:  # Secretion (production)
                            rate = self.calculate_spatial_production_rate(
                                spot_id=spot_id,
                                metabolite=substrate,
                                fba_secretion_flux=flux,
                                spatial_type=spatial_type,
                            )
                        
                        spatial_rates[spot_id][substrate] = rate
        
        console.print(f"Spatial rate matrix calculated for {len(spatial_rates)} spots")
        return spatial_rates

    def export_spatial_rates(
        self,
        spatial_rates: Dict[str, Dict[str, float]],
        output_path: str,
    ) -> None:
        """
        Export spatial rate matrix to CSV.
        
        Args:
            spatial_rates: Spatial rate matrix
            output_path: Path to save the data
        """
        from pathlib import Path
        import pandas as pd
        
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Convert to DataFrame
        data = []
        for spot_id, substrates in spatial_rates.items():
            for substrate, rate in substrates.items():
                data.append({
                    "spot_id": spot_id,
                    "substrate": substrate,
                    "spatial_rate_mol_L_s": rate,
                })
        
        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)
        
        console.print(f"Spatial rates exported to {output_path}")


def create_default_michaelis_menten_params() -> Dict[str, Dict[str, float]]:
    """
    Create default Michaelis-Menten parameters for common metabolites.
    
    Returns:
        Dictionary mapping metabolite names to MM parameters
    """
    return {
        "glucose": {"vmax": 1e-3, "km": 5e-4},      # mol·L⁻¹·s⁻¹, mol·L⁻¹
        "oxygen": {"vmax": 2e-3, "km": 2e-5},       # mol·L⁻¹·s⁻¹, mol·L⁻¹
        "glutamine": {"vmax": 5e-4, "km": 1e-4},    # mol·L⁻¹·s⁻¹, mol·L⁻¹
        "lactate": {"vmax": 8e-4, "km": 3e-4},      # mol·L⁻¹·s⁻¹, mol·L⁻¹
        "pyruvate": {"vmax": 6e-4, "km": 2e-4},     # mol·L⁻¹·s⁻¹, mol·L⁻¹
    }


# Convenience function for module-level access
def simulate_kinetics(
    metabolite: str,
    initial_concentration: float,
    time_points: np.ndarray,
    production_rate: float = 0.0,
    consumption_rate: float = 0.0,
    **kwargs
) -> np.ndarray:
    """
    Convenience function to simulate metabolite kinetics.
    
    Args:
        metabolite: Name of the metabolite
        initial_concentration: Initial concentration
        time_points: Time points for simulation
        production_rate: Production rate
        consumption_rate: Consumption rate
        **kwargs: Additional arguments
        
    Returns:
        Concentration time series
    """
    model = KineticsModel(metabolite=metabolite, **kwargs)
    return model.simulate_concentration(
        initial_concentration=initial_concentration,
        time_points=time_points,
        production_rate=production_rate,
        consumption_rate=consumption_rate
    )
