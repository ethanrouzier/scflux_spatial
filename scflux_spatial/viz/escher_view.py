"""
Escher viewer for metabolic flux visualization.

This module provides functionality to create and display Escher maps
for visualizing central carbon metabolism fluxes.
"""

from typing import Dict, List, Optional, Tuple, Union

try:
    import escher
    ESCHER_AVAILABLE = True
except ImportError:
    ESCHER_AVAILABLE = False
    escher = None

import numpy as np
import pandas as pd
from rich.console import Console

console = Console()


class EscherViewer:
    """Escher map viewer for metabolic flux visualization."""

    def __init__(self, model=None):
        """
        Initialize Escher viewer.

        Args:
            model: COBRApy model object (optional)
        """
        self.model = model
        self.escher_maps = {}
        self.flux_data = {}

    def create_central_carbon_map(
        self,
        flux_data: Dict[str, float],
        map_name: str = "central_carbon",
        title: str = "Central Carbon Metabolism",
    ) -> Optional[escher.Builder]:
        """
        Create Escher map for central carbon metabolism.

        Args:
            flux_data: Dictionary mapping reaction IDs to flux values
            map_name: Name for the map
            title: Title for the map

        Returns:
            Escher Builder object or None if Escher not available
        """
        if not ESCHER_AVAILABLE:
            console.print("Warning: Escher not available, cannot create map")
            return None
            
        console.print(f"Creating Escher map: {map_name}")

        # Filter flux data for central carbon reactions
        central_carbon_reactions = self._get_central_carbon_reactions()
        filtered_fluxes = {
            rxn_id: flux for rxn_id, flux in flux_data.items()
            if rxn_id in central_carbon_reactions
        }

        # Create Escher map
        builder = escher.Builder(
            map_name=map_name,
            model=self.model,
            reaction_data=filtered_fluxes,
            title=title,
            height=600,
            width=800
        )

        # Store map and data
        self.escher_maps[map_name] = builder
        self.flux_data[map_name] = filtered_fluxes

        console.print(f"Escher map created with {len(filtered_fluxes)} reactions")
        return builder

    def create_glycolysis_map(
        self,
        flux_data: Dict[str, float],
        map_name: str = "glycolysis",
        title: str = "Glycolysis Pathway",
    ) -> Optional[escher.Builder]:
        """
        Create Escher map for glycolysis pathway.

        Args:
            flux_data: Dictionary mapping reaction IDs to flux values
            map_name: Name for the map
            title: Title for the map

        Returns:
            Escher Builder object or None if Escher not available
        """
        if not ESCHER_AVAILABLE:
            console.print("Warning: Escher not available, cannot create map")
            return None
            
        console.print(f"Creating glycolysis map: {map_name}")

        # Filter flux data for glycolysis reactions
        glycolysis_reactions = self._get_glycolysis_reactions()
        filtered_fluxes = {
            rxn_id: flux for rxn_id, flux in flux_data.items()
            if rxn_id in glycolysis_reactions
        }

        # Create Escher map
        builder = escher.Builder(
            map_name=map_name,
            model=self.model,
            reaction_data=filtered_fluxes,
            title=title,
            height=400,
            width=600
        )

        # Store map and data
        self.escher_maps[map_name] = builder
        self.flux_data[map_name] = filtered_fluxes

        console.print(f"Glycolysis map created with {len(filtered_fluxes)} reactions")
        return builder

    def create_tca_map(
        self,
        flux_data: Dict[str, float],
        map_name: str = "tca_cycle",
        title: str = "TCA Cycle",
    ) -> Optional[escher.Builder]:
        """
        Create Escher map for TCA cycle.

        Args:
            flux_data: Dictionary mapping reaction IDs to flux values
            map_name: Name for the map
            title: Title for the map

        Returns:
            Escher Builder object or None if Escher not available
        """
        if not ESCHER_AVAILABLE:
            console.print("Warning: Escher not available, cannot create map")
            return None
            
        console.print(f"Creating TCA cycle map: {map_name}")

        # Filter flux data for TCA cycle reactions
        tca_reactions = self._get_tca_reactions()
        filtered_fluxes = {
            rxn_id: flux for rxn_id, flux in flux_data.items()
            if rxn_id in tca_reactions
        }

        # Create Escher map
        builder = escher.Builder(
            map_name=map_name,
            model=self.model,
            reaction_data=filtered_fluxes,
            title=title,
            height=500,
            width=500
        )

        # Store map and data
        self.escher_maps[map_name] = builder
        self.flux_data[map_name] = filtered_fluxes

        console.print(f"TCA cycle map created with {len(filtered_fluxes)} reactions")
        return builder

    def update_flux_data(
        self,
        map_name: str,
        new_flux_data: Dict[str, float],
    ) -> None:
        """
        Update flux data for an existing map.

        Args:
            map_name: Name of the map to update
            new_flux_data: New flux data
        """
        if map_name not in self.escher_maps:
            console.print(f"Warning: Map {map_name} not found")
            return

        console.print(f"Updating flux data for {map_name}")

        # Update the map
        builder = self.escher_maps[map_name]
        builder.reaction_data = new_flux_data
        self.flux_data[map_name] = new_flux_data

        console.print(f"Updated {map_name} with {len(new_flux_data)} reactions")

    def save_map(
        self,
        map_name: str,
        output_path: str,
        format: str = "html",
    ) -> None:
        """
        Save Escher map to file.

        Args:
            map_name: Name of the map to save
            output_path: Path to save the map
            format: File format ('html', 'json')
        """
        if map_name not in self.escher_maps:
            console.print(f"Warning: Map {map_name} not found")
            return

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        builder = self.escher_maps[map_name]

        if format == "html":
            builder.save_html(output_path)
        elif format == "json":
            builder.save_json(output_path)
        else:
            raise ValueError(f"Unsupported format: {format}")

        console.print(f"Map {map_name} saved to {output_path}")

    def display_map(self, map_name: str) -> None:
        """
        Display Escher map in Jupyter notebook.

        Args:
            map_name: Name of the map to display
        """
        if map_name not in self.escher_maps:
            console.print(f"Warning: Map {map_name} not found")
            return

        builder = self.escher_maps[map_name]
        return builder.display()

    def get_map_statistics(self, map_name: str) -> Dict[str, Union[int, float]]:
        """
        Get statistics for a map.

        Args:
            map_name: Name of the map

        Returns:
            Dictionary with map statistics
        """
        if map_name not in self.escher_maps:
            return {}

        flux_data = self.flux_data[map_name]
        flux_values = list(flux_data.values())

        stats = {
            "total_reactions": len(flux_data),
            "active_reactions": len([f for f in flux_values if abs(f) > 1e-6]),
            "max_flux": max(flux_values) if flux_values else 0,
            "min_flux": min(flux_values) if flux_values else 0,
            "mean_flux": np.mean(flux_values) if flux_values else 0,
            "total_flux": sum(flux_values) if flux_values else 0,
        }

        return stats

    def _get_central_carbon_reactions(self) -> List[str]:
        """Get list of central carbon metabolism reactions."""
        # Common central carbon metabolism reactions
        central_carbon_reactions = [
            # Glycolysis
            "HEX1", "PGI", "PFK", "FBA", "TPI", "GAPD", "PGK", "PGM", "ENO", "PYK",
            # Gluconeogenesis
            "PC", "PCK", "FBP", "G6PC",
            # Pentose Phosphate Pathway
            "G6PDH2r", "PGL", "GND", "RPE", "RPI", "TKT1", "TKT2", "TALA",
            # TCA Cycle
            "CS", "ACONT", "ICDHyr", "AKGDH", "SUCOAS", "SUCDi", "FUM", "MDH",
            # Pyruvate metabolism
            "PDH", "LDH_L", "ME1", "ME2",
            # Exchange reactions
            "EX_glc__D_e", "EX_lac__L_e", "EX_pyr_e", "EX_co2_e",
        ]
        return central_carbon_reactions

    def _get_glycolysis_reactions(self) -> List[str]:
        """Get list of glycolysis reactions."""
        glycolysis_reactions = [
            "HEX1", "PGI", "PFK", "FBA", "TPI", "GAPD", "PGK", "PGM", "ENO", "PYK",
            "EX_glc__D_e", "EX_pyr_e", "EX_atp_e", "EX_adp_e", "EX_nad_e", "EX_nadh_e"
        ]
        return glycolysis_reactions

    def _get_tca_reactions(self) -> List[str]:
        """Get list of TCA cycle reactions."""
        tca_reactions = [
            "CS", "ACONT", "ICDHyr", "AKGDH", "SUCOAS", "SUCDi", "FUM", "MDH",
            "EX_pyr_e", "EX_co2_e", "EX_nad_e", "EX_nadh_e", "EX_fad_e", "EX_fadh2_e",
            "EX_atp_e", "EX_adp_e", "EX_gdp_e", "EX_gtp_e"
        ]
        return tca_reactions

    def create_custom_map(
        self,
        flux_data: Dict[str, float],
        reaction_list: List[str],
        map_name: str,
        title: str,
        width: int = 800,
        height: int = 600,
    ) -> Optional[escher.Builder]:
        """
        Create custom Escher map with specified reactions.

        Args:
            flux_data: Dictionary mapping reaction IDs to flux values
            reaction_list: List of reaction IDs to include
            map_name: Name for the map
            title: Title for the map
            width: Map width
            height: Map height

        Returns:
            Escher Builder object or None if Escher not available
        """
        if not ESCHER_AVAILABLE:
            console.print("Warning: Escher not available, cannot create map")
            return None
            
        console.print(f"Creating custom map: {map_name}")

        # Filter flux data for specified reactions
        filtered_fluxes = {
            rxn_id: flux for rxn_id, flux in flux_data.items()
            if rxn_id in reaction_list
        }

        # Create Escher map
        builder = escher.Builder(
            map_name=map_name,
            model=self.model,
            reaction_data=filtered_fluxes,
            title=title,
            height=height,
            width=width
        )

        # Store map and data
        self.escher_maps[map_name] = builder
        self.flux_data[map_name] = filtered_fluxes

        console.print(f"Custom map created with {len(filtered_fluxes)} reactions")
        return builder

    def export_flux_summary(self, output_path: str) -> None:
        """
        Export flux summary for all maps.

        Args:
            output_path: Path to save the summary
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        console.print(f"Exporting flux summary to {output_path}")

        # Create summary data
        summary_data = []
        for map_name, flux_data in self.flux_data.items():
            stats = self.get_map_statistics(map_name)
            summary_data.append({
                "map_name": map_name,
                **stats
            })

        # Export summary
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(output_path / "flux_summary.csv", index=False)

        # Export detailed flux data for each map
        for map_name, flux_data in self.flux_data.items():
            flux_df = pd.DataFrame([
                {"reaction_id": rxn_id, "flux": flux}
                for rxn_id, flux in flux_data.items()
            ])
            flux_df.to_csv(output_path / f"{map_name}_fluxes.csv", index=False)

        console.print("Flux summary exported successfully")

    def get_available_maps(self) -> List[str]:
        """Get list of available maps."""
        return list(self.escher_maps.keys())

    def remove_map(self, map_name: str) -> None:
        """
        Remove a map.

        Args:
            map_name: Name of the map to remove
        """
        if map_name in self.escher_maps:
            del self.escher_maps[map_name]
            del self.flux_data[map_name]
            console.print(f"Removed map: {map_name}")
        else:
            console.print(f"Warning: Map {map_name} not found")

    def create_cell_type_flux_map(
        self,
        cell_type_fluxes: Dict[str, Dict[str, float]],
        map_type: str = "central_carbon",
        title_prefix: str = "Cell Type",
        width: int = 800,
        height: int = 600,
    ) -> Dict[str, Optional[escher.Builder]]:
        """
        Create Escher maps showing average pFBA fluxes by cell type/cluster.
        
        Args:
            cell_type_fluxes: Dictionary mapping cell_type -> reaction_id -> flux
            map_type: Type of map to create ("central_carbon", "glycolysis", "tca")
            title_prefix: Prefix for map titles
            width: Map width
            height: Map height
            
        Returns:
            Dictionary mapping cell_type -> Escher Builder or None if Escher not available
        """
        if not ESCHER_AVAILABLE:
            console.print("Warning: Escher not available for cell type flux maps")
            return {}
        
        console.print(f"Creating {map_type} flux maps for {len(cell_type_fluxes)} cell types")
        
        # Get reaction list based on map type
        if map_type == "central_carbon":
            reaction_list = self._get_central_carbon_reactions()
        elif map_type == "glycolysis":
            reaction_list = self._get_glycolysis_reactions()
        elif map_type == "tca":
            reaction_list = self._get_tca_reactions()
        else:
            raise ValueError(f"Unknown map type: {map_type}")
        
        cell_type_maps = {}
        
        for cell_type, fluxes in cell_type_fluxes.items():
            # Filter fluxes for this map type
            filtered_fluxes = {
                rxn_id: flux for rxn_id, flux in fluxes.items()
                if rxn_id in reaction_list
            }
            
            # Create title
            title = f"{title_prefix}: {cell_type} - {map_type.replace('_', ' ').title()}"
            
            # Create map name
            map_name = f"{cell_type}_{map_type}"
            
            # Create Escher map
            builder = escher.Builder(
                map_name=map_name,
                model=self.model,
                reaction_data=filtered_fluxes,
                title=title,
                height=height,
                width=width
            )
            
            cell_type_maps[cell_type] = builder
            
            # Store in main maps
            self.escher_maps[map_name] = builder
            self.flux_data[map_name] = filtered_fluxes
            
            console.print(f"  Created map for {cell_type}: {len(filtered_fluxes)} reactions")
        
        return cell_type_maps

    def create_flux_comparison_map(
        self,
        cell_type_fluxes: Dict[str, Dict[str, float]],
        map_type: str = "central_carbon",
        comparison_method: str = "difference",
        reference_type: Optional[str] = None,
        title: str = "Flux Comparison",
        width: int = 800,
        height: int = 600,
    ) -> Optional[escher.Builder]:
        """
        Create Escher map comparing fluxes between cell types.
        
        Args:
            cell_type_fluxes: Dictionary mapping cell_type -> reaction_id -> flux
            map_type: Type of map to create
            comparison_method: Method for comparison ("difference", "ratio", "zscore")
            reference_type: Reference cell type for difference/ratio comparisons
            title: Map title
            width: Map width
            height: Map height
            
        Returns:
            Escher Builder with comparison data or None if Escher not available
        """
        if not ESCHER_AVAILABLE:
            console.print("Warning: Escher not available for flux comparison maps")
            return None
        
        console.print(f"Creating flux comparison map: {comparison_method}")
        
        # Get reaction list
        if map_type == "central_carbon":
            reaction_list = self._get_central_carbon_reactions()
        elif map_type == "glycolysis":
            reaction_list = self._get_glycolysis_reactions()
        elif map_type == "tca":
            reaction_list = self._get_tca_reactions()
        else:
            raise ValueError(f"Unknown map type: {map_type}")
        
        # Calculate comparison fluxes
        comparison_fluxes = {}
        
        if comparison_method == "difference" and reference_type is not None:
            # Calculate difference from reference
            ref_fluxes = cell_type_fluxes.get(reference_type, {})
            for cell_type, fluxes in cell_type_fluxes.items():
                if cell_type != reference_type:
                    for rxn_id in reaction_list:
                        if rxn_id in fluxes and rxn_id in ref_fluxes:
                            comparison_fluxes[f"{cell_type}_vs_{reference_type}_{rxn_id}"] = (
                                fluxes[rxn_id] - ref_fluxes[rxn_id]
                            )
        
        elif comparison_method == "ratio" and reference_type is not None:
            # Calculate ratio to reference
            ref_fluxes = cell_type_fluxes.get(reference_type, {})
            for cell_type, fluxes in cell_type_fluxes.items():
                if cell_type != reference_type:
                    for rxn_id in reaction_list:
                        if rxn_id in fluxes and rxn_id in ref_fluxes:
                            if abs(ref_fluxes[rxn_id]) > 1e-6:
                                ratio = fluxes[rxn_id] / ref_fluxes[rxn_id]
                                comparison_fluxes[f"{cell_type}_vs_{reference_type}_{rxn_id}"] = ratio
        
        elif comparison_method == "zscore":
            # Calculate z-scores across cell types
            for rxn_id in reaction_list:
                flux_values = []
                for fluxes in cell_type_fluxes.values():
                    if rxn_id in fluxes:
                        flux_values.append(fluxes[rxn_id])
                
                if len(flux_values) > 1:
                    mean_flux = np.mean(flux_values)
                    std_flux = np.std(flux_values)
                    if std_flux > 1e-6:
                        for cell_type, fluxes in cell_type_fluxes.items():
                            if rxn_id in fluxes:
                                zscore = (fluxes[rxn_id] - mean_flux) / std_flux
                                comparison_fluxes[f"{cell_type}_{rxn_id}"] = zscore
        
        # Create comparison map
        map_name = f"comparison_{comparison_method}_{map_type}"
        builder = escher.Builder(
            map_name=map_name,
            model=self.model,
            reaction_data=comparison_fluxes,
            title=f"{title} - {comparison_method.title()}",
            height=height,
            width=width
        )
        
        # Store map
        self.escher_maps[map_name] = builder
        self.flux_data[map_name] = comparison_fluxes
        
        console.print(f"Comparison map created with {len(comparison_fluxes)} flux comparisons")
        return builder

    def create_flux_statistics_summary(
        self,
        cell_type_fluxes: Dict[str, Dict[str, float]],
        output_path: str,
    ) -> None:
        """
        Create and export comprehensive flux statistics summary.
        
        Args:
            cell_type_fluxes: Dictionary mapping cell_type -> reaction_id -> flux
            output_path: Path to save the summary
        """
        from pathlib import Path
        
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        console.print(f"Creating flux statistics summary...")
        
        # Collect all unique reactions
        all_reactions = set()
        for fluxes in cell_type_fluxes.values():
            all_reactions.update(fluxes.keys())
        
        # Create comprehensive statistics
        stats_data = []
        
        for cell_type, fluxes in cell_type_fluxes.items():
            flux_values = list(fluxes.values())
            
            # Basic statistics
            stats = {
                "cell_type": cell_type,
                "total_reactions": len(fluxes),
                "active_reactions": len([f for f in flux_values if abs(f) > 1e-6]),
                "max_flux": max(flux_values) if flux_values else 0,
                "min_flux": min(flux_values) if flux_values else 0,
                "mean_flux": np.mean(flux_values) if flux_values else 0,
                "std_flux": np.std(flux_values) if flux_values else 0,
                "total_flux": sum(flux_values) if flux_values else 0,
            }
            
            # Pathway-specific statistics
            for pathway in ["central_carbon", "glycolysis", "tca"]:
                if pathway == "central_carbon":
                    pathway_reactions = self._get_central_carbon_reactions()
                elif pathway == "glycolysis":
                    pathway_reactions = self._get_glycolysis_reactions()
                elif pathway == "tca":
                    pathway_reactions = self._get_tca_reactions()
                
                pathway_fluxes = [fluxes.get(rxn, 0) for rxn in pathway_reactions if rxn in fluxes]
                
                if pathway_fluxes:
                    stats[f"{pathway}_mean_flux"] = np.mean(pathway_fluxes)
                    stats[f"{pathway}_total_flux"] = sum(pathway_fluxes)
                    stats[f"{pathway}_active_reactions"] = len([f for f in pathway_fluxes if abs(f) > 1e-6])
            
            stats_data.append(stats)
        
        # Export statistics
        stats_df = pd.DataFrame(stats_data)
        stats_df.to_csv(output_path / "cell_type_flux_statistics.csv", index=False)
        
        # Export detailed flux matrix
        flux_matrix_data = []
        for cell_type, fluxes in cell_type_fluxes.items():
            for reaction, flux in fluxes.items():
                flux_matrix_data.append({
                    "cell_type": cell_type,
                    "reaction_id": reaction,
                    "flux": flux
                })
        
        flux_matrix_df = pd.DataFrame(flux_matrix_data)
        flux_matrix_df.to_csv(output_path / "flux_matrix.csv", index=False)
        
        # Export pathway summaries
        for pathway in ["central_carbon", "glycolysis", "tca"]:
            if pathway == "central_carbon":
                pathway_reactions = self._get_central_carbon_reactions()
            elif pathway == "glycolysis":
                pathway_reactions = self._get_glycolysis_reactions()
            elif pathway == "tca":
                pathway_reactions = self._get_tca_reactions()
            
            pathway_data = []
            for cell_type, fluxes in cell_type_fluxes.items():
                for reaction in pathway_reactions:
                    if reaction in fluxes:
                        pathway_data.append({
                            "cell_type": cell_type,
                            "reaction_id": reaction,
                            "flux": fluxes[reaction]
                        })
            
            if pathway_data:
                pathway_df = pd.DataFrame(pathway_data)
                pathway_df.to_csv(output_path / f"{pathway}_fluxes.csv", index=False)
        
        console.print(f"Flux statistics summary exported to {output_path}")
