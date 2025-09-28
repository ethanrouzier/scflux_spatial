"""
Human-GEM model loading and curation.

This module provides functionality to load and curate the Human-GEM
genome-scale metabolic model.
"""

import os
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import cobra
import pandas as pd
from rich.console import Console

console = Console()


class HumanGEM:
    """Human-GEM model loader and curator."""

    def __init__(
        self,
        model_url: Optional[str] = None,
        local_path: Optional[str] = None,
        biomass_reaction: str = "biomass_human",
        auto_download: bool = True,
    ):
        """
        Initialize Human-GEM loader.

        Args:
            model_url: URL to download Human-GEM model (if None, auto-detect latest)
            local_path: Local path to the model file
            biomass_reaction: Name of the biomass reaction
            auto_download: Whether to automatically download if model not found locally
        """
        self.model_url = model_url or self._get_latest_model_url()
        self.local_path = Path(local_path) if local_path else Path("data/models/Human-GEM.xml")
        self.biomass_reaction = biomass_reaction
        self.auto_download = auto_download
        self.model: Optional[cobra.Model] = None

    def _get_latest_model_url(self) -> str:
        """Get the URL for the latest Human-GEM model."""
        # Use the direct download link for the latest release
        return "https://github.com/SysBioChalmers/Human-GEM/raw/main/model/Human-GEM.xml"

    def download_model(self) -> Path:
        """
        Download Human-GEM model if not present locally.

        Returns:
            Path to the model file
        """
        if self.local_path.exists():
            console.print(f"Model already exists at {self.local_path}")
            return self.local_path

        if not self.auto_download:
            raise FileNotFoundError(
                f"Model not found at {self.local_path} and auto_download=False. "
                "Please provide a local path or enable auto_download."
            )

        console.print(f"Downloading Human-GEM model from {self.model_url}")
        self.local_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            urllib.request.urlretrieve(self.model_url, self.local_path)
            console.print(f"Model downloaded to {self.local_path}")
        except Exception as e:
            console.print(f"Error downloading model: {e}")
            raise

        return self.local_path

    def load_model(self) -> cobra.Model:
        """
        Load the Human-GEM model.

        Returns:
            COBRApy model object
        """
        if not self.local_path.exists():
            self.download_model()

        console.print(f"Loading Human-GEM model from {self.local_path}")

        try:
            self.model = cobra.io.read_sbml_model(str(self.local_path))
            console.print(f"Model loaded: {len(self.model.reactions)} reactions, "
                         f"{len(self.model.metabolites)} metabolites, "
                         f"{len(self.model.genes)} genes")
            
            # Harmonize compartments
            self._harmonize_compartments()
            
        except Exception as e:
            console.print(f"Error loading model: {e}")
            raise

        return self.model

    def _harmonize_compartments(self) -> None:
        """Harmonize metabolite compartments to standard notation."""
        console.print("Harmonizing compartments...")
        
        # Define compartment mapping
        compartment_mapping = {
            'c': 'c',      # Cytosol
            'm': 'c_mito', # Mitochondria
            'e': 'e',      # Extracellular
            'n': 'c',      # Nucleus -> Cytosol
            'r': 'c',      # Ribosome -> Cytosol
            'g': 'c',      # Golgi -> Cytosol
            'l': 'c',      # Lysosome -> Cytosol
            'x': 'c',      # Peroxisome -> Cytosol
        }
        
        # Update metabolite compartments
        for metabolite in self.model.metabolites:
            old_compartment = metabolite.compartment
            new_compartment = compartment_mapping.get(old_compartment, old_compartment)
            
            if old_compartment != new_compartment:
                metabolite.compartment = new_compartment
                metabolite.id = metabolite.id.replace(f'_{old_compartment}', f'_{new_compartment}')
        
        console.print("Compartment harmonization completed")

    def get_default_medium(self) -> Dict[str, float]:
        """
        Get default medium composition for human cells.
        
        Returns:
            Dictionary mapping exchange reaction IDs to bounds (negative for uptake)
        """
        return {
            # Essential nutrients
            'EX_glc__D_e': -10.0,      # Glucose uptake (mmol/gDW/h)
            'EX_o2_e': -1000.0,        # Oxygen uptake
            'EX_h2o_e': -1000.0,       # Water
            'EX_h_e': -1000.0,         # Protons
            'EX_na1_e': -1000.0,       # Sodium
            'EX_k_e': -1000.0,         # Potassium
            'EX_ca2_e': -1000.0,       # Calcium
            
            # Amino acids
            'EX_gln__L_e': -5.0,       # Glutamine
            'EX_ala__L_e': -1.0,       # Alanine
            'EX_arg__L_e': -1.0,       # Arginine
            'EX_asp__L_e': -1.0,       # Aspartate
            'EX_asn__L_e': -1.0,       # Asparagine
            'EX_cys__L_e': -1.0,       # Cysteine
            'EX_glu__L_e': -1.0,       # Glutamate
            'EX_gly_e': -1.0,          # Glycine
            'EX_his__L_e': -1.0,       # Histidine
            'EX_ile__L_e': -1.0,       # Isoleucine
            'EX_leu__L_e': -1.0,       # Leucine
            'EX_lys__L_e': -1.0,       # Lysine
            'EX_met__L_e': -1.0,       # Methionine
            'EX_phe__L_e': -1.0,       # Phenylalanine
            'EX_pro__L_e': -1.0,       # Proline
            'EX_ser__L_e': -1.0,       # Serine
            'EX_thr__L_e': -1.0,       # Threonine
            'EX_trp__L_e': -1.0,       # Tryptophan
            'EX_tyr__L_e': -1.0,       # Tyrosine
            'EX_val__L_e': -1.0,       # Valine
            
            # Vitamins and cofactors
            'EX_fol_e': -0.1,          # Folate
            'EX_ncam_e': -0.1,         # Nicotinamide
            'EX_pnto__R_e': -0.1,      # Pantothenate
            'EX_pydxn_e': -0.1,        # Pyridoxine
            'EX_ribflv_e': -0.1,       # Riboflavin
            'EX_thm_e': -0.1,          # Thiamine
            
            # Allow secretion
            'EX_lac__L_e': 1000.0,     # Lactate secretion
            'EX_co2_e': 1000.0,        # CO2 secretion
            'EX_hco3_e': 1000.0,       # Bicarbonate secretion
        }

    def curate_model(
        self,
        medium_composition: Optional[Dict[str, float]] = None,
        remove_blocked_reactions: bool = True,
        remove_unused_metabolites: bool = True,
        use_default_medium: bool = True,
    ) -> cobra.Model:
        """
        Curate the model for spatial flux analysis.

        Args:
            medium_composition: Dictionary of exchange reactions and bounds
            remove_blocked_reactions: Whether to remove blocked reactions
            remove_unused_metabolites: Whether to remove unused metabolites
            use_default_medium: Whether to use default medium if no custom medium provided

        Returns:
            Curated COBRApy model
        """
        if self.model is None:
            self.load_model()

        console.print("Curating Human-GEM model...")

        # Set medium composition
        if medium_composition:
            self._set_medium(medium_composition)
        elif use_default_medium:
            default_medium = self.get_default_medium()
            self._set_medium(default_medium)

        # Remove blocked reactions
        if remove_blocked_reactions:
            self._remove_blocked_reactions()

        # Remove unused metabolites
        if remove_unused_metabolites:
            self._remove_unused_metabolites()

        console.print("Model curation completed")
        return self.model

    def _set_medium(self, medium_composition: Dict[str, float]) -> None:
        """Set medium composition for the model."""
        console.print("Setting medium composition...")

        for reaction_id, bound in medium_composition.items():
            if reaction_id in self.model.reactions:
                reaction = self.model.reactions.get_by_id(reaction_id)
                if bound < 0:  # Uptake
                    reaction.lower_bound = bound
                else:  # Secretion
                    reaction.upper_bound = bound
            else:
                console.print(f"Warning: Reaction {reaction_id} not found in model")

    def _remove_blocked_reactions(self) -> None:
        """Remove blocked reactions from the model."""
        console.print("Removing blocked reactions...")

        # Find blocked reactions
        blocked_reactions = cobra.flux_analysis.find_blocked_reactions(self.model)
        console.print(f"Found {len(blocked_reactions)} blocked reactions")

        # Remove blocked reactions
        self.model.remove_reactions(blocked_reactions)
        console.print(f"Removed {len(blocked_reactions)} blocked reactions")

    def _remove_unused_metabolites(self) -> None:
        """Remove unused metabolites from the model."""
        console.print("Removing unused metabolites...")

        # Find unused metabolites
        unused_metabolites = []
        for metabolite in self.model.metabolites:
            if len(metabolite.reactions) == 0:
                unused_metabolites.append(metabolite)

        console.print(f"Found {len(unused_metabolites)} unused metabolites")

        # Remove unused metabolites
        self.model.remove_metabolites(unused_metabolites)
        console.print(f"Removed {len(unused_metabolites)} unused metabolites")

    def get_exchange_reactions(self) -> List[str]:
        """
        Get list of exchange reactions in the model.

        Returns:
            List of exchange reaction IDs
        """
        if self.model is None:
            self.load_model()

        exchange_reactions = []
        for reaction in self.model.reactions:
            # Check if it's an exchange reaction (single metabolite, coefficient 1 or -1)
            if len(reaction.metabolites) == 1:
                coeff = list(reaction.metabolites.values())[0]
                if abs(coeff) == 1:
                    exchange_reactions.append(reaction.id)

        return exchange_reactions

    def get_reaction_genes(self, reaction_id: str) -> List[str]:
        """
        Get genes associated with a reaction.

        Args:
            reaction_id: ID of the reaction

        Returns:
            List of gene IDs
        """
        if self.model is None:
            self.load_model()

        if reaction_id not in self.model.reactions:
            return []

        reaction = self.model.reactions.get_by_id(reaction_id)
        return [gene.id for gene in reaction.genes]

    def get_gene_reactions(self, gene_id: str) -> List[str]:
        """
        Get reactions associated with a gene.

        Args:
            gene_id: ID of the gene

        Returns:
            List of reaction IDs
        """
        if self.model is None:
            self.load_model()

        if gene_id not in self.model.genes:
            return []

        gene = self.model.genes.get_by_id(gene_id)
        return [reaction.id for reaction in gene.reactions]

    def get_central_carbon_reactions(self) -> List[str]:
        """
        Get reactions involved in central carbon metabolism.

        Returns:
            List of central carbon reaction IDs
        """
        if self.model is None:
            self.load_model()

        # Keywords for central carbon metabolism
        keywords = [
            "glycolysis", "glucose", "pyruvate", "lactate",
            "tca", "citrate", "oxaloacetate", "alpha_ketoglutarate",
            "pentose", "ribose", "gluconeogenesis"
        ]

        central_carbon_reactions = []
        for reaction in self.model.reactions:
            reaction_name = reaction.name.lower()
            reaction_id = reaction.id.lower()
            
            if any(keyword in reaction_name or keyword in reaction_id for keyword in keywords):
                central_carbon_reactions.append(reaction.id)

        return central_carbon_reactions

    def get_metabolite_info(self, metabolite_id: str) -> Dict[str, str]:
        """
        Get information about a metabolite.

        Args:
            metabolite_id: ID of the metabolite

        Returns:
            Dictionary with metabolite information
        """
        if self.model is None:
            self.load_model()

        if metabolite_id not in self.model.metabolites:
            return {}

        metabolite = self.model.metabolites.get_by_id(metabolite_id)
        return {
            "id": metabolite.id,
            "name": metabolite.name,
            "formula": metabolite.formula,
            "charge": metabolite.charge,
            "compartment": metabolite.compartment,
        }

    def export_reaction_network(self, output_path: str) -> None:
        """
        Export reaction network to CSV files.

        Args:
            output_path: Path to save the network files
        """
        if self.model is None:
            self.load_model()

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        console.print(f"Exporting reaction network to {output_path}")

        # Export reactions
        reactions_data = []
        for reaction in self.model.reactions:
            reactions_data.append({
                "id": reaction.id,
                "name": reaction.name,
                "equation": reaction.reaction,
                "lower_bound": reaction.lower_bound,
                "upper_bound": reaction.upper_bound,
                "objective_coefficient": reaction.objective_coefficient,
            })

        reactions_df = pd.DataFrame(reactions_data)
        reactions_df.to_csv(output_path / "reactions.csv", index=False)

        # Export metabolites
        metabolites_data = []
        for metabolite in self.model.metabolites:
            metabolites_data.append({
                "id": metabolite.id,
                "name": metabolite.name,
                "formula": metabolite.formula,
                "charge": metabolite.charge,
                "compartment": metabolite.compartment,
            })

        metabolites_df = pd.DataFrame(metabolites_data)
        metabolites_df.to_csv(output_path / "metabolites.csv", index=False)

        console.print("Reaction network exported successfully")
