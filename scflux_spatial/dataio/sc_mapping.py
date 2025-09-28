"""
Single-cell mapping interfaces for Tangram and cell2location.

This module provides optional interfaces to external tools for mapping
single-cell data to spatial coordinates.
"""

from typing import Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd
from rich.console import Console

console = Console()


class SCMapper:
    """Interface for single-cell to spatial mapping tools."""

    def __init__(self, method: str = "tangram"):
        """
        Initialize single-cell mapper.

        Args:
            method: Mapping method ('tangram' or 'cell2location')
        """
        self.method = method.lower()
        if self.method not in ["tangram", "cell2location"]:
            raise ValueError("Method must be 'tangram' or 'cell2location'")

    def map_cells_to_spatial(
        self,
        sc_data: np.ndarray,
        spatial_data: np.ndarray,
        sc_genes: np.ndarray,
        spatial_genes: np.ndarray,
        spatial_coords: np.ndarray,
        **kwargs
    ) -> Dict[str, np.ndarray]:
        """
        Map single-cell data to spatial coordinates.

        Args:
            sc_data: Single-cell expression matrix (n_cells, n_genes)
            spatial_data: Spatial expression matrix (n_spots, n_genes)
            sc_genes: Single-cell gene names
            spatial_genes: Spatial gene names
            spatial_coords: Spatial coordinates (n_spots, 2)
            **kwargs: Additional arguments for the mapping method

        Returns:
            Dictionary with mapping results
        """
        console.print(f"Mapping cells to spatial coordinates using {self.method}")

        if self.method == "tangram":
            return self._tangram_mapping(
                sc_data, spatial_data, sc_genes, spatial_genes, spatial_coords, **kwargs
            )
        elif self.method == "cell2location":
            return self._cell2location_mapping(
                sc_data, spatial_data, sc_genes, spatial_genes, spatial_coords, **kwargs
            )

    def _tangram_mapping(
        self,
        sc_data: np.ndarray,
        spatial_data: np.ndarray,
        sc_genes: np.ndarray,
        spatial_genes: np.ndarray,
        spatial_coords: np.ndarray,
        **kwargs
    ) -> Dict[str, np.ndarray]:
        """
        Tangram-based mapping (placeholder implementation).

        Note: This is a placeholder. In a real implementation, you would
        integrate with the actual Tangram package.
        """
        console.print("Using Tangram mapping (placeholder implementation)")

        # Find common genes
        common_genes = np.intersect1d(sc_genes, spatial_genes)
        sc_common_idx = np.where(np.isin(sc_genes, common_genes))[0]
        spatial_common_idx = np.where(np.isin(spatial_genes, common_genes))[0]

        # Subset data to common genes
        sc_subset = sc_data[:, sc_common_idx]
        spatial_subset = spatial_data[:, spatial_common_idx]

        # Simple correlation-based mapping (placeholder)
        n_sc_cells = sc_subset.shape[0]
        n_spatial_spots = spatial_subset.shape[0]

        # Compute correlations between each single cell and spatial spots
        mapping_probs = np.zeros((n_sc_cells, n_spatial_spots))
        for i in range(n_sc_cells):
            for j in range(n_spatial_spots):
                correlation = np.corrcoef(sc_subset[i], spatial_subset[j])[0, 1]
                mapping_probs[i, j] = max(0, correlation)  # Ensure non-negative

        # Normalize probabilities
        mapping_probs = mapping_probs / mapping_probs.sum(axis=1, keepdims=True)

        return {
            "mapping_probabilities": mapping_probs,
            "spatial_coordinates": spatial_coords,
            "common_genes": common_genes,
        }

    def _cell2location_mapping(
        self,
        sc_data: np.ndarray,
        spatial_data: np.ndarray,
        sc_genes: np.ndarray,
        spatial_genes: np.ndarray,
        spatial_coords: np.ndarray,
        **kwargs
    ) -> Dict[str, np.ndarray]:
        """
        Cell2location-based mapping (placeholder implementation).

        Note: This is a placeholder. In a real implementation, you would
        integrate with the actual cell2location package.
        """
        console.print("Using cell2location mapping (placeholder implementation)")

        # Find common genes
        common_genes = np.intersect1d(sc_genes, spatial_genes)
        sc_common_idx = np.where(np.isin(sc_genes, common_genes))[0]
        spatial_common_idx = np.where(np.isin(spatial_genes, common_genes))[0]

        # Subset data to common genes
        sc_subset = sc_data[:, sc_common_idx]
        spatial_subset = spatial_data[:, spatial_common_idx]

        # Simple deconvolution approach (placeholder)
        n_sc_cells = sc_subset.shape[0]
        n_spatial_spots = spatial_subset.shape[0]

        # Compute cell type abundances for each spatial spot
        cell_abundances = np.zeros((n_spatial_spots, n_sc_cells))
        for j in range(n_spatial_spots):
            spot_expr = spatial_subset[j]
            for i in range(n_sc_cells):
                cell_expr = sc_subset[i]
                # Simple dot product similarity
                similarity = np.dot(spot_expr, cell_expr) / (
                    np.linalg.norm(spot_expr) * np.linalg.norm(cell_expr) + 1e-8
                )
                cell_abundances[j, i] = max(0, similarity)

        # Normalize abundances
        cell_abundances = cell_abundances / cell_abundances.sum(axis=1, keepdims=True)

        return {
            "cell_abundances": cell_abundances,
            "spatial_coordinates": spatial_coords,
            "common_genes": common_genes,
        }

    def get_cell_type_profiles(
        self, mapping_results: Dict[str, np.ndarray]
    ) -> Dict[str, np.ndarray]:
        """
        Extract cell type profiles from mapping results.

        Args:
            mapping_results: Results from map_cells_to_spatial

        Returns:
            Dictionary with cell type profiles
        """
        console.print("Extracting cell type profiles from mapping results")

        if self.method == "tangram":
            # Extract dominant cell types for each spatial location
            mapping_probs = mapping_results["mapping_probabilities"]
            dominant_cells = np.argmax(mapping_probs, axis=0)
            return {"dominant_cell_types": dominant_cells}

        elif self.method == "cell2location":
            # Extract cell type abundances
            cell_abundances = mapping_results["cell_abundances"]
            return {"cell_type_abundances": cell_abundances}

    def validate_mapping(
        self, mapping_results: Dict[str, np.ndarray]
    ) -> Dict[str, float]:
        """
        Validate mapping results.

        Args:
            mapping_results: Results from map_cells_to_spatial

        Returns:
            Dictionary with validation metrics
        """
        console.print("Validating mapping results")

        validation_metrics = {}

        if self.method == "tangram":
            mapping_probs = mapping_results["mapping_probabilities"]
            # Check probability normalization
            prob_sums = mapping_probs.sum(axis=1)
            validation_metrics["probability_normalization"] = np.mean(
                np.abs(prob_sums - 1.0)
            )
            validation_metrics["max_probability"] = np.max(mapping_probs)
            validation_metrics["min_probability"] = np.min(mapping_probs)

        elif self.method == "cell2location":
            cell_abundances = mapping_results["cell_abundances"]
            # Check abundance normalization
            abundance_sums = cell_abundances.sum(axis=1)
            validation_metrics["abundance_normalization"] = np.mean(
                np.abs(abundance_sums - 1.0)
            )
            validation_metrics["max_abundance"] = np.max(cell_abundances)
            validation_metrics["min_abundance"] = np.min(cell_abundances)

        return validation_metrics
