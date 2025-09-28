"""
Spatial mapping and visualization tools.

This module provides functionality to create heatmaps and overlays
for spatial flux analysis using plotly.
"""

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from rich.console import Console

console = Console()


class SpatialMapper:
    """Spatial mapping and visualization tools."""

    def __init__(
        self,
        colormap: str = "viridis",
        figure_size: Tuple[int, int] = (10, 8),
        dpi: int = 300,
    ):
        """
        Initialize spatial mapper.

        Args:
            colormap: Colormap for visualizations
            figure_size: Figure size (width, height)
            dpi: DPI for saved figures
        """
        self.colormap = colormap
        self.figure_size = figure_size
        self.dpi = dpi

    def create_metabolite_heatmap(
        self,
        concentrations: np.ndarray,
        metabolite_name: str,
        spatial_coords: Optional[np.ndarray] = None,
        title: Optional[str] = None,
    ) -> go.Figure:
        """
        Create heatmap for metabolite concentrations.

        Args:
            concentrations: Concentration array (2D or 1D)
            metabolite_name: Name of the metabolite
            spatial_coords: Spatial coordinates (optional)
            title: Plot title (optional)

        Returns:
            Plotly figure object
        """
        console.print(f"Creating heatmap for {metabolite_name}")

        # Reshape concentrations if needed
        if concentrations.ndim == 1:
            # Assume square grid
            grid_size = int(np.sqrt(len(concentrations)))
            concentrations = concentrations.reshape(grid_size, grid_size)

        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=concentrations,
            colorscale=self.colormap,
            colorbar=dict(title=f"{metabolite_name} (mM)"),
            hoverongaps=False
        ))

        # Update layout
        if title is None:
            title = f"{metabolite_name} Concentration Heatmap"

        fig.update_layout(
            title=title,
            xaxis_title="X Position",
            yaxis_title="Y Position",
            width=self.figure_size[0] * 100,
            height=self.figure_size[1] * 100,
        )

        return fig

    def create_flux_heatmap(
        self,
        fluxes: Dict[str, float],
        spatial_coords: np.ndarray,
        reaction_name: str,
        title: Optional[str] = None,
    ) -> go.Figure:
        """
        Create heatmap for reaction fluxes.

        Args:
            fluxes: Dictionary mapping positions to flux values
            spatial_coords: Spatial coordinates
            reaction_name: Name of the reaction
            title: Plot title (optional)

        Returns:
            Plotly figure object
        """
        console.print(f"Creating flux heatmap for {reaction_name}")

        # Extract flux values and coordinates
        x_coords = []
        y_coords = []
        flux_values = []

        for i, coord in enumerate(spatial_coords):
            if i in fluxes:
                x_coords.append(coord[0])
                y_coords.append(coord[1])
                flux_values.append(fluxes[i])

        # Create scatter plot with color mapping
        fig = go.Figure(data=go.Scatter(
            x=x_coords,
            y=y_coords,
            mode='markers',
            marker=dict(
                size=10,
                color=flux_values,
                colorscale=self.colormap,
                colorbar=dict(title=f"{reaction_name} Flux"),
                showscale=True
            ),
            text=[f"Flux: {flux:.4f}" for flux in flux_values],
            hovertemplate="X: %{x}<br>Y: %{y}<br>Flux: %{text}<extra></extra>"
        ))

        # Update layout
        if title is None:
            title = f"{reaction_name} Flux Distribution"

        fig.update_layout(
            title=title,
            xaxis_title="X Position",
            yaxis_title="Y Position",
            width=self.figure_size[0] * 100,
            height=self.figure_size[1] * 100,
        )

        return fig

    def create_multi_metabolite_plot(
        self,
        metabolite_data: Dict[str, np.ndarray],
        grid_size: Tuple[int, int],
        title: Optional[str] = None,
    ) -> go.Figure:
        """
        Create subplot with multiple metabolites.

        Args:
            metabolite_data: Dictionary mapping metabolite names to concentration arrays
            grid_size: Grid size for subplots
            title: Overall title (optional)

        Returns:
            Plotly figure with subplots
        """
        console.print(f"Creating multi-metabolite plot with {len(metabolite_data)} metabolites")

        n_metabolites = len(metabolite_data)
        n_cols = min(3, n_metabolites)  # Max 3 columns
        n_rows = (n_metabolites + n_cols - 1) // n_cols

        # Create subplots
        fig = make_subplots(
            rows=n_rows,
            cols=n_cols,
            subplot_titles=list(metabolite_data.keys()),
            specs=[[{"type": "heatmap"} for _ in range(n_cols)] for _ in range(n_rows)]
        )

        # Add heatmaps for each metabolite
        for i, (metabolite_name, concentrations) in enumerate(metabolite_data.items()):
            row = i // n_cols + 1
            col = i % n_cols + 1

            # Reshape if needed
            if concentrations.ndim == 1:
                grid_side = int(np.sqrt(len(concentrations)))
                concentrations = concentrations.reshape(grid_side, grid_side)

            fig.add_trace(
                go.Heatmap(
                    z=concentrations,
                    colorscale=self.colormap,
                    showscale=False,
                    name=metabolite_name
                ),
                row=row,
                col=col
            )

        # Update layout
        if title is None:
            title = "Multi-Metabolite Concentration Maps"

        fig.update_layout(
            title=title,
            width=self.figure_size[0] * 100 * n_cols,
            height=self.figure_size[1] * 100 * n_rows,
        )

        return fig

    def create_time_series_plot(
        self,
        time_series: Dict[str, List[float]],
        time_points: List[float],
        title: Optional[str] = None,
    ) -> go.Figure:
        """
        Create time series plot for metabolite concentrations.

        Args:
            time_series: Dictionary mapping metabolite names to time series
            time_points: Time points
            title: Plot title (optional)

        Returns:
            Plotly figure object
        """
        console.print(f"Creating time series plot for {len(time_series)} metabolites")

        fig = go.Figure()

        # Add traces for each metabolite
        for metabolite_name, concentrations in time_series.items():
            fig.add_trace(go.Scatter(
                x=time_points,
                y=concentrations,
                mode='lines',
                name=metabolite_name,
                line=dict(width=2)
            ))

        # Update layout
        if title is None:
            title = "Metabolite Concentration Time Series"

        fig.update_layout(
            title=title,
            xaxis_title="Time",
            yaxis_title="Concentration (mM)",
            width=self.figure_size[0] * 100,
            height=self.figure_size[1] * 100,
            hovermode='x unified'
        )

        return fig

    def create_spatial_flux_overlay(
        self,
        spatial_coords: np.ndarray,
        flux_data: Dict[str, np.ndarray],
        tissue_image: Optional[np.ndarray] = None,
        title: Optional[str] = None,
    ) -> go.Figure:
        """
        Create spatial overlay of flux data.

        Args:
            spatial_coords: Spatial coordinates
            flux_data: Dictionary mapping reaction names to flux arrays
            tissue_image: Tissue image for background (optional)
            title: Plot title (optional)

        Returns:
            Plotly figure object
        """
        console.print(f"Creating spatial flux overlay with {len(flux_data)} reactions")

        fig = go.Figure()

        # Add tissue image if provided
        if tissue_image is not None:
            fig.add_trace(go.Image(z=tissue_image, opacity=0.7))

        # Add flux data as scatter plots
        for reaction_name, fluxes in flux_data.items():
            fig.add_trace(go.Scatter(
                x=spatial_coords[:, 0],
                y=spatial_coords[:, 1],
                mode='markers',
                marker=dict(
                    size=8,
                    color=fluxes,
                    colorscale=self.colormap,
                    opacity=0.8,
                    showscale=True,
                    colorbar=dict(title=f"{reaction_name} Flux")
                ),
                name=reaction_name,
                visible=False  # Initially hidden
            ))

        # Make first trace visible
        if fig.data:
            fig.data[0].visible = True

        # Add dropdown menu
        buttons = []
        for i, reaction_name in enumerate(flux_data.keys()):
            visibility = [False] * len(flux_data)
            visibility[i] = True
            buttons.append(
                dict(
                    label=reaction_name,
                    method="update",
                    args=[{"visible": visibility}]
                )
            )

        fig.update_layout(
            updatemenus=[{
                "buttons": buttons,
                "direction": "down",
                "showactive": True,
                "x": 0.1,
                "xanchor": "left",
                "y": 1.15,
                "yanchor": "top"
            }]
        )

        # Update layout
        if title is None:
            title = "Spatial Flux Overlay"

        fig.update_layout(
            title=title,
            xaxis_title="X Position",
            yaxis_title="Y Position",
            width=self.figure_size[0] * 100,
            height=self.figure_size[1] * 100,
        )

        return fig

    def create_flux_network_plot(
        self,
        flux_data: Dict[str, float],
        reaction_network: Dict[str, List[str]],
        title: Optional[str] = None,
    ) -> go.Figure:
        """
        Create network plot showing flux connections.

        Args:
            flux_data: Dictionary mapping reaction names to flux values
            reaction_network: Network structure (reaction -> connected reactions)
            title: Plot title (optional)

        Returns:
            Plotly figure object
        """
        console.print("Creating flux network plot")

        # Create network layout (simplified force-directed layout)
        nodes = list(flux_data.keys())
        n_nodes = len(nodes)
        
        # Simple circular layout
        angles = np.linspace(0, 2*np.pi, n_nodes, endpoint=False)
        x_coords = np.cos(angles)
        y_coords = np.sin(angles)

        # Create node traces
        node_trace = go.Scatter(
            x=x_coords,
            y=y_coords,
            mode='markers+text',
            marker=dict(
                size=20,
                color=list(flux_data.values()),
                colorscale=self.colormap,
                colorbar=dict(title="Flux Value"),
                showscale=True
            ),
            text=nodes,
            textposition="middle center",
            name="Reactions"
        )

        # Create edge traces
        edge_x = []
        edge_y = []
        
        for reaction, connected_reactions in reaction_network.items():
            if reaction in nodes:
                source_idx = nodes.index(reaction)
                for connected_reaction in connected_reactions:
                    if connected_reaction in nodes:
                        target_idx = nodes.index(connected_reaction)
                        
                        edge_x.extend([x_coords[source_idx], x_coords[target_idx], None])
                        edge_y.extend([y_coords[source_idx], y_coords[target_idx], None])

        edge_trace = go.Scatter(
            x=edge_x,
            y=edge_y,
            mode='lines',
            line=dict(width=1, color='gray'),
            hoverinfo='none',
            showlegend=False
        )

        # Create figure
        fig = go.Figure(data=[edge_trace, node_trace])

        # Update layout
        if title is None:
            title = "Flux Network"

        fig.update_layout(
            title=title,
            width=self.figure_size[0] * 100,
            height=self.figure_size[1] * 100,
            showlegend=False,
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        )

        return fig

    def save_figure(self, fig: go.Figure, output_path: str, format: str = "html") -> None:
        """
        Save figure to file.

        Args:
            fig: Plotly figure object
            output_path: Path to save the figure
            format: File format ('html', 'png', 'pdf', 'svg')
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if format == "html":
            fig.write_html(output_path)
        elif format == "png":
            fig.write_image(output_path, width=self.figure_size[0] * self.dpi, 
                           height=self.figure_size[1] * self.dpi)
        elif format == "pdf":
            fig.write_image(output_path)
        elif format == "svg":
            fig.write_image(output_path)
        else:
            raise ValueError(f"Unsupported format: {format}")

        console.print(f"Figure saved to {output_path}")

    def export_spatial_data(
        self,
        output_path: str,
        spatial_data: Dict[str, np.ndarray],
        spatial_coords: Optional[np.ndarray] = None,
    ) -> None:
        """
        Export spatial data for external visualization.

        Args:
            output_path: Path to save the data
            spatial_data: Dictionary mapping names to spatial arrays
            spatial_coords: Spatial coordinates (optional)
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        console.print(f"Exporting spatial data to {output_path}")

        # Export each spatial dataset
        for name, data in spatial_data.items():
            if data.ndim == 1:
                # 1D data - save as CSV
                df = pd.DataFrame({"value": data})
                if spatial_coords is not None:
                    df["x"] = spatial_coords[:, 0]
                    df["y"] = spatial_coords[:, 1]
                df.to_csv(output_path / f"{name}_spatial.csv", index=False)
            else:
                # 2D data - save as numpy array
                np.save(output_path / f"{name}_spatial.npy", data)

        console.print("Spatial data exported successfully")

    def heatmap_concentration(
        self,
        adata,
        C_field: np.ndarray,
        field_name: str = "concentration",
        overlay_spots: bool = True,
        spot_size: int = 8,
        spot_opacity: float = 0.7,
        title: Optional[str] = None,
    ) -> go.Figure:
        """
        Create concentration heatmap with spot overlays.
        
        Args:
            adata: AnnData object with spatial coordinates
            C_field: 2D concentration field array
            field_name: Name of the field for display
            overlay_spots: Whether to overlay spot positions
            spot_size: Size of spot markers
            spot_opacity: Opacity of spot markers
            title: Plot title
            
        Returns:
            Plotly figure with heatmap and spot overlays
        """
        console.print(f"Creating concentration heatmap for {field_name}")
        
        # Create base heatmap
        fig = go.Figure(data=go.Heatmap(
            z=C_field,
            colorscale=self.colormap,
            colorbar=dict(title=f"{field_name} (mol·L⁻¹)"),
            hoverongaps=False,
            name="Concentration Field"
        ))
        
        # Add spot overlays if requested
        if overlay_spots and hasattr(adata, 'obs') and 'x' in adata.obs and 'y' in adata.obs:
            # Get spot coordinates
            x_spots = adata.obs['x'].values
            y_spots = adata.obs['y'].values
            
            # Normalize coordinates to field dimensions
            x_norm = (x_spots - x_spots.min()) / (x_spots.max() - x_spots.min()) * (C_field.shape[1] - 1)
            y_norm = (y_spots - y_spots.min()) / (y_spots.max() - y_spots.min()) * (C_field.shape[0] - 1)
            
            # Add spot scatter overlay
            fig.add_trace(go.Scatter(
                x=x_norm,
                y=y_norm,
                mode='markers',
                marker=dict(
                    size=spot_size,
                    color='white',
                    line=dict(width=1, color='black'),
                    opacity=spot_opacity
                ),
                name="Visium Spots",
                hovertemplate="Spot<br>X: %{x}<br>Y: %{y}<extra></extra>"
            ))
        
        # Update layout
        if title is None:
            title = f"{field_name.title()} Concentration Field"
        
        fig.update_layout(
            title=title,
            xaxis_title="X Position",
            yaxis_title="Y Position",
            width=self.figure_size[0] * 100,
            height=self.figure_size[1] * 100,
        )
        
        return fig

    def field_quiver(
        self,
        C_field: np.ndarray,
        field_name: str = "concentration",
        skip: int = 5,
        scale: float = 1.0,
        title: Optional[str] = None,
    ) -> go.Figure:
        """
        Create quiver plot showing gradient field.
        
        Args:
            C_field: 2D concentration field array
            field_name: Name of the field for display
            skip: Number of points to skip for arrow density
            scale: Scale factor for arrow length
            title: Plot title
            
        Returns:
            Plotly figure with quiver plot
        """
        console.print(f"Creating quiver plot for {field_name}")
        
        # Calculate gradients
        grad_y, grad_x = np.gradient(C_field)
        
        # Create meshgrid for positions
        y_coords, x_coords = np.mgrid[0:C_field.shape[0], 0:C_field.shape[1]]
        
        # Sample points for arrows
        x_quiver = x_coords[::skip, ::skip]
        y_quiver = y_coords[::skip, ::skip]
        u_quiver = grad_x[::skip, ::skip] * scale
        v_quiver = grad_y[::skip, ::skip] * scale
        
        # Create figure
        fig = go.Figure()
        
        # Add background heatmap
        fig.add_trace(go.Heatmap(
            z=C_field,
            colorscale='gray',
            opacity=0.7,
            showscale=False,
            name="Concentration"
        ))
        
        # Add quiver arrows
        for i in range(x_quiver.shape[0]):
            for j in range(x_quiver.shape[1]):
                if i == 0 and j == 0:  # Only add legend for first arrow
                    fig.add_trace(go.Scatter(
                        x=[x_quiver[i, j], x_quiver[i, j] + u_quiver[i, j]],
                        y=[y_quiver[i, j], y_quiver[i, j] + v_quiver[i, j]],
                        mode='lines+markers',
                        line=dict(color='red', width=2),
                        marker=dict(size=4, color='red'),
                        name="Gradient Vectors"
                    ))
                else:
                    fig.add_trace(go.Scatter(
                        x=[x_quiver[i, j], x_quiver[i, j] + u_quiver[i, j]],
                        y=[y_quiver[i, j], y_quiver[i, j] + v_quiver[i, j]],
                        mode='lines+markers',
                        line=dict(color='red', width=2),
                        marker=dict(size=4, color='red'),
                        showlegend=False
                    ))
        
        # Update layout
        if title is None:
            title = f"{field_name.title()} Gradient Field"
        
        fig.update_layout(
            title=title,
            xaxis_title="X Position",
            yaxis_title="Y Position",
            width=self.figure_size[0] * 100,
            height=self.figure_size[1] * 100,
        )
        
        return fig

    def spot_metric_map(
        self,
        adata,
        metric_name: str,
        metric_values: Optional[np.ndarray] = None,
        metric_key: Optional[str] = None,
        colormap: str = "viridis",
        title: Optional[str] = None,
        show_clusters: bool = False,
    ) -> go.Figure:
        """
        Create spatial map of spot-level metrics (e.g., ATP flux, lactate export).
        
        Args:
            adata: AnnData object with spatial coordinates
            metric_name: Name of the metric for display
            metric_values: Array of metric values (if not in adata.obs)
            metric_key: Key in adata.obs to use for metric values
            colormap: Colormap for the metric
            title: Plot title
            show_clusters: Whether to show cluster boundaries
            
        Returns:
            Plotly figure with spot metric map
        """
        console.print(f"Creating spot metric map for {metric_name}")
        
        # Get metric values
        if metric_values is not None:
            values = metric_values
        elif metric_key is not None and metric_key in adata.obs:
            values = adata.obs[metric_key].values
        else:
            raise ValueError("Must provide either metric_values or valid metric_key")
        
        # Get spatial coordinates
        if 'x' in adata.obs and 'y' in adata.obs:
            x_coords = adata.obs['x'].values
            y_coords = adata.obs['y'].values
        else:
            raise ValueError("Spatial coordinates (x, y) not found in adata.obs")
        
        # Create scatter plot
        fig = go.Figure(data=go.Scatter(
            x=x_coords,
            y=y_coords,
            mode='markers',
            marker=dict(
                size=12,
                color=values,
                colorscale=colormap,
                colorbar=dict(title=f"{metric_name}"),
                showscale=True,
                line=dict(width=1, color='black')
            ),
            text=[f"{metric_name}: {val:.3f}" for val in values],
            hovertemplate="X: %{x}<br>Y: %{y}<br>%{text}<extra></extra>",
            name=metric_name
        ))
        
        # Add cluster boundaries if requested
        if show_clusters and 'leiden' in adata.obs:
            clusters = adata.obs['leiden'].unique()
            colors = px.colors.qualitative.Set1[:len(clusters)]
            
            for i, cluster in enumerate(clusters):
                cluster_mask = adata.obs['leiden'] == cluster
                cluster_x = x_coords[cluster_mask]
                cluster_y = y_coords[cluster_mask]
                
                fig.add_trace(go.Scatter(
                    x=cluster_x,
                    y=cluster_y,
                    mode='markers',
                    marker=dict(
                        size=8,
                        color=colors[i],
                        opacity=0.3,
                        symbol='circle-open'
                    ),
                    name=f"Cluster {cluster}",
                    showlegend=True
                ))
        
        # Update layout
        if title is None:
            title = f"{metric_name.title()} Spatial Distribution"
        
        fig.update_layout(
            title=title,
            xaxis_title="X Position",
            yaxis_title="Y Position",
            width=self.figure_size[0] * 100,
            height=self.figure_size[1] * 100,
        )
        
        return fig

    def create_multi_field_comparison(
        self,
        adata,
        field_dict: Dict[str, np.ndarray],
        overlay_spots: bool = True,
        title: Optional[str] = None,
    ) -> go.Figure:
        """
        Create comparison plot of multiple concentration fields.
        
        Args:
            adata: AnnData object with spatial coordinates
            field_dict: Dictionary mapping field names to 2D arrays
            overlay_spots: Whether to overlay spot positions
            title: Plot title
            
        Returns:
            Plotly figure with subplots
        """
        console.print(f"Creating multi-field comparison with {len(field_dict)} fields")
        
        n_fields = len(field_dict)
        n_cols = min(3, n_fields)
        n_rows = (n_fields + n_cols - 1) // n_cols
        
        # Create subplots
        fig = make_subplots(
            rows=n_rows,
            cols=n_cols,
            subplot_titles=list(field_dict.keys()),
            specs=[[{"type": "heatmap"} for _ in range(n_cols)] for _ in range(n_rows)]
        )
        
        # Add heatmaps for each field
        for i, (field_name, field_data) in enumerate(field_dict.items()):
            row = i // n_cols + 1
            col = i % n_cols + 1
            
            fig.add_trace(
                go.Heatmap(
                    z=field_data,
                    colorscale=self.colormap,
                    showscale=(i == 0),  # Only show colorbar for first subplot
                    colorbar=dict(title="Concentration (mol·L⁻¹)") if i == 0 else None,
                    name=field_name
                ),
                row=row,
                col=col
            )
            
            # Add spot overlays if requested
            if overlay_spots and hasattr(adata, 'obs') and 'x' in adata.obs and 'y' in adata.obs:
                x_spots = adata.obs['x'].values
                y_spots = adata.obs['y'].values
                
                # Normalize coordinates
                x_norm = (x_spots - x_spots.min()) / (x_spots.max() - x_spots.min()) * (field_data.shape[1] - 1)
                y_norm = (y_spots - y_spots.min()) / (y_spots.max() - y_spots.min()) * (field_data.shape[0] - 1)
                
                fig.add_trace(
                    go.Scatter(
                        x=x_norm,
                        y=y_norm,
                        mode='markers',
                        marker=dict(size=6, color='white', line=dict(width=1, color='black')),
                        showlegend=False,
                        name=f"Spots_{field_name}"
                    ),
                    row=row,
                    col=col
                )
        
        # Update layout
        if title is None:
            title = "Multi-Field Concentration Comparison"
        
        fig.update_layout(
            title=title,
            width=self.figure_size[0] * 100 * n_cols,
            height=self.figure_size[1] * 100 * n_rows,
        )
        
        return fig
