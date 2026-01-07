"""
Advanced plotting functions for ADCIRC mesh and output data.
These functions handle common challenges like dry areas (NaN values) and
provide fallback methods when traditional plotting approaches fail.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.cm as cm
from matplotlib.colors import BoundaryNorm
import matplotlib.animation as animation
from pathlib import Path
import xarray as xr
import pandas as pd
import warnings
from typing import Tuple, List, Dict, Optional, Union
import os


def plot_full_mesh(nodes, elements, font_size=16, savefile=None, figsize=(10, 8), 
                   color="gray", show=True):
    """
    Plot a full ADCIRC mesh.
    
    Args:
        nodes (pd.DataFrame): DataFrame with at least 'x' and 'y' columns
        elements (pd.DataFrame): DataFrame with at least 'node_1', 'node_2', 'node_3' columns
        font_size (int, optional): Font size for labels. Defaults to 16.
        savefile (str, optional): Path to save the figure. Defaults to None.
        figsize (tuple, optional): Figure size. Defaults to (10, 8).
        color (str, optional): Color for the mesh lines. Defaults to "gray".
        show (bool, optional): Whether to show the plot. Defaults to True.
        
    Returns:
        tuple: Figure and axis objects
    """
    # Set global font size
    plt.rcParams.update({"font.size": font_size})

    # Extract data
    x = nodes["x"].values
    y = nodes["y"].values
    triangles = elements[["node_1", "node_2", "node_3"]].values - 1

    # Build figure & plot mesh
    fig, ax = plt.subplots(figsize=figsize)
    triang = mtri.Triangulation(x, y, triangles)
    ax.triplot(triang, color=color, linewidth=0.5, alpha=0.5, label="Mesh")

    # Labels, limits, ticks
    ax.set_title("ADCIRC Mesh", fontsize=font_size)
    ax.set_xlabel("Longitude", fontsize=font_size)
    ax.set_ylabel("Latitude", fontsize=font_size)
    ax.set_aspect("equal", "box")
    ax.tick_params(axis="both", which="major", labelsize=font_size)

    # Legend
    legend = ax.legend(loc="upper left", fontsize=font_size, frameon=False)
    legend_lines = legend.get_lines()
    legend_lines[0].set_linewidth(1.0)

    plt.tight_layout() 
    if savefile:
        plt.savefig(savefile, dpi=300)
    if show:
        plt.show()  
    return fig, ax


def plot_subgrid_overlay(
    full_mesh,
    sub_mesh,
    global_color="gray",
    subgrid_color="blue",
    font_size=16,
    zoom_region=None,
    savefile=None,
    figsize=(10, 8),
    show=True
):
    """
    Plot a full mesh with a submesh overlay.
    
    Args:
        full_mesh: Full mesh object with nodes and elements attributes
        sub_mesh: Submesh object with nodes and elements attributes
        global_color (str, optional): Color for the global mesh. Defaults to "gray".
        subgrid_color (str, optional): Color for the submesh. Defaults to "blue".
        font_size (int, optional): Font size for labels. Defaults to 16.
        zoom_region (List[float], optional): Region to zoom to [lon_min, lon_max, lat_min, lat_max].
        savefile (str, optional): Path to save the figure. Defaults to None.
        figsize (tuple, optional): Figure size. Defaults to (10, 8).
        show (bool, optional): Whether to show the plot. Defaults to True.
        
    Returns:
        tuple: Figure and axis objects
    """
    # Set global font size
    plt.rcParams.update({"font.size": font_size})

    # Extract data from full mesh
    nodes = full_mesh.nodes
    elements = full_mesh.elements
    x = nodes["x"].values
    y = nodes["y"].values
    triangles = elements[["node_1", "node_2", "node_3"]].values - 1

    # Build figure & plot meshes
    fig, ax = plt.subplots(figsize=figsize)
    global_triang = mtri.Triangulation(x, y, triangles)
    ax.triplot(global_triang, color=global_color, linewidth=0.5, alpha=0.5, label="Full Mesh")
    
    # Extract data from submesh
    sub_nodes = sub_mesh.nodes
    sub_x = sub_nodes["x"].values
    sub_y = sub_nodes["y"].values
    sub_triangles = sub_mesh.elements[["node_1", "node_2", "node_3"]].values - 1
    sub_triang = mtri.Triangulation(sub_x, sub_y, sub_triangles)
    ax.triplot(sub_triang, color=subgrid_color, linewidth=1.0, label="Submesh")

    # Labels, limits, ticks
    ax.set_title("Submesh", fontsize=font_size)
    ax.set_xlabel("Longitude", fontsize=font_size)
    ax.set_ylabel("Latitude", fontsize=font_size)
    
    if zoom_region is not None:
        # Set the limits to the zoom region
        lon_min, lon_max, lat_min, lat_max = zoom_region
        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
    else:
        # Set the limits to the full mesh bounds with a small margin
        lon_min, lon_max = nodes.x.min(), nodes.x.max()
        lat_min, lat_max = nodes.y.min(), nodes.y.max()
        margin = 0.05 * max(lon_max - lon_min, lat_max - lat_min)
        ax.set_xlim(lon_min - margin, lon_max + margin)
        ax.set_ylim(lat_min - margin, lat_max + margin)

    ax.set_aspect("equal", "box")
    ax.tick_params(axis="both", which="major", labelsize=font_size)

    # Legend with a fatter proxy for Global Mesh
    legend = ax.legend(loc="upper left", fontsize=font_size, frameon=False)
    legend_lines = legend.get_lines()
    legend_lines[0].set_linewidth(1.0)  # heavier stroke for global

    plt.tight_layout()
    
    if savefile:
        plt.savefig(savefile, dpi=300)
    
    if show:
        plt.show()
    
    return fig, ax


def plot_elevation(mesh, data, time_index=0, variable_name=None, title=None, 
                  cmap='viridis', mask_dry=True, dry_color='lightgray',
                  savefile=None, figsize=(10, 8), show=True, diagnostic=False):
    """
    Plot water surface elevation or other nodal data on a mesh, handling NaN values.
    
    Args:
        mesh: Mesh object with nodes and elements attributes
        data: Either xarray Dataset or numpy array with data values
        time_index (int, optional): Time index to plot. Defaults to 0.
        variable_name (str, optional): Variable name if data is xarray. Defaults to None.
        title (str, optional): Plot title. If None, auto-generated. Defaults to None.
        cmap (str, optional): Colormap name. Defaults to 'viridis'.
        mask_dry (bool, optional): Whether to mask dry areas (NaN values). Defaults to True.
        dry_color (str, optional): Color for dry areas. Defaults to 'lightgray'.
        savefile (str, optional): Path to save the figure. Defaults to None.
        figsize (tuple, optional): Figure size. Defaults to (10, 8).
        show (bool, optional): Whether to show the plot. Defaults to True.
        diagnostic (bool, optional): Whether to print diagnostic information. Defaults to False.
        
    Returns:
        tuple: Figure and axis objects, and a dictionary with plot data
    """
    # Extract mesh data
    x = mesh.nodes["x"].values
    y = mesh.nodes["y"].values
    triangles = mesh.elements[["node_1", "node_2", "node_3"]].values - 1
    
    # Extract data values
    if isinstance(data, xr.Dataset):
        if variable_name is None:
            # Try to find elevation variables
            elev_vars = [var for var in data.data_vars 
                        if var in ['zeta', 'eta', 'elev', 'elevation']]
            if not elev_vars:
                # Try to find any variable with time and node dimensions
                for var in data.data_vars:
                    if len(data[var].dims) == 2 and 'time' in data[var].dims and 'node' in data[var].dims:
                        elev_vars = [var]
                        if diagnostic:
                            print(f"Using '{var}' as variable.")
                        break
            
            if not elev_vars:
                raise ValueError("Could not identify a suitable variable in the dataset")
            
            variable_name = elev_vars[0]
        
        values = data[variable_name][time_index].values
        
        # Auto-generate title if not provided
        if title is None:
            time_str = ""
            if 'time' in data.coords:
                try:
                    time_val = data.time[time_index].values
                    if hasattr(time_val, 'item'):  # Convert numpy scalar to Python type
                        time_val = time_val.item()
                    time_str = f" at {time_val}"
                except:
                    time_str = f" at Time Step {time_index}"
            
            title = f"{variable_name.title()}{time_str}"
    else:
        # Assume data is a numpy array
        values = data
        
        # Auto-generate title if not provided
        if title is None:
            if variable_name:
                title = f"{variable_name.title()} at Time Step {time_index}"
            else:
                title = f"Data Visualization at Time Step {time_index}"
    
    # Check for dimension mismatch
    if len(x) != len(values):
        raise ValueError(f"Dimension mismatch: mesh has {len(x)} nodes, but data has {len(values)} points")
    
    # Create triangulation
    triang = mtri.Triangulation(x, y, triangles)
    
    # Check for NaN values
    nan_count = np.isnan(values).sum()
    
    # Print diagnostics if requested
    if diagnostic:
        print(f"Non-NaN values: {len(values) - nan_count} out of {len(values)} ({(len(values) - nan_count) / len(values):.1%})")
        if nan_count > 0:
            print(f"Found {nan_count} NaN values (likely dry areas)")
            print(f"Min value: {np.nanmin(values)}, Max value: {np.nanmax(values)}")
        
        # Count valid vs. NaN triangles
        valid_triangles = 0
        nan_triangles = 0
        for i in range(len(triangles)):
            triangle_vals = values[triangles[i]]
            if np.any(np.isnan(triangle_vals)):
                nan_triangles += 1
            else:
                valid_triangles += 1
        print(f"Valid triangles: {valid_triangles}, NaN triangles: {nan_triangles}")
        
        # Create histogram of values
        plt.figure(figsize=(8, 5))
        valid_values = values[~np.isnan(values)]
        if len(valid_values) > 0:
            plt.hist(valid_values, bins=30)
            plt.title(f"Histogram of {variable_name if variable_name else 'data'} (excluding NaNs)")
            plt.xlabel("Value")
            plt.ylabel("Count")
            if savefile:
                hist_file = savefile.replace(".png", "_histogram.png")
                plt.savefig(hist_file, dpi=300)
            if show:
                plt.show()
    
    # Handle NaN values for plotting
    result_data = {}
    
    if nan_count > 0 and mask_dry:
        # Create a mask for triangles with NaN nodes
        triangles_z = np.take(values, triangles)
        mask = np.any(np.isnan(triangles_z), axis=1)
        triang.set_mask(mask)
        
        # Custom color map with special color for NaN values
        custom_cmap = cm.get_cmap(cmap).copy()
        custom_cmap.set_bad(dry_color)
        
        # Create normalized levels
        min_val = np.nanmin(values)
        max_val = np.nanmax(values)
        levels = np.linspace(min_val, max_val, 21)
        norm = BoundaryNorm(levels, custom_cmap.N)
        
        result_data['levels'] = levels
        result_data['min_val'] = min_val
        result_data['max_val'] = max_val
    else:
        custom_cmap = cmap
        norm = None
    
    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)
    
    try:
        # Try contour plot first
        if nan_count > 0 and mask_dry:
            contour = ax.tricontourf(triang, values, cmap=custom_cmap, 
                                    norm=norm, levels=levels)
        else:
            contour = ax.tricontourf(triang, values, cmap=custom_cmap, levels=20)
        
        result_data['plot_type'] = 'contour'
        result_data['contour'] = contour
        
        # Add mesh lines
        ax.triplot(triang, 'k-', alpha=0.2, linewidth=0.5)
        
        # Add colorbar
        cbar = plt.colorbar(contour, ax=ax)
        if variable_name:
            cbar.set_label(variable_name)
        else:
            cbar.set_label('Value')
    except Exception as e:
        if diagnostic:
            print(f"Error in tricontourf: {str(e)}")
            print("Falling back to scatter plot...")
        
        # Fall back to scatter plot if contouring fails
        valid_mask = ~np.isnan(values)
        scatter = ax.scatter(x[valid_mask], y[valid_mask], c=values[valid_mask], 
                           cmap=cmap, s=25, edgecolor='none', alpha=0.8)
        
        result_data['plot_type'] = 'scatter'
        result_data['scatter'] = scatter
        
        # Add mesh lines
        ax.triplot(x, y, triangles, 'k-', alpha=0.2, linewidth=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        if variable_name:
            cbar.set_label(variable_name)
        else:
            cbar.set_label('Value')
        
        # Update title
        title = f"{title} (Scatter Plot)"
    
    # Add dry area note if needed
    if nan_count > 0 and mask_dry:
        if title:
            title = f"{title}\n(Dry areas in {dry_color})"
    
    # Set title and labels
    ax.set_title(title)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_aspect('equal')
    
    # Save if requested
    if savefile:
        plt.savefig(savefile, dpi=300)
    
    # Show if requested
    if show:
        plt.show()
    
    return fig, ax, result_data


def compare_elevation_data(original_ds, submesh_ds, mapping_df, node_indices=None, num_nodes=5):
    """
    Compare elevation data between original and submesh datasets for the same physical nodes in a time series plot.
    
    Args:
        original_ds: xarray Dataset from original fort.63
        submesh_ds: xarray Dataset from processed submesh fort.63
        mapping_df: DataFrame with node ID mapping (mesh_node_id -> subset_node_id)
        node_indices: List of specific node indices to compare (if None, random selection)
        num_nodes: Number of nodes to compare if node_indices is None
    """
    # Identify elevation variable
    elev_vars = [var for var in original_ds.data_vars if var in ['zeta', 'eta', 'elev', 'elevation']]
    if not elev_vars:
        for var in original_ds.data_vars:
            if len(original_ds[var].dims) == 2 and 'time' in original_ds[var].dims and 'node' in original_ds[var].dims:
                elev_vars = [var]
                break
    
    if not elev_vars:
        raise ValueError("Could not identify elevation variable")
    
    elev_var = elev_vars[0]
    
    # Handle time coordinates
    times_orig = original_ds.time.values
    times_sub = submesh_ds.time.values
    
    if np.issubdtype(times_orig.dtype, np.datetime64):
        time_hours_orig = np.array([(t - times_orig[0]) for t in times_orig], dtype='timedelta64[h]').astype(float)
        time_hours_sub = np.array([(t - times_sub[0]) for t in times_sub], dtype='timedelta64[h]').astype(float)
        time_label = 'Hours since start'
    else:
        time_hours_orig = times_orig * 24
        time_hours_sub = times_sub * 24
        time_label = 'Hours'
    
    # Select nodes to compare
    if node_indices is None:
        # Select random nodes that exist in both datasets
        available_indices = mapping_df.index[:min(len(mapping_df), 100)]  # First 100 for performance
        node_indices = np.random.choice(available_indices, min(num_nodes, len(available_indices)), replace=False)
    
    fig, axes = plt.subplots(len(node_indices), 1, figsize=(12, 3*len(node_indices)))
    if len(node_indices) == 1:
        axes = [axes]
    
    for i, orig_node_idx in enumerate(node_indices):
        # Get corresponding submesh node index
        sub_node_idx = mapping_df.loc[orig_node_idx, 'subset_node_id'] - 1  # Convert to 0-based
        
        # Extract data for both datasets
        orig_data = original_ds[elev_var][:, orig_node_idx].values
        sub_data = submesh_ds[elev_var][:, sub_node_idx].values
        
        # Plot comparison
        ax = axes[i]
        ax.plot(time_hours_orig, orig_data, 'b-', linewidth=2, label=f'Original (Node {orig_node_idx})', alpha=0.8)
        ax.plot(time_hours_sub, sub_data, 'r--', linewidth=2, label=f'Submesh (Node {sub_node_idx})', alpha=0.8)
        
        ax.set_xlabel(time_label)
        ax.set_ylabel('Water Surface Elevation (m)')
        ax.set_title(f'Elevation Comparison: Original Node {orig_node_idx} vs Submesh Node {sub_node_idx}')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        rmse = np.sqrt(np.mean((orig_data - sub_data)**2))
        max_diff = np.max(np.abs(orig_data - sub_data))
        ax.text(0.02, 0.95, f'RMSE: {rmse:.4f} m\nMax Diff: {max_diff:.4f} m', 
                transform=ax.transAxes, bbox=dict(boxstyle='round', alpha=0.1))
    
    plt.tight_layout()
    plt.savefig('elevation_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return fig, axes


def compare_velocity_data(original_ds, submesh_ds, mapping_df, node_indices=None, num_nodes=3):
    """
    Compare velocity data between original and submesh datasets for the same physical nodes.
    
    Args:
        original_ds: xarray Dataset from original fort.64
        submesh_ds: xarray Dataset from processed submesh fort.64
        mapping_df: DataFrame with node ID mapping
        node_indices: List of specific node indices to compare
        num_nodes: Number of nodes to compare if node_indices is None
    """
    # Identify velocity variables
    uvel_vars = [var for var in original_ds.data_vars if var in ['u_vel', 'u-vel', 'uvel', 'x_vel']]
    vvel_vars = [var for var in original_ds.data_vars if var in ['v_vel', 'v-vel', 'vvel', 'y_vel']]
    
    if not uvel_vars or not vvel_vars:
        raise ValueError("Could not identify velocity variables")
    
    uvel_var = uvel_vars[0]
    vvel_var = vvel_vars[0]
    
    # Handle time coordinates
    times_orig = original_ds.time.values
    times_sub = submesh_ds.time.values
    
    if np.issubdtype(times_orig.dtype, np.datetime64):
        time_hours_orig = np.array([(t - times_orig[0]) for t in times_orig], dtype='timedelta64[h]').astype(float)
        time_hours_sub = np.array([(t - times_sub[0]) for t in times_sub], dtype='timedelta64[h]').astype(float)
        time_label = 'Hours since start'
    else:
        time_hours_orig = times_orig * 24
        time_hours_sub = times_sub * 24
        time_label = 'Hours'
    
    # Select nodes to compare
    if node_indices is None:
        available_indices = mapping_df.index[:min(len(mapping_df), 50)]
        node_indices = np.random.choice(available_indices, min(num_nodes, len(available_indices)), replace=False)
    
    fig, axes = plt.subplots(len(node_indices), 3, figsize=(18, 4*len(node_indices)))
    if len(node_indices) == 1:
        axes = axes.reshape(1, -1)
    
    for i, orig_node_idx in enumerate(node_indices):
        # Get corresponding submesh node index
        sub_node_idx = mapping_df.loc[orig_node_idx, 'subset_node_id'] - 1  # Convert to 0-based
        
        # Extract velocity data
        u_orig = original_ds[uvel_var][:, orig_node_idx].values
        v_orig = original_ds[vvel_var][:, orig_node_idx].values
        u_sub = submesh_ds[uvel_var][:, sub_node_idx].values
        v_sub = submesh_ds[vvel_var][:, sub_node_idx].values
        
        # Calculate magnitudes
        mag_orig = np.sqrt(u_orig**2 + v_orig**2)
        mag_sub = np.sqrt(u_sub**2 + v_sub**2)
        
        # Plot U-velocity
        axes[i, 0].plot(time_hours_orig, u_orig, 'b-', linewidth=2, label=f'Original (Node {orig_node_idx})', alpha=0.8)
        axes[i, 0].plot(time_hours_sub, u_sub, 'r--', linewidth=2, label=f'Submesh (Node {sub_node_idx})', alpha=0.8)
        axes[i, 0].set_xlabel(time_label)
        axes[i, 0].set_ylabel('U-Velocity (m/s)')
        axes[i, 0].set_title(f'U-Velocity: Node {orig_node_idx} vs {sub_node_idx}')
        axes[i, 0].legend()
        axes[i, 0].grid(True, alpha=0.3)
        
        # Plot V-velocity
        axes[i, 1].plot(time_hours_orig, v_orig, 'b-', linewidth=2, label=f'Original (Node {orig_node_idx})', alpha=0.8)
        axes[i, 1].plot(time_hours_sub, v_sub, 'r--', linewidth=2, label=f'Submesh (Node {sub_node_idx})', alpha=0.8)
        axes[i, 1].set_xlabel(time_label)
        axes[i, 1].set_ylabel('V-Velocity (m/s)')
        axes[i, 1].set_title(f'V-Velocity: Node {orig_node_idx} vs {sub_node_idx}')
        axes[i, 1].legend()
        axes[i, 1].grid(True, alpha=0.3)
        
        # Plot velocity magnitude
        axes[i, 2].plot(time_hours_orig, mag_orig, 'b-', linewidth=2, label=f'Original (Node {orig_node_idx})', alpha=0.8)
        axes[i, 2].plot(time_hours_sub, mag_sub, 'r--', linewidth=2, label=f'Submesh (Node {sub_node_idx})', alpha=0.8)
        axes[i, 2].set_xlabel(time_label)
        axes[i, 2].set_ylabel('Velocity Magnitude (m/s)')
        axes[i, 2].set_title(f'Velocity Magnitude: Node {orig_node_idx} vs {sub_node_idx}')
        axes[i, 2].legend()
        axes[i, 2].grid(True, alpha=0.3)
        
        # Add statistics to magnitude plot
        rmse_mag = np.sqrt(np.mean((mag_orig - mag_sub)**2))
        max_diff_mag = np.max(np.abs(mag_orig - mag_sub))
        axes[i, 2].text(0.02, 0.95, f'RMSE: {rmse_mag:.4f} m/s\nMax Diff: {max_diff_mag:.4f} m/s', 
                       transform=axes[i, 2].transAxes, bbox=dict(boxstyle='round', alpha=0.1))
    
    plt.tight_layout()
    plt.savefig('velocity_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return fig, axes


def create_elevation_animation(submesh, ds, output_filename='elevation_animation.gif', 
                             elev_var=None, fps=2, dpi=100, figsize=(10, 8)):
    """
    Create an animated GIF of water surface elevation over time for the submesh.
    
    Args:
        submesh: Mesh object with nodes and elements
        ds: xarray Dataset with elevation data
        output_filename: Name of output GIF file
        elev_var: Name of elevation variable (auto-detected if None)
        fps: Frames per second for the animation
        dpi: Resolution of the animation
        figsize: Figure size (width, height)
    """
    # Identify elevation variable
    if elev_var is None:
        elev_vars = [var for var in ds.data_vars if var in ['zeta', 'eta', 'elev', 'elevation']]
        if not elev_vars:
            for var in ds.data_vars:
                if len(ds[var].dims) == 2 and 'time' in ds[var].dims and 'node' in ds[var].dims:
                    elev_vars = [var]
                    break
        
        if not elev_vars:
            raise ValueError("Could not identify elevation variable")
        
        elev_var = elev_vars[0]
    
    # Extract mesh geometry
    x = submesh.nodes["x"].values
    y = submesh.nodes["y"].values
    triangles = submesh.elements[["node_1", "node_2", "node_3"]].values - 1
    
    # Create triangulation
    triang = mtri.Triangulation(x, y, triangles)
    
    # Get all elevation data
    elevation_data = ds[elev_var].values
    times = ds.time.values
    
    # Handle time formatting
    if np.issubdtype(times.dtype, np.datetime64):
        time_labels = [f"Time: {pd.to_datetime(t).strftime('%Y-%m-%d %H:%M')}" for t in times]
    else:
        time_labels = [f"Time: {t:.2f} hours" for t in times * 24]
    
    # Calculate global min/max for consistent color scale
    global_min = np.nanmin(elevation_data)
    global_max = np.nanmax(elevation_data)
    
    # Create colormap and normalization
    cmap = cm.get_cmap('viridis').copy()
    cmap.set_bad('lightgray')
    levels = np.linspace(global_min, global_max, 21)
    norm = BoundaryNorm(levels, cmap.N)
    
    # Set up the figure and axis
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(x.min() - 0.01, x.max() + 0.01)
    ax.set_ylim(y.min() - 0.01, y.max() + 0.01)
    ax.set_aspect('equal')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    
    # Initialize the plot elements
    def init():
        ax.clear()
        ax.set_xlim(x.min() - 0.01, x.max() + 0.01)
        ax.set_ylim(y.min() - 0.01, y.max() + 0.01)
        ax.set_aspect('equal')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        return []
    
    def animate(frame):
        ax.clear()
        ax.set_xlim(x.min() - 0.01, x.max() + 0.01)
        ax.set_ylim(y.min() - 0.01, y.max() + 0.01)
        ax.set_aspect('equal')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        
        # Get data for this time step
        elev_frame = elevation_data[frame]
        
        # Handle NaN values by masking triangles
        triangles_z = np.take(elev_frame, triangles)
        mask = np.any(np.isnan(triangles_z), axis=1)
        triang_masked = mtri.Triangulation(x, y, triangles)
        triang_masked.set_mask(mask)
        
        try:
            # Plot the elevation data
            contour = ax.tricontourf(triang_masked, elev_frame, cmap=cmap, norm=norm, levels=levels)
            
            # Add mesh lines
            ax.triplot(triang_masked, 'k-', alpha=0.1, linewidth=0.3)
            
        except Exception as e:
            print(f"Warning: Contour plot failed for frame {frame}, using scatter plot: {e}")
            # Fallback to scatter plot
            valid_mask = ~np.isnan(elev_frame)
            scatter = ax.scatter(x[valid_mask], y[valid_mask], c=elev_frame[valid_mask], 
                               cmap=cmap, norm=norm, s=10, edgecolor='none')
        
        # Set title with time information
        ax.set_title(f'Water Surface Elevation\n{time_labels[frame]}')
        
        return []
    
    # Create animation
    print(f"Creating animation with {len(times)} frames...")
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(times), 
                                 interval=1000//fps, blit=False, repeat=True)
    
    # Add colorbar (do this outside the animation loop)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8)
    cbar.set_label('Water Surface Elevation (m)')
    
    # Save as GIF
    print(f"Saving animation to {output_filename}...")
    try:
        # Try to use pillow writer first (better compression)
        writer = animation.PillowWriter(fps=fps)
        anim.save(output_filename, writer=writer, dpi=dpi)
    except Exception as e:
        print(f"PillowWriter failed ({e}), trying ImageMagickWriter...")
        try:
            writer = animation.ImageMagickWriter(fps=fps)
            anim.save(output_filename, writer=writer, dpi=dpi)
        except Exception as e2:
            print(f"ImageMagickWriter also failed ({e2}), saving individual frames...")
            # Save individual frames as fallback
            for i in range(len(times)):
                animate(i)
                plt.savefig(f'frame_{i:03d}.png', dpi=dpi, bbox_inches='tight')
            print(f"Saved {len(times)} individual frames as frame_XXX.png")
    
    plt.close(fig)
    print(f"Animation saved successfully!")
    
    return anim



