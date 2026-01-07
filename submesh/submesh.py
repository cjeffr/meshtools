"""
ADCIRC Mesh Subsetting Module
This module p[rovides tools for extract subsets of ADCIRC mesh files, typically used in coastal ocean modeling.
It allows users to create smaller meeshes based on goegraphic constraints like rectangles, circles or arbitrary polygons.
The module includes classes for different mask generators, a Mesh class for handling mesh data, and functions for
reading/writing mesh files in the ADCIRC format.
The module also provides functionality for subsetting output data from ADCIRC simulations based on the generated masks.

Key components:
- MaskGenerator: Abstract base class for creating masks to select subsets of mesh nodes.
- Mesh class for manipulating mesh data (fort.14 format)
- Functions for processing and subsetting related output files (fort.63, fort.64)

Typical usage:
```python
# Load a full mesh
mesh = Mesh.from_fort14("path/to/fort.14")
# Create a mask for a rectangular area
mesh_subset  = mesh.subset(RectangleMask(x1, y1, x2, y2))
# Save the subset mesh to a new fort.14 file
mesh_subset.to_fort14("path/to/subset.fort.14")
```"""

import sys
import os
from pathlib import Path
from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import dataclass, field
from typing import List, Dict, Optional
import pandas as pd
import numpy as np
import xarray as xr

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from .io import read_fort14 as sio


class MaskGenerator(ABC):
    """Abstract base for any 'mask' that selects a subset of mesh nodes.
    Subclasses implement make_mask(nodes)-> pd.Series['bool'].
    """

    @abstractmethod
    def make_mask(self, nodes: pd.DataFrame) -> pd.Series:
        """Given the full 'nodes' DataFrame with at least 'x' and 'y' columns,
        return a boolean Series indexed indentically, where True marks nodes
        to keep.

        Args:
            nodes (pd.DataFrame): Dataframe of mesh nodes, 'x', 'y', 'z' with index as full mesh
            id numbers

        Raises:
            ValueError: _description_

        Returns:
            pd.Series: Boolean mask for mesh nodes included in subset of original mesh.
        """


class RectangleMask(MaskGenerator):
    """Create a boolean mask based on corners of a rectangle for a mesh.

    Args:
        MaskGenerator (_type_): _description_
    """

    def __init__(self, x1, y1, x2, y2):
        """
        Extracts a rectangular subgrid from a fort.14 mesh file.

        Args:
            fort14_path (str): Path to the fort.14 file.
            x1, y1, x2, y2 (float): Coordinates defining the rectangular region.

        Retuns:
            subgrid_nodes (DataFrame): renumbered subgrid nodes
            subgrid_elements (DataFrame): Elements using the new node IDs
            node_id_mapping (dict): Mapping from original node IDs -> new subgrid Ids
        """

        p1 = np.array([x1, y1])
        p2 = np.array([x2, y2])
        self.center = (p1 + p2) / 2.0
        self.half_width = np.abs(p1 - p2) / 2.0

    def make_mask(self, nodes: pd.DataFrame) -> pd.Series:
        """
        nodes must have 'x' and 'y' columns.
        This method checks if the nodes are within the rectangle defined by
        (x1, y1) and (x2, y2).
        It returns a boolean Series where True indicates that the node is
        within the rectangle.
        """

        return nodes[["x", "y"]].sub(self.center).abs().le(self.half_width).all(axis=1)


class ArbitrMask(MaskGenerator):
    def __init__(self, mesh_filename: str, eps: float = 1.0e-12):
        """
        A MaskGenerator that reads a second mesh (SMS/Fort.14) of arbitrary triangle shapes
        and then, for any target mesh, returns True for nodes that lie inside any of those triangles.
        """

        data = sio(mesh_filename)
        # nodes: DataFrame with at least ('node', 'x', 'y')
        # elements: DataFrame with at least ('node_1', 'node_2', 'node_3')
        poly_nodes = data["nodes"]
        poly_elements = data["elements"]

        # Convert to numpy for speed and compute 0-based connectivity
        # shape xx: (2, n_poly_nodes)
        self.xx = poly_nodes[["x", "y"]].to_numpy().T
        # shape node_conn: (3, n_poly_elements), indices 0...n_nodes-1
        conn = poly_elements[["node_1", "node_2", "node_3"]].to_numpy()
        self.node_conn = (conn - 1).T.astype(int)

        # Precompute the bounding box of the polygon mesh
        xmax = self.xx.max(axis=1)
        xmin = self.xx.min(axis=1)
        self.xa = 0.5 * (xmax + xmin)
        self.aa = 0.5 * np.abs(xmax - xmin)
        self.eps = eps

    def make_mask(self, nodes: pd.DataFrame) -> pd.Series:
        """For *each* row in 'nodes' (must have 'x' and 'y'), returns True
        if that point lies within the SMS triangles loaded at init
        """
        coords = nodes[["x", "y"]].to_numpy()
        xc, yc = coords[:, 0], coords[:, 1]
        n_pts = len(nodes)

        # First filter by bounding box
        outside = (np.abs(xc - self.xa[0]) > self.aa[0]) & (
            np.abs(yc - self.xa[1]) > self.aa[1]
        )
        to_test = np.where(~outside)[0]

        if len(to_test) == 0:
            # No points to test, return all False
            return pd.Series(np.zeros(n_pts, dtype=bool), index=nodes.index)

        # Extract triangle vertices
        i1, i2, i3 = self.node_conn
        x1, y1 = self.xx[0, i1], self.xx[1, i1]
        x2, y2 = self.xx[0, i2], self.xx[1, i2]
        x3, y3 = self.xx[0, i3], self.xx[1, i3]

        # Compute vectors for barycentric coordinates
        # Precompute denominator components
        a0 = (x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)

        # Create mask to avoid division by zero
        valid_triangles = np.abs(a0) > self.eps
        a0_safe = np.where(valid_triangles, a0, 1.0)  # Avoid division by zero

        inside_any = np.zeros(n_pts, dtype=bool)
        batch_size = 1000

        for batch_start in range(0, len(to_test), batch_size):
            batch_end = min(batch_start + batch_size, len(to_test))
            batch_indices = to_test[batch_start:batch_end]
            n_test = len(batch_indices)

            # Reshape for broadcasting
            batch_size_actual = len(batch_indices)
            batch_x = xc[batch_indices]
            batch_y = yc[batch_indices]

            x_pts = batch_x[:, np.newaxis]
            y_pts = batch_y[:, np.newaxis]

            # Compute barycentric coordinates for all points against all triangles
            a1 = ((x_pts - x2) * (y_pts - y3) - (x_pts - x3) * (y_pts - y2)) / a0_safe
            a2 = ((x1 - x_pts) * (y1 - y3) - (x1 - x3) * (y1 - y_pts)) / a0_safe
            a3 = ((x1 - x2) * (y1 - y_pts) - (x1 - x_pts) * (y1 - y2)) / a0_safe
            # A point is inside a triangle if all barycentric coordinates are >= 0
            is_inside = (
                (a1 >= -self.eps)
                & (a2 >= -self.eps)
                & (a3 >= -self.eps)
                & valid_triangles
            )

            inside_batch = np.any(is_inside, axis=1)
            inside_any[batch_indices] = inside_batch

        return pd.Series(inside_any, index=nodes.index)


class CircleMask(MaskGenerator):
    """Extracts a mask for a circular subgrid from a fort.14 mesh file.
    Uses the center of the circle in lat/lon coordinates and the diameter
    of the circle in decimal degrees.
    """

    def __init__(self, center: tuple, diameter: float):
        """
        Args:
            center (tuple): Center of the circle in lat/lon coordinates.
            diameter (float): Diameter of the circle in decimal degrees.
        """
        self.cx = center[0]
        self.cy = center[1]
        self.diameter = diameter
        self.radius = diameter / 2.0

    def make_mask(self, nodes: pd.DataFrame) -> pd.Series:
        """nodes must have 'x' and 'y' columns in the same units decimal degrees
        returns a boolean Series where True indicatres that the node lies within the circle
        """
        dx = nodes["x"] - self.cx
        dy = nodes["y"] - self.cy
        # Calculate the distance from the center of the circle
        inside = dx**2 + dy**2 <= (self.radius) ** 2
        return pd.Series(inside.values, index=nodes.index)


class ShapefileMask(MaskGenerator):
    """A maskgenerator that reads one (or more) polygons from a shapefile
    and then, for any target mesh, returns True for nodes that lie inside
    any of those polygons.
    """

    def __init__(self, shapefile: str, eps: float = 1.0e-12):
        """
        A MaskGenerator that reads a shapefile of arbitrary polygon shapes
        and then, for any target mesh, returns True for nodes that lie inside any of those polygons.
        """
        import geopandas as gpd
        from shapely.geometry import Polygon

        self.gdf = gpd.read_file(shapefile)
        # 2) Reproject into geographic (lon/lat)
        if self.gdf.crs != "EPSG:4326":
            self.gdf = self.gdf.to_crs(epsg=4326)
        if not all(self.gdf.geometry.geom_type.isin(["Polygon", "MultiPolygon"])):
            raise ValueError("Only Polygon and MultiPolygon geometries are supported.")

        self.mask_geom = self.gdf.unary_union

    def make_mask(self, nodes: pd.DataFrame) -> pd.Series:
        """For *each* row in 'nodes' (must have 'x' and 'y'), returns True
        if that point lies within the polygons loaded at init
        """
        from shapely.geometry import Point

        coords = nodes[["x", "y"]].to_numpy()
        points = [Point(x, y) for x, y in coords]
        mask = [self.mask_geom.contains(point) for point in points]
        return pd.Series(mask, index=nodes.index)


_MASK_REGISTRY = {
    "rectangle": RectangleMask,
    "arbitrary": ArbitrMask,
    "shapefile": ShapefileMask,
    "circle": CircleMask,
}


def get_mask_generator(mask_type: str, *args, **kwargs) -> MaskGenerator:
    """
    Factory function to create a mask generator based on the specified type.

    Args:
        mask_type (str): Type of mask to create ('rectangle', 'arbitrary').
        **params: Additional parameters for the mask generator.
    """
    cls = _MASK_REGISTRY.get(mask_type)
    if cls is None:
        raise ValueError(
            f"Unknown mask type: {mask_type}. Available types: {list(_MASK_REGISTRY.keys())}"
        )
    return cls(*args, **kwargs)


@dataclass
class Mesh:
    nodes: pd.DataFrame
    elements: pd.DataFrame
    boundaries: Dict[str, List[Dict[str, any]]]
    node_props: Optional[pd.Series] = None
    description: str = "ADCIRC mesh"

    @staticmethod
    def from_fort14(filepath: str) -> "Mesh":
        """Read a fort.14 file and return a Mesh object."""
        data = sio(filepath)
        return Mesh(
            nodes=data["nodes"],
            elements=data["elements"],
            boundaries=data["boundaries"],
            node_props=data.get("node_props"),
            description=data.get("description", "ADCIRC Mesh"),
        )

    def subset(self, mask: MaskGenerator | str, folder: str = '.', **mask_params) -> "Mesh":
        """
        Creates a subset of the mesh based on coordinates of a rectangle with corners
        (x1,y1), (x2,y2)
        Args:
            mask (pd.Series): mask of all nodes that are contained within the rectangle

        Returns:
            Mesh: subset mesh dataclass
        """
        if isinstance(mask, str):
            mask = get_mask_generator(mask, **mask_params)
        bool_series = mask.make_mask(self.nodes)
        self.nodes["subgrid_mask"] = bool_series
        return self._build_submesh(folder=folder)

    def _build_submesh(self, folder: str = '.') -> "Mesh":
        """Submesh creation from the original mesh with a boolean mask of nodes to not include

        Returns:
            Mesh: submesh
        """
        self.expand_element_mask()

        # Extract subset of nodes
        nodes_sub = self.nodes[self.nodes["subgrid_mask"]].copy()
        nodes_sub["new_index"] = np.arange(1, len(nodes_sub) + 1)

        # ✅ Create Series-based mapping (faster for .map())
        node_id_map = pd.Series(
            nodes_sub.new_index.values,
            index=nodes_sub.index.values,
            name="subset_node_id",
        )
        Path(folder).mkdir(parents=True, exist_ok=True)
        node_id_map.to_csv(Path(folder) / "node_id_map.csv", index_label="mesh_node_id")

        # Update elements
        node_cols = ["node_1", "node_2", "node_3"]
        elements_sub = self.elements[self.elements["subgrid_mask"]].copy()
        # ✅ Use vectorized .map instead of applymap
        elements_sub[node_cols] = elements_sub[node_cols].apply(
            lambda col: col.map(node_id_map)
        )
        # SAVE ELEMENT MASK TOO!
        element_mask = self.elements["subgrid_mask"].copy()
        element_mask.to_csv(Path(folder) / "element_mask.csv", index_label="mesh_element_id", header=["in_subgrid"])

        # Continue as before
        nprops_sub = self.compute_nprops()
        boundaries_sub = self.parse_subgrid_boundaries(node_id_map)

        return Mesh(
            nodes=nodes_sub,
            elements=elements_sub,
            boundaries=boundaries_sub,
            description="Subgrid",
        )

    def remap_node_ids(self) -> Dict[int, int]:
        """Return a mapping from original node ID to new 1-based_index."""
        return dict(zip(self.nodes["id"], self.nodes.index))

    def expand_element_mask(self) -> None:
        """Expand the node mask to include all nodes from selected elements, and update the element mask."""
        node_ids_in_mask = self.nodes.index[self.nodes["subgrid_mask"]]
        node_cols = ["node_1", "node_2", "node_3"]

        self.elements["subgrid_mask"] = (
            self.elements[node_cols].isin(node_ids_in_mask).any(axis=1)
        )
        selected_elements = self.elements[self.elements["subgrid_mask"]]
        expanded_node_ids = pd.unique(selected_elements[node_cols].values.ravel())
        self.nodes["subgrid_mask"] = self.nodes.index.isin(expanded_node_ids)

    def compute_nprops(self) -> pd.Series:
        """
        Compute nodal property array (nprops) used in ADCIRC style meshes

        This method preserves the global node indexing and computes the sum of
        directional 'delta' values foe each node across all elements.

        Args:
            df_elements (DataFrame): Element table with 'node_1', 'node_2', 'node_3'
            num_nodes (int): Total number of nodes in the global mesh

        Returns:
            Series: nprops values indexed by global node ID (1-based)
        """
        node_1 = self.elements["node_1"]
        node_2 = self.elements["node_2"]
        node_3 = self.elements["node_3"]

        delta_1 = node_2 - node_3
        delta_2 = node_3 - node_1
        delta_3 = node_1 - node_2

        updates = pd.concat(
            [
                pd.DataFrame({"node": node_1, "delta": delta_1}),
                pd.DataFrame({"node": node_2, "delta": delta_2}),
                pd.DataFrame({"node": node_3, "delta": delta_3}),
            ]
        )

        return (
            updates.groupby("node", sort=False)["delta"]
            .sum()
            .reindex(self.nodes.index, fill_value=0)
        )

    def count_boundary_nodes(self) -> pd.Series:
        """
        Vectorized counting of boundary nodes based on the element connectivity table.

        Args:
            df_elements (pd.DataFrame): DataFrame with columns ['node1', 'node2', 'node3'].

        Returns:
            pd.Series: A series where index is node_id and values are the count of occurrences.
        """
        # Step 1: Convert from wide format (3 columns) to long format (single node column)
        node_series = self.elements.melt(value_name="node_id", var_name="node_column")[
            ["node_id"]
        ]
        # Step 2: Count occurrences of each node
        return node_series["node_id"].value_counts()

    def get_element_edges(self) -> np.ndarray:
        """
        Extracts and counts edges from a triangle mesh, then identifies boundary edges.

        Args:
            elements_df (DataFrame): Triangle elemtsn with 'node_1', 'node_2', 'node_3'

        Returns:
            nbl (int): Number boundary edges
            nblnc (ndarray): 2xnbl array of boundary edge node pairs
        """

        # Extract the node columns into separate arrays
        n1 = self.elements["node_1"].to_numpy()
        n2 = self.elements["node_2"].to_numpy()
        n3 = self.elements["node_3"].to_numpy()

        # Build all edges as (min, max) tuples
        edges = np.vstack(
            [
                np.stack([np.minimum(n1, n2), np.maximum(n1, n2)], axis=1),
                np.stack([np.minimum(n2, n3), np.maximum(n2, n3)], axis=1),
                np.stack([np.minimum(n3, n1), np.maximum(n3, n1)], axis=1),
            ]
        )

        # Convert to list of tuples for Counter
        edge_tuples = list(map(tuple, edges))

        # Count and filter
        edge_counts = Counter(edge_tuples)
        boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
        return boundary_edges

    def remap_boundary_edges(self, node_id_map):
        """
        From element connectivity, extract unique boundary edges and convert them to subgrid node indices.

        Args:
            node_id_map (pd.Series): Pandas Series mapping original node ID ➜ new subgrid index

        Returns:
            boundary_segments (List[List[int]]): List of node pairs defining boundary edges
            boundary_lengths (List[int]): Lengths of each boundary segment
            mapped_edges (np.ndarray): 2 x n array of subgrid indices for boundary edges
        """
        # Get original boundary edges
        edge_array = self.get_element_edges()  # List of (n1, n2) tuples

        # Convert to 2 x N array for mapping
        edge_df = pd.DataFrame(edge_array, columns=["n1", "n2"])

        # ✅ Use .map on each column
        mapped = (
            edge_df.apply(lambda col: col.map(node_id_map), axis=0).dropna().astype(int)
        )

        # Transpose into 2 x nbl shape for compatibility
        mapped_edges = mapped.to_numpy().T

        # Optional: build boundary segments as 2-node lists
        boundary_segments = mapped.values.tolist()
        boundary_lengths = [2] * len(boundary_segments)

        return boundary_segments, boundary_lengths, mapped_edges

    def parse_subgrid_boundaries(self, node_id_map: Dict[int, int]) -> Dict:
        """Parse and rebuild both land and open boundary structures for a subgrid mesh.

        This method remaps all boundary segments using the provided node ID mapping.
        Boundary segments that have no remaining nodes after remapping are discarded.

        Args:
            node_id_map: Mapping from original node IDs to new submesh IDs

        Returns:
            Dictionary with 'land' and 'open' keys containing lists of BoundarySegment objects
        """
        bound_segments, bound_lengths, mapped_edges = self.remap_boundary_edges(
            node_id_map
        )

        # Remap all boundaries using the new unified approach
        remapped_boundaries = {'open': [], 'land': []}

        for category in ['open', 'land']:
            segments = self.boundaries.get(category, [])

            for segment in segments:
                remapped_segment = segment.remap(node_id_map)

                # Only keep segments that have at least one node remaining
                if remapped_segment is not None:
                    remapped_boundaries[category].append(remapped_segment)

        return remapped_boundaries


    def to_fort14(self, filepath: str):
        """Write the mesh to an ADCIRC fort.14 file.

        This method writes the current mesh to a fort.14 file that can be used
        by ADCIRC or other compatible software.

        Args:
            filepath (str): Path to the output fort.14 file

        Note:
            The file format follows the ADCIRC conventions with specific
            formatting requirements. The output includes node coordinates,
            element connectivity, and boundary information.

        Example:
            ```python
            mesh.to_fort14("output_mesh.fort.14")
            ```
        """
        lines = []
        # Write the header
        lines.append("grid")
        # elements and nodes
        lines.append(f"{len(self.elements)} {len(self.nodes)}")
        # Write node coorindates (ID, x, y, z)
        for _, row in self.nodes.iterrows():
            lines.append(
                f"{row['new_index']} {row['x']:18.10f} {row['y']:18.10f} {row['z']:18.10E}"
            )
        # Write element connectivity (ID, num_nodes_per_elements, node_1, node_2, node_3)
        for i, (_, row) in enumerate(self.elements.iterrows(), start=1):
            lines.append(
                f"{i:10d} {3:10d} {row['node_1']:10d} {row['node_2']:10d} {row['node_3']:10d}"
            )

        # Open boundaries
        open_segments = self.boundaries.get("open", [])
        lines.append(f"{len(open_segments):10d}")

        # Total number of open boundary nodes
        total_nope = sum(len(seg.nodes) for seg in open_segments)
        lines.append(f"{total_nope:10d}")

        for segment in open_segments:
            lines.append(f"{len(segment.nodes):10d}")
            for node in segment.nodes:
                lines.append(f"{node.node_id:10d}")

        # Land boundaries
        land_segments = self.boundaries.get("land", [])
        lines.append(f"{len(land_segments):10d}")

        # Total number of land boundary nodes
        total_nvell = sum(len(seg.nodes) for seg in land_segments)
        lines.append(f"{total_nvell:10d}")

        for segment in land_segments:
            ibtype = segment.ibtype
            lines.append(f"{len(segment.nodes):10d} {ibtype:10d}")

            # Write each node based on its attributes
            for node in segment.nodes:
                if node.is_simple():
                    lines.append(f"{node.node_id:10d}")

                elif node.is_barrier():
                    h1, h2 = node.barrier_heights
                    lines.append(f"{node.node_id:10d} {h1:18.10f} {h2:18.10f}")

                elif node.is_connected():
                    v1, v2, v3 = node.barrier_heights
                    lines.append(f"{node.node_id:10d} {node.connected_to:10d} {v1:18.10f} {v2:18.10f} {v3:18.10f}")

        # Write all lines in one go, no trailing newline
        with open(filepath, "w") as f:
            f.write("\n".join(lines))

        print(f"✅ fort.14 written to {filepath}")


def extract_fort_variables(nc_path, variables=None):
    ds = xr.open_dataset(nc_path)
    ds_trimmed = ds[variables]
    ds_trimmed.to_netcdf(
        "tmp_output.nc", encoding={v: {"zlib": False} for v in ds_trimmed.data_vars}
    )


def subset_output(mapping_df, output_file, savefile=None):
    """Subset ADCIRC output data based on a node mapping.

    This function automatically determines the input file type (.nc for NetCDF,
    .63/.64 for ASCII) and processes it accordingly, subsetting the data based
    on the provided node mapping and saving the result as a NetCDF file.

    Args:
        mapping_df (str): Path to CSV file containing node ID mapping
        output_file (str): Path to the input file (NetCDF or ASCII)
        savefile (str, optional): Path to save the subsetted NetCDF file.
                                 If None, a name will be generated based on the input file.

    Returns:
        str: Path to the saved output file

    Raises:
        ValueError: If mapping_df is None
        ValueError: If the file extension is not recognized

    Note:
        This function creates a node-based subset of the data and preserves
        all other dimensions and variables. For NetCDF inputs, the node ordering
        from the mapping file will be maintained in the output.
    """
    from .io import output_to_xarray, save_netcdf

    # Validate input
    if mapping_df is None:
        raise ValueError("Subgrid mapping is required")

    # Read the node ID mapping
    df_map = pd.read_csv(mapping_df, engine="python")
    mesh_ids = df_map["mesh_node_id"].values

    # Determine input file type and process accordingly
    file_ext = os.path.splitext(output_file)[1].lower()

    # Generate default output filename if none provided
    if savefile is None:
        base_name = os.path.splitext(os.path.basename(output_file))[0]
        savefile = f"{base_name}_subset.nc"

    # Process based on file extension
    elif file_ext == ".nc":
        print(f"📊 Processing NetCDF file: {output_file}")
        
        ds = xr.open_dataset(output_file).load()
        ds = ds.assign(node_id=('node', ds.node.values + 1))# Add 1-based node_id coordinate
        # But select using converted indices
        mask = ds.node_id.isin(mesh_ids)
        indices = np.where(mask)[0]

        # Step 3: Use those indices to subset
        ds_subset = ds.isel(node=indices)
        
        # Check what the node dimension looks like in the subset
        print(f"ds_subset.node values (first 5): {ds_subset.node.values[:5]}")
        ds_concise = ds_subset.drop_vars(['element', 'nvdll', 'ibtypee', 'nbvv', 'adcirc_mesh', 'nbdv'])  # Drop the node coordinate for clarity
        ds_concise.to_netcdf(savefile, encoding={v: {"zlib": False} for v in ds_concise.data_vars})
    # if file_ext == ".nc":
    #     # NetCDF input
        
    #     print(f"📊 Processing NetCDF file: {output_file}")

    #     # Load the dataset
    #     ds = xr.open_dataset(output_file).load()
    #     mesh_indices = mesh_ids - 1  # Convert to 0-based indices
    #     # Subset by mesh_node_id
    #     ds_subset = ds.isel(node=mesh_indices)
    #     # # Reorder to match fort.14 node order
    #     sort_order = df_map.sort_values("subset_node_id").index
    #     # print(sort_order)
    #     ds_ordered = ds_subset.isel(node=sort_order)

    #     # Save result
    #     ds_ordered.to_netcdf(
    #         savefile, encoding={v: {"zlib": False} for v in ds_ordered.data_vars}
    #     )

    elif file_ext in [".63", ".64"]:
        # ASCII input
        print(f"📝 Processing ASCII file: {output_file}")

        # Convert ASCII to xarray
        ds = output_to_xarray(output_file)

        # Subset the data
        ds_subset = ds.isel(node=mesh_ids)

        # Save the subset as NetCDF
        save_netcdf(ds_subset, savefile)

    else:
        raise ValueError(
            f"Unrecognized file extension: {file_ext}. Supported formats: .nc, .63, .64"
        )

    print(f"✅ Subset data saved to {savefile}")
    return savefile


if __name__ == "__main__":
    mesh = Mesh.from_fort14("/Users/crj85568/Projects/KingsBay/fort.14")
    mask = mesh.subset("rectangle" ** [-74.052297, 32.189362, -67.020014, 24.487077])
    mesh = mesh.subset(mask)
    mesh.to_fort14("ocean_fort.14")
