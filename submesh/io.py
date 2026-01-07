"""
ADCIRC file I/O Module
This module provides functions to read and write ADCIRC mesh and output files.
It includes various formats (ASCII, NetCDF, XMDF) and handles unstructured meshes.
    - Mesh files fort.14
    - Water elevation files fort.63, maxele.63
    - Velocity files fort.64, maxvel.64
    - Format conversion between ADCIRC ascii and NetCDF, and XMDF

    Typical Usage:
        1. Read a mesh file:
            nodes, elements, boundaries = read_fort14('fort.14')
        2. Read output files:
            ds = output_to_xarray('fort.63')
        3. Write to NetCDF:
            save_netcdf(ds, 'output.nc')
        4. Convert fort.14 to NetCDF:
            fort14_to_netcdf(fort14_dict, 'output.nc')
"""
import os
import re
import math
import uuid
import logging
import tempfile
from io import StringIO
from typing import Dict, List, Tuple, Optional, Union, Any

import pandas as pd
import dask.array as da
import xarray as xr
import numpy as np
import h5py

# Format templates
NODE_FMT = "{:10d} {:18.10f} {:18.10f} {:18.10e}\n"
ELEMENT_FMT = "{:10d} 3 {:10d} {:10d} {:10d}\n"
BOUNDARY_COUNT_FMT = "{:10d}\n"
_HEADER_LINE1 = "  {rundes:<32}  {runid:<24} {agrid:<24}\n"
_HEADER_LINE2 = "{NDSETSE:11d}{NP:11d}{formatted_dtdp.rjust(16)}{NSPOOLGE:9d}{IRTYPE:6d} FileFmtVersion:    1050624\n"

# Data-line formats
_TIME_FMT   = "{time:22.10E}{iter:15d}\n"
_DATA_FMT63 = "{idx:10d}{elev:22.10f}\n"
_DATA_FMT64 = "{idx:10d}{uvel:22.10f}{vvel:22.10f}\n"


def format_header(ds, irtype):
    rundes = "SAB"; runid="Tides"; agrid="grid"
    NP       = ds.dims["node"]
    times0   = ds.time.values[0]
    NDSETSE  = ds.dims["time"]
    formatted = custom_sci_format(times0, 7)
    return (
        _HEADER_LINE1.format(rundes=rundes, runid=runid, agrid=agrid),
        _HEADER_LINE2.format(
            NDSETSE=NDSETSE, NP=NP,
            formatted_dtdp=formatted,
            NSPOOLGE=int(times0*10), IRTYPE=irtype
        )
    )

def format_time_line(t, iteration):
    return _TIME_FMT.format(time=t, iter=iteration)

def format_data_block63(node_ids, elevations):
    # returns an array of strings or a single big string
    return "\n".join(
        _DATA_FMT63.format(idx=i, elev=e)
        for i, e in zip(node_ids, elevations)
    ) + "\n"

def format_data_block64(node_ids, uvel, vvel):
    return "\n".join(
        _DATA_FMT64.format(idx=i, uvel=u, vvel=v)
        for i, u, v in zip(node_ids, uvel, vvel)
    ) + "\n"

def read_fort14(filepath):
    """
    Reads a fort.14 ADCIRC mesh file.
    Args: 
        filepath (str): Path to the fort.14 file
    Returns:
    dict: Dictionary containing mesh data with keys:
        description: Description string from the file
        nodes: Dataframe with columns ['node_id', 'x', 'y', 'z']
        elements: Dataframe with columns ['element_id', 'node_1', 'node_2', 'node_3']
        boundaries: dict with keys: ['open', 'land'], each key is a list of the nodes for that boundary

    Raises:
        FileNotFoundError: If the file does not exist
        ValueError: If the file format is incorrect or if there are parsing errors
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File {filepath} not found.")
    
    output = {'description': None, 'nodes': None, 'elements': None, 'boundaries': {}}
    try:
        with open(filepath) as f:
            description = f.readline().strip()
            size_line = f.readline().strip()
            try:
                num_ele, num_nodes = [int(value) for value in size_line.split()]
            except ValueError as e:
                raise ValueError(f"Invalid mesh size line: {size_line}. Expected two integers")
            nodes = _read_nodes(f, num_nodes)
            elements = _read_elements(f, num_ele)
            open_boundaries = _read_open_bounds(f)
            land_boundaries = _read_land_bounds(f)

        return {
            'description': description,
            'nodes': nodes,
            'elements': elements,
            'boundaries': {'open': open_boundaries,
                        'land': land_boundaries}}
    except Exception as e:
        raise ValueError(f"Error reading fort.14 file: {e}")

def _read_nodes(f, num_nodes):
    """
    Reads the node data from the fort.14
    
    Args:
        f: Open file handle position at the start of the node data
        num_nodes(int): Number of nodes to read
    Returns:
        DataFrame: DataFrame containing node data with columns ['id', 'x', 'y', 'z']
    """
    # Read all node lines at once
    node_lines = ''.join(f.readline() for _ in range(num_nodes))
    nodes = pd.read_csv(StringIO(node_lines),
                        sep=r'\s+',
                        index_col=0,
                        names=['id', 'x', 'y', 'z'],
                        dtype={'id': int, 'x': float, 'y': float, 'z': float}, # Ensure correct types)
    )
    return nodes

def _read_elements(f, num_ele):
    """
    Reads element data from fort.14
    Args: 
        f: Open file handle position at the start of the element data
        num_ele(int): Number of elements to read
    Returns:
        DataFrame: DataFrame containing element data with columns ['id', 'n', 'node_1', 'node_2', 'node_3']
    """
    element_lines = ''.join(f.readline() for _ in range(num_ele))
    elements = pd.read_csv(
        StringIO(element_lines),
        sep=r"\s+",
        index_col=0,
        names=['id', 'n', 'node_1', 'node_2', 'node_3'],
        dtype={'id': int, 'n': int, 'node_1': int, 'node_2': int, 'node_3': int} # Ensure correct types
    )
    return elements
    
def _read_open_bounds(f):
    """
    Read open boundary data from an open fort.14 file handle.

    Args:
        f: Open file handle position at the start of the open boundary data

    Returns:
        List[BoundarySegment]: List of open boundary segments
    Raises:
        ValueError: If the file format is incorrect or if there are parsing errors
    """
    from submesh.boundaries import BoundaryNode, BoundarySegment

    open_bounds = []
    try:
        num_open_boundaries = int(f.readline().split(maxsplit=1)[0])
        f.readline()  # Skip total boundary nodes line

        for _ in range(num_open_boundaries):
            num_boundary_nodes = int(f.readline().split(maxsplit=1)[0])

            # Read node IDs and create BoundaryNode objects
            nodes = tuple(
                BoundaryNode(node_id=int(f.readline().split(maxsplit=1)[0].strip()))
                for _ in range(num_boundary_nodes)
            )

            segment = BoundarySegment(
                nodes=nodes,
                category='open',
                ibtype=None
            )
            open_bounds.append(segment)

    except (IndexError, ValueError) as e:
        raise ValueError(f"Error parsing open boundaries: {e}")

    return open_bounds

def _read_land_bounds(f):
    """
    Read land boundary data from an open fort.14 file handle.

    Args:
        f: Open file handle position at the start of the land boundary data
    Returns:
        List[BoundarySegment]: List of land boundary segments
    Raises:
        ValueError: If the land boundary data is malformed or contains unsupported types
    """
    from submesh.boundaries import (
        BoundaryNode, BoundarySegment,
        IBTYPE_SIMPLE, IBTYPE_BARRIER, IBTYPE_CONNECTED
    )

    land_bounds = []
    try:
        num_land_boundaries = int(f.readline().split()[0])
        total_land_boundaries = int(f.readline().split()[0])

        for _ in range(num_land_boundaries):
            num_nodes, ibtype = _parse_header_line(f.readline())
            nodes = []

            for _ in range(num_nodes):
                tokens = f.readline().split()
                node_id = int(tokens[0])

                if ibtype in IBTYPE_SIMPLE:
                    # Simple boundary: just node ID
                    nodes.append(BoundaryNode(node_id=node_id))

                elif ibtype in IBTYPE_BARRIER:
                    # Barrier boundary: node + 2 barrier heights
                    barrier = (float(tokens[1]), float(tokens[2]))
                    nodes.append(BoundaryNode(
                        node_id=node_id,
                        barrier_heights=barrier
                    ))

                elif ibtype in IBTYPE_CONNECTED:
                    # Connected boundary: node + connected_to + 3 barrier values
                    connected_to = int(tokens[1])
                    barrier = (float(tokens[2]), float(tokens[3]), float(tokens[4]))
                    nodes.append(BoundaryNode(
                        node_id=node_id,
                        connected_to=connected_to,
                        barrier_heights=barrier
                    ))

                else:
                    raise ValueError(f"Unsupported IBTYPE: {ibtype}")

            segment = BoundarySegment(
                nodes=tuple(nodes),
                category='land',
                ibtype=ibtype
            )
            land_bounds.append(segment)

    except (IndexError, ValueError) as e:
        raise ValueError(f'Error while reading land boundaries: {e}')

    return land_bounds

def _parse_header_line(line):
    """
    Parse a healder line containing integers

    Args:
        line (str): Line to parse

    Returns:
        tuple: Pair of integers extracted from the line

    Raises:
        ValueError: If the line doesn't contain two integers
    """
    # Split the line into tokens and convert to integers
    tokens = line.strip().split()
    int_tokens = [int(tok) for tok in tokens if tok.isdigit()]

    if len(int_tokens) < 2:
        raise ValueError(f"Expected two integers in line, got: {line}")
    
    return int_tokens[0], int_tokens[1]


def _parse_header(f):
    line = f.readline()
    return re.findall(r"(\S+)\s*!.*?(\S+)\s*!", line)[0]

def _parse_metadata(f):
    line = f.readline()
    # Regex pattern to extract the first five values before "FileFmtVersion"
    pattern = r"^\s*(\d+)\s+(\d+)\s+([\d.Ee+-]+)\s+(\d+)\s+(\d+)\s+FileFmtVersion"
    match = re.search(pattern, line)
    if match:
        # Convert extracted values to their respective types
        return {
            'num_datasets' : int(match.group(1)),
            'num_nodes' : int(match.group(2)),
            'timestep' : float(match.group(3)),
            'record_type' : int(match.group(5)),
            'grid_id': None}
        
def _get_format_from_extension(filepath: str) -> str:
    """ Determine file format from file extension"""
    _, ext = os.path.splitext(filepath)
    ext = ext.lower()

    if ext == '.nc':
        return 'nc'
    elif ext in ['.63', '.64', '.14']:
        return 'ascii'
    else:
        return 'unknown'

def _parse_elev_datasets(f, num_nodes, num_datasets):
    columns = ['id', 'eta']
    timesteps = []
    nodes = []
    for _ in range(num_datasets):
        line = f.readline()
        current_ts = line.strip().split()[0]
        timesteps.append(current_ts)
        # Read num_nodes lines for dataset
        with StringIO("\n".join(f.readline() for _ in range(num_nodes))) as nodes_stream:
            nodes.append(pd.read_csv(nodes_stream, sep=r'\s+',
                                                index_col=0, names=columns))
    return timesteps, nodes

def _parse_velo_datasets(f, num_nodes, num_datasets):
    columns = ['node', 'u-vel', 'v-vel']
    timesteps = []
    nodes = []
    for time in range(num_datasets):
        current_ts = f.readline().strip().split()[0]
        timesteps.append(current_ts)
        with StringIO("\n".join(f.readline() for _ in range(num_nodes))) as nodes_stream:
            # nodes.append(pd.read_csv(nodes_stream, sep=r'\s+',
            #                          index_col=0, names=columns))
            df = pd.read_csv(nodes_stream, sep=r'\s+',index_col=0, names=columns)
            
            df['mag'] = np.hypot(df['u-vel'], df['v-vel'])
            ds = df.to_xarray()
            ds = ds.expand_dims(time=[float(current_ts)])
            if not os.path.exists('temp_nc'):
                os.mkdir('temp_nc')
            ds.to_netcdf(f"temp_nc/fort.64.{time}.nc")
        # for df in nodes:
        #     df['mag'] = np.hypot(df['u-vel'], df['v-vel'])
    return timesteps, nodes


def read_output(f):
    description, run_id = _parse_header(f)
    metadata = _parse_metadata(f)
    if metadata['record_type'] == 1:
        timesteps, nodes = _parse_elev_datasets(f, metadata['num_nodes'], metadata['num_datasets'])
    elif metadata['record_type'] == 2:
        timesteps, nodes = _parse_velo_datasets(f, metadata['num_nodes'], metadata['num_datasets'])
    else:
        raise ValueError("Unsupported record type")
        
    return {'description': description,
            'run_id': run_id,
            'grid_id': metadata.get('grid_id',None),
            'nodes': nodes,
            'timesteps': timesteps,
            'record_type': metadata['record_type']}

def output_to_xarray(file):
    """
    Reads an adcirc output (fort.63, fort.64, etc) file object and converst to an xarray
    dataset, assumes each dataset in 'nodes' is a dataframe with node data
    and that 'timesteps' provides the time coordinate
    """
    with open(file) as f:
        data = read_output(f)
    times = [float(t) for t in data['timesteps']]
    record_cols = {
        1: ['eta'],
        2: ['u-vel', 'v-vel', 'mag']
    }
    try:
        cols = record_cols[data['record_type']]
    except KeyError:
        raise ValueError(f"Unsupported record_type: {data['record_type']}")
    if cols == ['u-vel', 'v-vel', 'mag']:
        ds = xr.open_mfdataset(
            f"temp_nc/fort.64.*.nc",
            combine='by_coords',
            parallel=True,
            chunks={'time': 1, 'node': -1}
        )
    else:
        ds_list = []
        for time, df in zip(times, data['nodes']):
            vars_dict = {
                col: (('node',), da.from_array(df[col].values, chunks=(len(df),)))
                for col in cols
            }
            ds_temp = xr.Dataset(
                data_vars=vars_dict,
                coords = {'node': df.index.values}
            )
            ds_temp = ds_temp.expand_dims(time=[time])
            ds_list.append(ds_temp)
        ds = xr.concat(ds_list, dim='time')
        ds.attrs.update({
                'description': data['description'],
                'run_id': data['run_id'],
                'record_type': data['record_type']})
    return ds

def combine_adcirc_data(fort14_file, fort63_file, fort64_file=None, maxele=None, maxvel=None):
    # Read fort.14 and fort.63 files.
    fort14_data = read_fort14(fort14_file)
    ds_out = output_to_xarray(fort63_file)
    
    
    # Extract node coordinates from fort.14.
    nodes_df = fort14_data['nodes']
    # Assumes the index is the node ID.
    x_coords = nodes_df['x'].values
    y_coords = nodes_df['y'].values
    node_ids = nodes_df.index.values

    ds_out = ds_out.assign_coords(x=(("node",),x_coords))
    ds_out = ds_out.assign_coords(y=(("node",), y_coords))
    # Build the connectivity array from fort.14 elements.
    # Each row is a triangle defined by its 3 node IDs.
    elements_df = fort14_data['elements']
    connectivity = elements_df[['node_1', 'node_2', 'node_3']].values

    ds_out['connectivity'] = (("element", "vertex"), connectivity)
    if fort64_file:
        ds_velo = output_to_xarray(fort64_file)
        ds_out['x_vel'] = ds_velo['xvel']
        ds_out['y_vel']
    elif maxele:
        max_elev = output_to_xarray(maxele)
        ds_out['max_ele'] = max_elev['max_ele']
    elif maxvel:
        # max_velo = 
        pass
    
    ds_out = ds_out.assign_coords(
        element=np.arange(connectivity.shape[0]),
        vertex=np.arange(connectivity.shape[1]))
    
    return ds_out

def save_netcdf(ds, savepath='adcirc_output.nc'):
    # Write out to a netCDF file.
    # You can add the 'chunks' parameter to ds.to_netcdf if you wish to leverage Dask for writing.
    ds.to_netcdf(savepath, encoding={v: {"zlib": False} for v in ds.data_vars})
    # ds.to_netcdf(savepath, mode='w', engine="netcdf4")



def fort14_to_netcdf(fort14_dict, savepath='fort14_mesh.nc'):
    """
    Converts a parsed fort.14 dictionary (from read_fort14) into an xarray.Dataset
    and saves it to NetCDF.
    """
    nodes_df = fort14_dict['nodes']
    elements_df = fort14_dict['elements']

    ds = xr.Dataset(
        data_vars={
            'z': (('node',), nodes_df['z'].values),
            'connectivity': (('element', 'vertex'), elements_df[['node_1', 'node_2', 'node_3']].values)
        },
        coords={
            'node': nodes_df.index.values,
            'x': ('node', nodes_df['x'].values),
            'y': ('node', nodes_df['y'].values),
            'element': elements_df.index.values,
            'vertex': [0, 1, 2]
        },
        attrs={
            'description': fort14_dict.get('description', 'ADCIRC Mesh')
        }
    )

    ds.to_netcdf(savepath)
    print(f"✅ Saved mesh to {savepath}")

def write_fort14(filepath: str, mesh: 'Mesh'):
    """
    Write a fort.14-style mesh file from a Mesh dataclass.

    Args:
        filepath (str): Path to output file.
        mesh (Mesh): Mesh dataclass instance.
    """
    with open(filepath, 'w') as f:
        f.write(mesh.description + "\n")
        f.write(f"{len(mesh.elements)} {len(mesh.nodes)}\n")

        for node_id, row in mesh.nodes.iterrows():
            f.write(f"{node_id} {row['x']:18.10f} {row['y']:18.10f} {row['z']:18.10e}\n")

        for i, row in mesh.elements.iterrows():
            f.write(f"{i} 3 {row['node_1']} {row['node_2']} {row['node_3']}\n")

        # Write open boundaries
        open_segments = mesh.boundaries.get('open', [])
        f.write(f"{len(open_segments)}\n")

        total_nope = sum(len(seg.nodes) for seg in open_segments)
        f.write(f"{total_nope}\n")

        for segment in open_segments:
            f.write(f"{len(segment.nodes)}\n")
            for node in segment.nodes:
                f.write(f"{node.node_id}\n")

        # Write land boundaries
        land_segments = mesh.boundaries.get('land', [])
        f.write(f"{len(land_segments)}\n")

        # Total land boundary nodes
        total_nvell = sum(len(seg.nodes) for seg in land_segments)
        f.write(f"{total_nvell}\n")

        for segment in land_segments:
            ibtype = segment.ibtype
            f.write(f"{len(segment.nodes)} {ibtype}\n")

            # Write each node based on its attributes
            for node in segment.nodes:
                if node.is_simple():
                    f.write(f"{node.node_id}\n")

                elif node.is_barrier():
                    h1, h2 = node.barrier_heights
                    f.write(f"{node.node_id} {h1:18.10f} {h2:18.10f}\n")

                elif node.is_connected():
                    v1, v2, v3 = node.barrier_heights
                    f.write(f"{node.node_id} {node.connected_to} "
                           f"{v1:18.10f} {v2:18.10f} {v3:18.10f}\n")

    print(f"✅ fort.14 written to {filepath}")

def read_fort63_ascii(file):
    output = {'description': None,
              'run_id': None,
              'grid_id': None,
              'datasets': None,
              'nodes': [],
              'timesteps': [],
              'record_type': None}
    columns = ['id', 'eta']
    
    
    with open(file) as f:
        line = f.readline()
        output['description'], output['run_id'] = re.findall(r"(\S+)\s*!.*?(\S+)\s*!", line)[0]
        # Extract everything before "FileFmtVersion"
        line = f.readline()
       # Regex pattern to extract the first five values before "FileFmtVersion"
        pattern = r"^\s*(\d+)\s+(\d+)\s+([\d.Ee+-]+)\s+(\d+)\s+(\d+)\s+FileFmtVersion"

        match = re.search(pattern, line)
        if match:
            # Convert extracted values to their respective types
            num_datasets = int(match.group(1))
            num_nodes = int(match.group(2))
            timestep = float(match.group(3))
            total_time = int(match.group(4))
            output['record_type'] = int(match.group(5))
        for time in range(num_datasets):
            current_ts = f.readline().strip().split()[0]
            output['timesteps'].append(current_ts)
            with StringIO("\n".join(f.readline() for _ in range(num_nodes))) as nodes_stream:
                output['nodes'].append(pd.read_csv(nodes_stream, sep=r'\s+',
                                                    index_col=0, names=['id', f'{current_ts}']))
    return output


def read_fort63_nc(filepath):
    ds = xr.open_dataset(filepath)
    return ds

def load_subgrid_mapping(mapping_file):
    """Loads the subgrid node id mapping from original -> subgrid

    Args:
        mapping_file (pd.Series): Mapping with subgrid_node ID as the column, global_id as index
    """
    df = pd.read_csv(mapping_file)
    return df

def write_fort63_file(output_file, subset_file="subset_fort63.63"):
    """Writes subset elevation data to fort.63 format."""
    if output_file.endswith('.nc'):
        ds = xr.open_dataset(output_file, decode_times=False)
        rundes = "SAB"
        runid = "Tides"
        agrid = 'grid'
        NP = ds.zeta.shape[1]
        DTDP = ds.time[0]
        NSPOOLGE = int(ds.time[0].values * 10)
        IRTYPE = 1
        NDSETSE = len(ds.time)
        formatted_dtdp = custom_sci_format(DTDP, 7)
        buf = np.empty((NP, 2), dtype=np.float64)
        fmt = "%10d%22.10f"
        with open(subset_file, "w") as f:
            # Write metadata (modify if necessary)
            f.write(f"  {rundes:<32}  {runid:<24} {agrid:<24}\n")
        
            f.write(f"{NDSETSE:11d}{NP:11d}{formatted_dtdp.rjust(16)}{NSPOOLGE:9d}{IRTYPE:6d} FileFmtVersion:    1050624\n")
            times = ds.time.values
            num_timesteps = (times*10).astype(int)
            eta = ds.zeta.values
            time_headers = [f"{t:22.10E}{nt:15d}\n" for t, nt in zip(times, num_timesteps)]
            # 3) loop over each time‐step chunk
            for i, (t, eta) in enumerate(zip(times, eta)):
                #    TIME (as float) , IT (tidal iteration—often 0)
                f.write(time_headers[i])
                #    then NP lines: node_index, elevation
                buf[:, 0] = np.arange(1, NP+1)
                buf[:, 1] = eta
            
                np.savetxt(f, buf, fmt=fmt)

    else:
        raise ValueError("Output file must be a NetCDF file (.nc)")


    print(f"✅ Subset data saved to {subset_file}")

def custom_sci_format(value, decimals=7):
    import math

    if value == 0:
        return f"0.{''.join(['0']*decimals)}E+000"

    exponent = int(math.log10(abs(value))) + 1
    mantissa = value / (10 ** exponent)
    format_str = f"{{:.{decimals}f}}E+{exponent:03d}"
    return format_str.format(mantissa)


def write_fort64_file(output_file, subset_file):
    """Writes subset elevation data to fort.63 format."""
    if output_file.endswith('.nc'):
        ds = xr.open_dataset(output_file, decode_times=False)
        ds.load()
        rundes = "SAB"
        runid = "Tides"
        agrid = 'grid'
        NP = ds['u-vel'].shape[1]
        DTDP = ds.time[0]
        NSPOOLGE = int(ds.time[0].values * 10)
        IRTYPE = 2
        NDSETSE = len(ds.time)
        # right after opening ds:
        u_data = np.nan_to_num(ds['u-vel'].values, nan=-99999.0)
        v_data = np.nan_to_num(ds['v-vel'].values, nan=-99999.0)
        times  = ds.time.values

        formatted_dtdp = custom_sci_format(DTDP, 7)
        print(formatted_dtdp)
        node_ids = np.arange(1, NP+1)
        buf = np.empty((NP, 3), dtype=np.float64)
        fmt = "%10d%22.10f%22.10f"
        with open(subset_file, "w", buffering=4*1024*1024) as f:
            # Write metadata (modify if necessary)
            f.write(f"  {rundes:<32}  {runid:<24} {agrid:<24}\n")
        
            f.write(f"{NDSETSE:11d}{NP:11d}{formatted_dtdp.rjust(16)}{NSPOOLGE:9d}{IRTYPE:6d} FileFmtVersion:    1050624\n")
            
            # 3) loop over each time‐step chunk
            for t, uvel, vvel in zip(times, u_data, v_data):
                #    TIME (as float) , IT (tidal iteration—often 0)ls
                f.write(f"{t:22.10E}{int(t*10):15d}\n")
               # before the t-loop, once:
                buf[:, 0] = node_ids
                buf[:, 1] = uvel
                buf[:, 2] = vvel
                np.savetxt(f, buf, fmt=fmt, delimiter="")
                # #    then NP lines: node_index, elevation
                # for k, vel in enumerate(uvel, start=1):
                #     v = vvel[k-1]
                   
                #     f.write(f"{k:10d}{vel:22.10f}{v:22.10f}\n")
    else:
        raise ValueError("Output file must be a NetCDF file (.nc)")

def create_xmdf_hdf5(file_path, dataarrays):
    """
    Writes one or more xarray Dataset/DataArray objects to an XMDF-style HDF5 file.
    
    Parameters:
    - file_path: str, path to the output .h5 file.
    - dataarrays: dict mapping variable name to either an xarray DataArray or
      an xarray Dataset (with 1 or more data variables).
    """
    with h5py.File(file_path, 'w') as f:
        # Root-level metadata using bytes_ strings
        f.create_dataset('File Type',    data=np.bytes_('XMDF'))
        f.create_dataset('File Version', data=np.bytes_('1.0'))
        f.create_dataset('Guid',         data=np.bytes_(str(uuid.uuid4())))
        
        # Top-level group for variables
        dsg = f.create_group('Datasets')
        
        for name, da in dataarrays.items():
            # If a Dataset is passed:
            # … inside your loop over dataarrays …
            if isinstance(da, xr.Dataset):
                # 1) pick only the vars that look like time‐series on nodes
                time_vars = [
                    v for v in da.data_vars
                    if da[v].ndim == 2
                    and da[v].dims[0] == 'time'
                    and da[v].dims[1] == 'node'
                ]
                if len(time_vars) == 0:
                    raise ValueError(f"No time×node variables found in Dataset '{name}'")
                elif len(time_vars) == 1:
                    # just one: treat it as a scalar
                    da = da[time_vars[0]]
                elif len(time_vars) == 2:
                    # exactly two: stack into a vector
                    arr1, arr2 = da[time_vars[0]], da[time_vars[1]]
                    # 2) align them on exactly the same coords
                    a1, a2 = xr.align(arr1, arr2, join='inner')
                    # 3) build a new DataArray with a 'component' axis
                    vals = np.stack([a1.values, a2.values], axis=-1)
                    da = xr.DataArray(
                        data=vals,
                        dims=('time','node','component'),
                        coords={
                            'time': a1.coords['time'],
                            'node': a1.coords['node'],
                            'component': time_vars
                        }
                    )
                else:
                    raise ValueError(
                        f"Found {len(time_vars)} time×node variables in '{name}', "
                        "but expected at most 2."
                    )

            
            # Now da is a DataArray
            values = da.values
            times = da.coords['time'].values
            
            # Determine node count for naming
            node_count = int(np.prod(values.shape[1:])) if values.ndim > 1 else 1
            grp_name = f"{name} ({node_count})"
            ds_grp = dsg.create_group(grp_name)
            
            # PROPERTIES subgroup
            props = ds_grp.create_group('PROPERTIES')
            props.create_dataset('Object Type', data=np.bytes_('Transient Dataset'))
            nullval = da.attrs.get('nullvalue', None)
            if nullval is not None:
                props.create_dataset('nullvalue', data=np.array(nullval).astype('f'))
            
            # Core datasets
            ds_grp.create_dataset('Times',  data=times)
            ds_grp.create_dataset('Values', data=values)
            ds_grp.create_dataset('Maxs',   data=np.nanmax(values, axis=tuple(range(1, values.ndim))))
            ds_grp.create_dataset('Mins',   data=np.nanmin(values, axis=tuple(range(1, values.ndim))))
            
            # For vector fields, add magnitude
            if values.ndim > 2:
                mag_grp = ds_grp.create_group('mag')
                mag_props = mag_grp.create_group('PROPERTIES')
                mag_props.create_dataset('Object Type', data=np.bytes_('Transient Dataset'))
                
                mag_vals = np.linalg.norm(values, axis=-1)
                mag_grp.create_dataset('Times',  data=times)
                mag_grp.create_dataset('Values', data=mag_vals)
                mag_grp.create_dataset('Maxs',   data=np.nanmax(mag_vals, axis=1))
                mag_grp.create_dataset('Mins',   data=np.nanmin(mag_vals, axis=1))

# Usage example:
# write_xmdf('datasets.h5', {
#     'Velocity': xr.open_dataset('TaylorBOunds/fort.64.nc'),
#     'Water Elev': xr.open_dataset('TaylorBOunds/fort.63.nc')
# })
