"""Mesh merging module for recombining refined submeshes with global mesh.
This module provides robust mesh merging that handles:
1. Unmodified nodes (fast CSV mapping lookup)
2. Refined regions (spatial KD-tree matching)
3. Boundary reconstruction
4. Element selection and renumbering

Key features:
- Two-phase node matching (CSV map + spatial index)
- Provenance tracking (original vs refined nodes)
- Handles resolution changes in submesh
- Preserves boundaries correctly
"""
import numpy as np
import pandas as pd
from submesh.submesh import Mesh
from submesh.boundaries import IBTYPE_CONNECTED
from submesh.boundaries import BoundaryNode, BoundarySegment
from submesh.io import write_fort14
import time
from pathlib import Path
from collections import Counter
from numba import njit


@njit
def mark_boundaries(nblg, nblncg, nbl, nblnc, npline):
    global_edges = set()
    for i in range(nblg):
        a = nblncg[0, i]
        b = nblncg[1, i]
        global_edges.add((min(a,b), max(a,b)))

    for i in range(nbl):
        a = nblnc[0, i]
        b = nblnc[1, i]
        if (min(a,b), max(a,b)) in global_edges:
            npline[i] = 0
    return npline, global_edges

from collections import defaultdict

def boundary_nodes_fortran_df(df):
    nprops = defaultdict(int)
    tris = df[["node_1", "node_2", "node_3"]]

    for a, b, c in tris.itertuples(index=False):
        nprops[a] += b - c
        nprops[b] += c - a
        nprops[c] += a - b

    nbn_max = sum(1 for v in nprops.values() if v != 0)
    
    return nbn_max, nprops


def merge_domain(submesh, full_mesh, path='.', file_header='grid'):
    """Merge a refined submesh back into the full mesh.

    Args:
        submesh: Mesh object for the refined submesh
        full_mesh: Mesh object for the original full mesh
        path: Directory path where node_id_map.csv is located (default: current directory)
        file_header: Header for output file (default: 'grid')

    Returns:
        None (writes merged.14 to disk)
    """
    print("\n" + "="*80)
    print("MERGE DOMAIN - Progress Tracking")
    print("="*80)

    num_nodes = len(full_mesh.nodes)
    num_elements = len(full_mesh.elements)
    print(f"Mesh size: {num_nodes:,} nodes, {num_elements:,} elements\n")

    # Step 1: Find boundary of entire global grid
    print("[Step 1/10] Finding global grid boundary (netable=all 1s)...")
    nbn_max, nprop1 = boundary_nodes_fortran_df(full_mesh.elements)
    nbn_max = 20
    t0 = time.time()
    netable = np.ones(num_elements, dtype=np.int32)
    nblg, nblncg = mkeline(num_nodes, num_elements, full_mesh.elements, netable, nbn_max)
    print(f"            Found {nblg:,} global boundary edges ({time.time()-t0:.1f}s)\n")

    # Step 2: Read netable from element mask file
    # Reads the elements from the global grid which are contained 
    # in the subgrid (True) or not (False)
    t0 = time.time()
    element_mask_path = Path(path) / 'element_mask.csv'
    print(f"[Step 2/10] Reading element mask from {element_mask_path}...")
    element_mask_df = pd.read_csv(element_mask_path, index_col=0)
    netable = element_mask_df['in_subgrid'].astype(np.int32).values
    nptable = [0] * (num_nodes + 1)
    nprops = np.zeros(num_elements + 1)
    for m, elem in full_mesh.elements.iterrows():
        if netable[m-1] == 1:
            a = elem['node_1']
            b = elem['node_2']
            c = elem['node_3']
            nptable[a -1] = 1
            nptable[b- 1] = 1
            nptable[c -1] = 1
            nprops[a] += b - c
            nprops[b] += c - a
            nprops[c] += a - b

    nbn_max = 20
    for i in range(num_nodes):
        if nprops[i] !=0:
            nbn_max += 1

    # Step 3: Find boundary around subgrid
    t0 = time.time()
    print("[Step 3/10] Finding subgrid boundary edges...")
    
    nbl, nblnc = mkeline(num_nodes, num_elements, full_mesh.elements, netable, nbn_max)
    print(f"            Found {nbl:,} subgrid boundary edges ({time.time()-t0:.1f}s)\n")
    # Step 4: Compare boundary edges and mark differences
    t0 = time.time()
    print("[Step 4/10] Comparing global vs subgrid boundaries...")
    # Create set of global boundary edges for O(1) lookup
    npline = np.ones(nbl)
    # npline, edge_lookup = mark_boundaries(nblg, nblncg, nbl, nblnc, npline)
    for m in range(nbl):
        m1 = nblnc[0, m]
        m2 = nblnc[1, m]
        for n in range(nblg):
            n1 = nblncg[0, n]
            n2 = nblncg[1,n]
            if (m1==n1) and (n2==m2):
                npline[m] = 0
                break
    for segments in full_mesh.boundaries.values():
        for segment in segments:
            ibtype = segment.ibtype
            if (ibtype !=4) and (ibtype!=24):
                continue
            nodes = segment.nodes
            nvell = len(nodes)
            for j in range(1, nvell+1):
                n1 = nodes[j-1].node_id
                n2 = nodes[j].node_id
                l1 = nodes[j-1].connected_to or n1
                l2 = nodes[j].connected_to or n2
                
                nchk = 0
                for m in range(1, nbl+1):
                    m1 = nblnc[1][m]
                    m2 = nblnc[2][m]
                    if (m1==n1 and m2 == n2) or (m1 == n2 and m2==n1):
                        nchk = m
                        break
                lchk = 0
                for m in range(1, nbl+1):
                    m1 = nblnc[1][m]
                    m2 = nblnc[2][m]
                    if (m1==l1 and m2==l2) or (m1==l2 and m2==l1):
                        lchk = m
                        break
                if (nchk * lchk == 0) and (nchk + lchk != 0):
                    npline[nchk + lchk] = 1

    nbd = int(np.sum(npline))
    nbdnc = np.zeros((2, nbd))
    nbd = 0
    for n in range(nbl):
        if(npline[n] == 0):
            continue
        
        for i in range(2):
            nbdnc[i, nbd] = nblnc[i,n]
        nbd += 1
   

    nprops = np.zeros(len(submesh.elements) + 1)
    for m, elem in submesh.elements.iterrows():
        a = elem['node_1']
        b = elem['node_2']
        c = elem['node_3']
        nprops[a] += b - c
        nprops[b] += c - a
        nprops[c] += a - b

    nbns_max = 20
    for i in range(len(submesh.nodes)):
        if nprops[i] !=0:
            nbns_max += 1
    mmaps = np.ones(len(submesh.elements))
    nbls, nblncs = mkeline(len(submesh.nodes), len(submesh.elements), submesh.elements, mmaps, nbns_max)
    
    # n_connected_segs = sum(1 for seg in full_mesh.boundaries['land']
    #                       if seg.ibtype in IBTYPE_CONNECTED)
    # for segment in full_mesh.boundaries['land']:
    #     if segment.ibtype not in IBTYPE_CONNECTED:
    #         continue
    #     # Check consecutive node pairs along this boundary
    #     for j in range(len(segment.nodes) -1):
    #         node_current = segment.nodes[j]
    #         node_next = segment.nodes[j + 1]
    #         n1 = node_current.node_id
    #         n2 = node_next.node_id
    #         l1 = node_current.connected_to
    #         l2 = node_next.connected_to

    #         # check if edge(n1, n2) is in the subgrid boundary edges (O(1) lookup)
    #         nchk = edge_lookup.get((n1, n2), -1)
    #         lchk = edge_lookup.get((l1, l2), -1)

    #         if (nchk * lchk == -1) and (nchk + lchk != -1):

    #             idx = nchk if nchk != -1 else lchk
    #             npline[idx] = 1

    # print(f"            Processed {n_connected_segs} connected boundary segments ({time.time()-t0:.1f}s)\n")
    # Step 6: Create final boundary edge array
    t0 = time.time()
    print("[Step 6/10] Creating final boundary edge array...")
    # nbd = int(np.sum(npline))
    # nbdnc = np.zeros((2, nbd), dtype=np.int32)
    # nbd_counter = 0
    # for n in range(nbl):
    #     if npline[n] == 0:
    #         continue
    #     nbdnc[0, nbd_counter] = nblnc[0, n]
    #     nbdnc[1, nbd_counter] = nblnc[1,n]
    #     nbd_counter += 1
    # print(f"            Created {nbd:,} boundary edges ({time.time()-t0:.1f}s)\n")
    # print(f"            Created {nbd:,} boundary edges ({time.time()-t0:.1f}s)\n")

    # # ADD THIS:
    # print("            Boundary edge node IDs:")
    # for i in range(min(nbd, 10)):  # Show first 10
    #     n1, n2 = nbdnc[0, i], nbdnc[1, i]
    #     x1, y1 = full_mesh.nodes.loc[n1, 'x'], full_mesh.nodes.loc[n1, 'y']
    #     x2, y2 = full_mesh.nodes.loc[n2, 'x'], full_mesh.nodes.loc[n2, 'y']
    #     print(f"              Edge {i}: ({n1}, {n2}) - ({x1:.2f}, {y1:.2f}) to ({x2:.2f}, {y2:.2f})")
    # # Step 7: Find submesh boundaries
    # t0 = time.time()
    # print("[Step 7/10] Finding submesh boundary edges...")
    # mmaps = np.ones(len(submesh.elements), dtype=np.int32)
    # nbls, nblncs = mkeline(len(submesh.nodes), len(submesh.elements), submesh.elements, mmaps)
    # print(f"            Found {nbls:,} submesh boundary edges ({time.time()-t0:.1f}s)\n")


    # Step 8: Match boundary nodes between meshes
    t0 = time.time()
    print("[Step 8/10] Matching boundary nodes between global and submesh...")
    print(" NOW CHECKING! the pair nodes of broders to submesh")
    # Create position-to-ID and ID-to-position mappings for both meshes
    global_node_ids = full_mesh.nodes.index.to_numpy()
    submesh_node_ids = submesh.nodes.index.to_numpy()
    global_elem_ids = full_mesh.elements.index.to_numpy()

    # Create reverse mappings: ID -> position (0-indexed)
    global_node_id_to_pos = {node_id: pos for pos, node_id in enumerate(global_node_ids)}
    submesh_node_id_to_pos = {node_id: pos for pos, node_id in enumerate(submesh_node_ids)}
    global_elem_id_to_pos = {elem_id: pos for pos, elem_id in enumerate(global_elem_ids)}

    # Arrays indexed by POSITION (0-indexed)
    nmap = np.zeros(num_nodes, dtype=int)
    nmaps = np.ones(len(submesh.nodes), dtype=int) * -1
    mmap = np.zeros(num_elements, dtype=np.int32)
    nprop1 = np.zeros(num_nodes, dtype=int)  # Store submesh POSITIONS, not IDs
    
    eps = 1e-6
    j = 0
    
    for n in range(nbd):
       
        n1 = int(nbdnc[0, n])+ 1
        n2 = int(nbdnc[1, n])+ 1
        xn1 = full_mesh.nodes.loc[n1, 'x']
        yn1 = full_mesh.nodes.loc[n1, 'y']
        xn2 = full_mesh.nodes.loc[n2, 'x']
        yn2 = full_mesh.nodes.loc[n2, 'y']
        nptable[n1 -1] = 2
        nptable[n2 -1] = 2
        i = 0
        for m in range(nbls):
            
            m1 = nblncs[0, m] + 1
            m2 = nblncs[1, m] + 1
            dx1 = xn1 - submesh.nodes.loc[m1, 'x']
            dy1 = yn1 - submesh.nodes.loc[m1, 'y']
            dx2 = xn2 - submesh.nodes.loc[m2, 'x']
            dy2 = yn2 - submesh.nodes.loc[m2, 'y']
            if (np.sqrt(dx1) <= eps and np.sqrt(dx2) <= eps):
                i += 1
                nmaps[m1-1] = n1
                nmaps[m2-1] = n2
                nprop1[n1-1] = m1
                nprop1[n2-1] = m2
                break
        if (i == 0):
            j += 1
            print(' Node on border of grids are not matching', n1, n2, 'not found in subgrid')
        if (j != 0):
            'Hit enter key to stop'
        nprops = np.zeros(len(submesh.nodes))
        nd = 0
        for n in range(len(submesh.nodes)):
            if (nmaps[n] == -1):
                nprops[nd] = n
                nd += 1
    
    t0 = time.time()
    print("[Step 9/10] Creating node and element ID mappings...")
    n1 = 0
    n2 = 0
    m = 0
    npn = 0
    nptable = np.array(nptable)
    
    for n in range(num_nodes):
        if nptable[n] == 0:
            m += 1
            nmap[n] = m
        elif nptable[n] == 1:
            if (n1 + 1 <= nd):
                m += 1
                nmaps[int(nprops[n1])] = m
                n1 += 1
        elif nptable[n] == 2:
            m += 1
            nmap[n] = m
            nmaps[nprop1[n]-1] = m
   

    for n in range(n1+1, nd+1):
        m += 1
        nmaps[int(nprops[n])] = num_nodes + n - n1

    npn = m

    # Element mapping
    print('NOW STARTING ELEMENT MAPPING')
    n = 0
    for m in range(num_elements):
        if netable[m] == 1:
            if n < len(submesh.elements):  # Boundary check from Fortran line 1246
                mmaps[n] = m
                n += 1
    nd = mmaps[n-1] if n > 0 else 0  # Use n-1 since n was incremented
    for i in range(n+1, len(submesh.elements)+1):
        mmaps[i] = nd + i - n
    nen = len(submesh.elements)
    n = 0
    for m in range(num_elements):
        if netable[m] == 0:
            nen += 1
            if m <= nd:
                mmap[m] = m
            else:
                n += 1
                mmap[m] = mmaps[len(submesh.elements) - 1] + n  # Fix: -1 for 0-based indexing
    print(f"            Created mappings for {npn:,} nodes, {nen:,} elements ({time.time()-t0:.1f}s)\n")
    
    print('NOW CREATING GLOBAL GRID')
    # allocate(xydn(3,npn), nmn(nen,3))

   
    xydn = [[0.0] * (npn + 1) for _ in range(3)]
    nmn  = [[0]   * 4          for _ in range(nen + 1)]  # index 1..3 used
    for n in range(len(submesh.nodes)):
        new_n = nmaps[n]
        xydn[0][new_n] = submesh.nodes.loc[n+1, 'x']
        xydn[1][new_n] = submesh.nodes.loc[n+1, 'y']
        xydn[2][new_n] = submesh.nodes.loc[n+1, 'z']
    
    for n in range(num_nodes):
        if nptable[n] == 0:
            new_n = nmap[n]
            xydn[0][new_n] = full_mesh.nodes.loc[n+1, 'x']
            xydn[1][new_n] = full_mesh.nodes.loc[n+1, 'y']
            xydn[2][new_n] = full_mesh.nodes.loc[n+1, 'z']

    new_nodes_df = pd.DataFrame(
        {
            'x': [xydn[0][n] for n in range(1, npn + 1)],
            'y': [xydn[1][n] for n in range(1, npn + 1)],
            'z': [xydn[2][n] for n in range(1, npn + 1)]
        },
        index=range(1, npn + 1)
    )
    for m in range(len(submesh.elements)):
        new_m = int(mmaps[m])
        nmn[new_m][0] = mmaps[ submesh.elements.loc[m+1, 'node_1'] ]
        nmn[new_m][1] = mmaps[ submesh.elements.loc[m+1, 'node_2'] ]
        nmn[new_m][2] = mmaps[ submesh.elements.loc[m+1, 'node_3'] ]
    for m in range(num_elements):
        if netable[m] == 0:
            new_m = mmap[m]

            nmn[new_m][0] = mmap[ full_mesh.elements.loc[m+1, 'node_1'] ]
            nmn[new_m][1] = mmap[ full_mesh.elements.loc[m+1, 'node_2'] ]
            nmn[new_m][2] = mmap[ full_mesh.elements.loc[m+1, 'node_3'] ]
    elements_df = pd.DataFrame(
        {
            'node_1': [int(nmn[m][0]) for m in range(nen)],
            'node_2': [int(nmn[m][1]) for m in range(nen)],
            'node_3': [int(nmn[m][2]) for m in range(nen)],
        },
        index=range(1, nen + 1)
    )
    elements_df['n'] = 3
    elements_df = elements_df.iloc[:, [3, 0, 1, 2]]
   

    # Boundary reconstruction
    t0 = time.time()
    print("[Boundaries] Reconstructing open boundaries...")
    # preliminary_segments = []
    # for segment in full_mesh.boundaries['open']:
    #     current_segment_nodes = []
    #     for boundary_node in segment.nodes:
    #         old_node_id = boundary_node.node_id
    #         old_node_pos = global_node_id_to_pos[old_node_id]
    #         if nptable[old_node_pos] == 0:  # Node preserved in merged mesh
    #             new_node_id = int(nmap[old_node_pos])
    #             current_segment_nodes.append(BoundaryNode(node_id=new_node_id))
    #         else:  # Node replaced, break segment
    #             if current_segment_nodes:
    #                 preliminary_segments.append(current_segment_nodes)
    #                 current_segment_nodes = []
    #     if current_segment_nodes:
    #         preliminary_segments.append(current_segment_nodes)

    #     # Phase 2: Reconnect segments by matching tails to heads
    # merged_open_bounds = []
    # used_segments = set()

    # while len(used_segments) < len(preliminary_segments):
    #     # Find an unused segment with no predecessor (head of chain)
    #     chain_start_idx = None
    #     for n in range(len(preliminary_segments)):
    #         if n in used_segments:
    #             continue
            
    #         segment_first_node = preliminary_segments[n][0].node_id
            
    #         # Check if any other unused segment ends with this node
    #         has_predecessor = False
    #         for m in range(len(preliminary_segments)):
    #             if m == n or m in used_segments:
    #                 continue
    #             segment_last_node = preliminary_segments[m][-1].node_id
    #             if segment_first_node == segment_last_node:
    #                 has_predecessor = True
    #                 break
            
    #         if not has_predecessor:
    #             chain_start_idx = n
    #             break
        
    #     if chain_start_idx is None:
    #         # This shouldn't happen with valid boundaries, but handle gracefully
    #         # Take any remaining unused segment
    #         for n in range(len(preliminary_segments)):
    #             if n not in used_segments:
    #                 chain_start_idx = n
    #                 break
        
    #     # Build the chain starting from chain_start_idx
    #     chain_nodes = list(preliminary_segments[chain_start_idx])
    #     used_segments.add(chain_start_idx)
    #     current_end_node = chain_nodes[-1].node_id
        
    #     # Keep extending the chain
    #     while True:
    #         found_next = False
    #         for n in range(len(preliminary_segments)):
    #             if n in used_segments:
    #                 continue
                
    #             segment_first_node = preliminary_segments[n][0].node_id
    #             if current_end_node == segment_first_node:
    #                 # Found continuation, append nodes (skip first since it duplicates last)
    #                 chain_nodes.extend(preliminary_segments[n][1:])
    #                 used_segments.add(n)
    #                 current_end_node = chain_nodes[-1].node_id
    #                 found_next = True
    #                 break
            
    #         if not found_next:
    #             break  # No more segments to connect
    # ============================================================
    # Boundary reconstruction (Fortran-faithful)
    # ============================================================

    preliminary_segments = []
    # Each entry:
    #   {
    #     'nodes': List[BoundaryNode],
    #     'ibtype': int | None,
    #     'polarity': +1 | -1
    #   }
    # ------------------------------------------------------------
    # Phase 1: split boundaries into preliminary segments
    # ------------------------------------------------------------
    for boundary_name, segments in full_mesh.boundaries.items():

        for segment in segments:

            current_nodes = []
            polarity = 1
            ibtype = segment.ibtype  # None for open, int for land

            for j, bnode in enumerate(segment.nodes):

                old_id = bnode.node_id
                old_pos = global_node_id_to_pos[old_id]

                # if( nptable(nbvv(k,j)) == 1 ) cycle
                if nptable[old_pos] == 1:
                    if current_nodes:
                        preliminary_segments.append({
                            'nodes': current_nodes,
                            'ibtype': ibtype,
                            'polarity': polarity
                        })
                        current_nodes = []
                        polarity = 1
                    continue

                new_node_id = int(nmap[old_pos])

                # ---- connected_to handling (ibtype 4/24 etc) ----
                new_connected = None
                if bnode.connected_to is not None:
                    old_conn_pos = global_node_id_to_pos[bnode.connected_to]
                    if nptable[old_conn_pos] == 1:
                        if current_nodes:
                            preliminary_segments.append({
                                'nodes': current_nodes,
                                'ibtype': ibtype,
                                'polarity': polarity
                            })
                            current_nodes = []
                            polarity = 1
                        continue
                    new_connected = int(nmap[old_conn_pos])

                current_nodes.append(
                    BoundaryNode(
                        node_id=new_node_id,
                        connected_to=new_connected,
                        barrier_heights=bnode.barrier_heights
                    )
                )

                # ---- polarity detection (Fortran nplseg = -1) ----
                if j + 1 < len(segment.nodes):
                    next_old_id = segment.nodes[j + 1].node_id
                    next_old_pos = global_node_id_to_pos[next_old_id]

                    if nptable[next_old_pos] != 1:
                        n1 = old_id
                        n2 = next_old_id
                        for m in range(nbl):
                            m1 = nblnc[0, m]
                            m2 = nblnc[1, m]
                            if (n1 == m2) and (n2 == m1):
                                polarity = -1
                                break

            if current_nodes:
                preliminary_segments.append({
                    'nodes': current_nodes,
                    'ibtype': ibtype,
                    'polarity': polarity
                })
    
    # ------------------------------------------------------------
    # Phase 2: chain preliminary segments
    # ------------------------------------------------------------
    merged_open_bounds = []
    merged_land_bounds = []

    used = [0] * len(preliminary_segments)

    while sum(used) < len(preliminary_segments):

        # ---- find head ----
        head_idx = None
        # ---- find head ----
        for n, seg in enumerate(preliminary_segments):
            if used[n]:
                continue

            nsta = seg['nodes'][0]
            found_predecessor = False

            for m, other in enumerate(preliminary_segments):
                if m == n or used[m]:
                    continue
                if seg['ibtype'] != other['ibtype']:
                    continue
                    continue

                nend = other['nodes'][-1]

                io = ii = 0
                if (nsta.node_id == nend.node_id and
                    nsta.connected_to == nend.connected_to):
                    io = 1
                if (nsta.node_id == nend.connected_to and
                    nsta.connected_to == nend.node_id):
                    ii = 1

                if io + ii == 2:
                    raise RuntimeError("Fatal boundary ambiguity (Basis104)")

                if io + ii == 1:
                    found_predecessor = True
                    break

            if not found_predecessor:
                head_idx = n
                break

        if head_idx is None:
            head_idx = used.index(0)  # Fortran-safe fallback

        # ---- start chain ----
        seg = preliminary_segments[head_idx]
        chain_nodes = list(seg['nodes'])
        chain_ibtype = seg['ibtype']
        used[head_idx] = 1

        nend = chain_nodes[-1]

        # ---- extend chain ----
        extended = True
        while extended:
            extended = False
            for n, seg in enumerate(preliminary_segments):
                if used[n] or seg['ibtype'] != chain_ibtype:
                    continue

                first = seg['nodes'][0]

                io = ii = 0
                if (nend.node_id == first.node_id and
                    nend.connected_to == first.connected_to):
                    io = 1
                if (nend.node_id == first.connected_to and
                    nend.connected_to == first.node_id):
                    ii = 1

                if io + ii == 2:
                    raise RuntimeError("Fatal boundary ambiguity (Basis105)")

                if io + ii == 1:
                    if io == 1:
                        chain_nodes.extend(seg['nodes'][1:])
                    else:
                        # swapped orientation
                        for node in seg['nodes'][1:]:
                            chain_nodes.append(
                                BoundaryNode(
                                    node_id=node.connected_to,
                                    connected_to=node.node_id,
                                    barrier_heights=node.barrier_heights
                                )
                            )

                    nend = chain_nodes[-1]
                    used[n] = 1
                    extended = True
                    break

        # ---- append EXACTLY ONCE per chain ----
        if chain_ibtype is None:
            merged_open_bounds.append(
                BoundarySegment(
                    nodes=tuple(chain_nodes),
                    category='open',
                    ibtype=None
                )
            )
        else:
            merged_land_bounds.append(
                BoundarySegment(
                    nodes=tuple(chain_nodes),
                    category='land',
                    ibtype=chain_ibtype
                )
            )

    # ------------------------------------------------------------
    # Final result (drop-in compatible)
    # ------------------------------------------------------------
    # merged_mesh.boundaries = {
    #     'open': merged_open_bounds,
    #     'land': merged_land_bounds
    # }



        # # Create the final boundary segment from the complete chain
        # merged_open_bounds.append(BoundarySegment(
        #     nodes=tuple(chain_nodes),
        #     category='open',
        #     ibtype=None
        # ))
    # print(f"             Reconnected {len(preliminary_segments)} fragments into {len(merged_open_bounds)} boundaries ({time.time()-t0:.1f}s)\n")

    # t0 = time.time()
    # print("[Boundaries] Reconstructing land boundaries...", len(full_mesh.boundaries['land']))
    # merged_land_bounds = reconstruct_flow_boundaries(full_mesh.boundaries['land'], nmap, nptable, global_node_id_to_pos)
    # print(f"             Created {len(merged_land_bounds)} land boundaries ({time.time()-t0:.1f}s)\n")

    # t0 = time.time()
    # print("[Final] Creating merged mesh object and writing to file...")
    merged_mesh = Mesh(nodes=new_nodes_df, elements=elements_df, boundaries={'open': merged_open_bounds, 'land': merged_land_bounds})
    
    # merged_mesh.nodes = merged_nodes
    # merged_mesh.elements = merged_ele
    # merged_mesh.boundaries = {'open': merged_open_bounds, 'land': merged_land_bounds}
    write_fort14('merged.14', merged_mesh)
    print(f"        Wrote merged.14 ({time.time()-t0:.1f}s)")

    print("\n" + "="*80)
    print("MERGE COMPLETE")
    print("="*80)
    
            
def reconstruct_flow_boundaries(global_land_boundaries, nmap, nptable, global_node_id_to_pos):
    """Reconstruct land/flow boundaries for merged mesh with segment reconnection.

    Args:
        global_land_boundaries (List): lists of boundary segments
        nmap: node_id_mapping (global -> merged), indexed by position
        nptable: node classification (0=preserved, 1=replaced, 2=border), indexed by position
        global_node_id_to_pos: dict mapping node IDs to array positions
    """
    from submesh.boundaries import IBTYPE_SIMPLE, IBTYPE_BARRIER, IBTYPE_CONNECTED, BoundaryNode, BoundarySegment
    
    # Phase 1: Break into preliminary segments grouped by ibtype
    preliminary_segments = []  # List of dicts with 'nodes', 'ibtype'
    
    for segment in global_land_boundaries:
        current_segment_nodes = []
        ibtype = segment.ibtype
        
        for boundary_node in segment.nodes:
            old_node_id = boundary_node.node_id
            old_node_pos = global_node_id_to_pos[old_node_id]

            if nptable[old_node_pos] == 0:  # Node preserved
                new_node_id = int(nmap[old_node_pos])
                new_connected_to = None

                # Handle connected_to remapping
                if boundary_node.connected_to is not None:
                    old_connected_id = boundary_node.connected_to
                    old_connected_pos = global_node_id_to_pos[old_connected_id]
                    if nptable[old_connected_pos] == 0:
                        new_connected_to = int(nmap[old_connected_pos])
                    else:
                        # Connected node was replaced, break segment
                        if current_segment_nodes:
                            preliminary_segments.append({
                                'nodes': current_segment_nodes,
                                'ibtype': ibtype
                            })
                            current_segment_nodes = []
                        continue
                
                new_node = BoundaryNode(
                    node_id=new_node_id,
                    connected_to=new_connected_to,
                    barrier_heights=boundary_node.barrier_heights
                )
                current_segment_nodes.append(new_node)
            else:
                # Node replaced, break segment
                if current_segment_nodes:
                    preliminary_segments.append({
                        'nodes': current_segment_nodes,
                        'ibtype': ibtype
                    })
                    current_segment_nodes = []
        
        if current_segment_nodes:
            preliminary_segments.append({
                'nodes': current_segment_nodes,
                'ibtype': ibtype
            })
    
    # Phase 2: Reconnect segments by matching tails to heads (same ibtype only)
    merged_land_bounds = []
    used_segments = set()
    
    while len(used_segments) < len(preliminary_segments):
        # Find an unused segment with no predecessor (head of chain)
        chain_start_idx = None
        for n in range(len(preliminary_segments)):
            if n in used_segments:
                continue
            
            seg_n = preliminary_segments[n]
            first_node = seg_n['nodes'][0]
            ibtype_n = seg_n['ibtype']
            
            # Check if any other unused segment (same ibtype) ends with this node
            has_predecessor = False
            for m in range(len(preliminary_segments)):
                if m == n or m in used_segments:
                    continue
                if preliminary_segments[m]['ibtype'] != ibtype_n:
                    continue
                
                seg_m = preliminary_segments[m]
                last_node = seg_m['nodes'][-1]
                
                # Check if nodes match (considering connected_to swapping)
                if ibtype_n in IBTYPE_CONNECTED:
                    # For connected boundaries, check both orientations
                    match_normal = (first_node.node_id == last_node.node_id and 
                                   first_node.connected_to == last_node.connected_to)
                    match_swapped = (first_node.node_id == last_node.connected_to and 
                                    first_node.connected_to == last_node.node_id)
                    if match_normal or match_swapped:
                        has_predecessor = True
                        break
                else:
                    # Simple or barrier boundaries
                    if first_node.node_id == last_node.node_id:
                        has_predecessor = True
                        break
            
            if not has_predecessor:
                chain_start_idx = n
                break
        
        if chain_start_idx is None:
            # Take any remaining unused segment
            for n in range(len(preliminary_segments)):
                if n not in used_segments:
                    chain_start_idx = n
                    break
        
        # Build the chain starting from chain_start_idx
        chain_segment = preliminary_segments[chain_start_idx]
        chain_nodes = list(chain_segment['nodes'])
        chain_ibtype = chain_segment['ibtype']
        used_segments.add(chain_start_idx)
        current_end_node = chain_nodes[-1]
        
        # Keep extending the chain
        while True:
            found_next = False
            for n in range(len(preliminary_segments)):
                if n in used_segments:
                    continue
                
                seg_n = preliminary_segments[n]
                if seg_n['ibtype'] != chain_ibtype:
                    continue
                
                first_node = seg_n['nodes'][0]
                
                # Check for connection (considering orientation)
                normal_match = False
                swapped_match = False
                
                if chain_ibtype in IBTYPE_CONNECTED:
                    # Check both orientations for connected boundaries
                    normal_match = (current_end_node.node_id == first_node.node_id and
                                   current_end_node.connected_to == first_node.connected_to)
                    swapped_match = (current_end_node.node_id == first_node.connected_to and
                                    current_end_node.connected_to == first_node.node_id)
                else:
                    normal_match = (current_end_node.node_id == first_node.node_id)
                
                if normal_match:
                    # Normal orientation: append nodes, skip first (duplicate)
                    chain_nodes.extend(seg_n['nodes'][1:])
                    used_segments.add(n)
                    current_end_node = chain_nodes[-1]
                    found_next = True
                    break
                elif swapped_match:
                    # Swapped orientation: swap node_id and connected_to for all nodes
                    for node in seg_n['nodes'][1:]:  # Skip first (duplicate)
                        swapped_node = BoundaryNode(
                            node_id=node.connected_to,
                            connected_to=node.node_id,
                            barrier_heights=node.barrier_heights
                        )
                        chain_nodes.append(swapped_node)
                    used_segments.add(n)
                    current_end_node = chain_nodes[-1]
                    found_next = True
                    break
            
            if not found_next:
                break  # No more segments to connect
        
        # Create the final boundary segment from the complete chain
        merged_land_bounds.append(BoundarySegment(
            nodes=tuple(chain_nodes),
            category='land',
            ibtype=chain_ibtype
        ))
    
    return merged_land_bounds


def create_netable_from_node_mapping(full_mesh, node_id_map_path='node_id_map.csv'):
    """Create element mask using node ID mapping.
    Args:
        global_mesh: Mesh object with .elements DataFrame
        node_id_map_path: Path to csv with columns (mesh_node_id, subset_node_id)
        
    Returns:
        netable: numpy array indexed by global element ID (0 or 1)
    """
    # Load the node_mapping
    node_map_df = pd.read_csv(node_id_map_path)
    
    # Get set of global node_id that are in the subgrid
    subgrid_global_nodes = set(node_map_df['mesh_node_id'])

    # Initialize netable (all zeros = not in subgrid)
    num_elements = len(full_mesh.elements)
    netable = np.zeros(num_elements, dtype=np.int32)

    # For each element, check if ALL its nodes are in teh subgrid
    for idx, elem_row in full_mesh.elements.iterrows():
        elem_nodes = [elem_row['node_1'], elem_row['node_2'], elem_row['node_3']]

        if all(node_id in subgrid_global_nodes for node_id in elem_nodes):
            netable[idx] = 1
    return netable

def mkeline(num_nodes, num_elements, nodes_per_element, netable, nbn_max):
    """Direct translation of Fortran mkeline subroutine.
    
    Finds boundary edges of the subgrid region by:
    1. Building node-to-elements adjacency (all elements)
    2. For each edge in subgrid elements, checking if the neighboring 
       element sharing that edge is also in the subgrid
    """
    import time
    print("Active elements:", np.sum(netable == 1))
    print("Inactive elements:", np.sum(netable != 1))

    t_start = time.time()
    nm = nodes_per_element[['node_1', 'node_2', 'node_3']]
    # Convert to numpy array (0-indexed)
    if hasattr(nodes_per_element, 'values'):
        nm = nm.values.astype(np.int32)-1 #nodes_per_element.values.astype(np.int32)
    else:
        nm = np.asarray(nm, dtype=np.int32)-1
   
    netable = np.asarray(netable, dtype=np.int32)
    
    # Count active elements
    n_active = int(np.sum(netable == 1))
    print(f"              mkeline: Processing {n_active:,}/{num_elements:,} elements")
    
    # Search Elements around a node (Fortran lines 1785-1809)
    # This builds nean(node, :) = list of elements containing that node
    nean = [[] for _ in range(num_nodes)]
    for m in range(num_elements):
        for i in range(3):
            n = nm[m][i]
            nean[n].append(m)
            
    
    # Search the border of global grid and sub-grid (Fortran lines 1811-1840)
    nbl = 0
    nblnc_list = []
    
    for m in range(num_elements):
        if netable[m] != 1:
            continue
        
        # Check each of the 3 edges (Fortran line 1815)
        for i in range(3): 
            # Call eline to get edge nodes (Fortran line 1816)
            n1, n2 = eline(i, nm, m)
            is_share_active = False
            # Check if this edge has a neighbor in the subgrid (Fortran lines 1817-1832)
            icheck = 0
            share_elements = set(nean[n1]).intersection(nean[n2])
            for neighbor in share_elements:
                if neighbor != m:
                    if netable[neighbor] == 1:
                        is_share_active = True
                        break
            if not is_share_active:
                nblnc_list.append([n1, n2])
            
    
    print(f"              mkeline: Found {len(nblnc_list),} boundary edges ({time.time()-t_start:.1f}s total)")
    
    nbl = len(nblnc_list)
    nblnc = np.array(nblnc_list).T
    
    return nbl, nblnc


def eline(i, nm, m):
    """
    Returns the two nodes (n1, n2) defining edge 'i' of element 'm'.
    
    Parameters:
    i  : Edge index (1, 2, or 3) from your loop
    nm : The connectivity matrix (num_elements x 3)
    m  : The current element index
    """
    if i == 0:
        # Edge 1: First and Second node (indices 0 and 1)
        n1 = nm[m, 0]
        n2 = nm[m, 1]
    elif i == 1:
        # Edge 2: Second and Third node (indices 1 and 2)
        n1 = nm[m, 1]
        n2 = nm[m, 2]
    elif i == 2:
        # Edge 3: Third and First node (indices 2 and 0)
        n1 = nm[m, 2]
        n2 = nm[m, 0]
    else:
        raise ValueError(f"Invalid edge index i={i}. Expected 1, 2, or 3.")
        
    return n1, n2

 # from scipy.spatial import cKDTree

    # t0 = time.time()
    # print("[Step 8/10] Matching boundary nodes between global and submesh...")

    # # Create position-to-ID and ID-to-position mappings for both meshes
    # # Fortran uses 1-indexed positional arrays, Python uses node IDs as DataFrame index
    # global_node_ids = full_mesh.nodes.index.to_numpy()
    # submesh_node_ids = submesh.nodes.index.to_numpy()
    # global_elem_ids = full_mesh.elements.index.to_numpy()

    # # Create reverse mappings: ID -> position (0-indexed)
    # global_node_id_to_pos = {node_id: pos for pos, node_id in enumerate(global_node_ids)}
    # submesh_node_id_to_pos = {node_id: pos for pos, node_id in enumerate(submesh_node_ids)}
    # global_elem_id_to_pos = {elem_id: pos for pos, elem_id in enumerate(global_elem_ids)}

    # # Arrays indexed by POSITION (Fortran-style, but 0-indexed in Python)
    # nmap = np.zeros(len(full_mesh.nodes), dtype=int)
    # nmaps = np.ones(len(submesh.nodes), dtype=int) * -1
    # mmap = np.zeros(len(full_mesh.elements), dtype=np.int32)
    # nprop1 = np.zeros(len(full_mesh.nodes), dtype=int)
    # nptable = np.zeros(len(full_mesh.nodes), dtype=int)

    # # Mark all nodes in subgrid elements with nptable = 1 (to be replaced by submesh)
    # for elem_pos in range(len(full_mesh.elements)):
    #     if netable[elem_pos] == 1:
    #         elem_id = global_elem_ids[elem_pos]
    #         n1_id = full_mesh.elements.loc[elem_id, 'node_1']
    #         n2_id = full_mesh.elements.loc[elem_id, 'node_2']
    #         n3_id = full_mesh.elements.loc[elem_id, 'node_3']
    #         nptable[global_node_id_to_pos[n1_id]] = 1
    #         nptable[global_node_id_to_pos[n2_id]] = 1
    #         nptable[global_node_id_to_pos[n3_id]] = 1

    # eps = 1e-6

    # # --- Collect unique boundary nodes ---
    # global_bnodes = np.unique(nbdnc)
    # sub_bnodes = np.unique(nblncs)

    # # --- Coordinates ---
    # global_xy = full_mesh.nodes.loc[global_bnodes, ['x', 'y']].to_numpy()
    # sub_xy    = submesh.nodes.loc[sub_bnodes,    ['x', 'y']].to_numpy()

    # # --- KD-tree ---
    # tree = cKDTree(sub_xy)
    # dist, idx = tree.query(global_xy, distance_upper_bound=eps)

    # matched = dist <= eps
    # n_matched = matched.sum()

    # # --- Populate mappings ---
    # for g_id, s_id in zip(global_bnodes[matched], sub_bnodes[idx[matched]]):
    #     s_pos = submesh_node_id_to_pos[s_id]
    #     g_pos = global_node_id_to_pos[g_id]
    #     nmaps[s_pos] = g_id
    #     nprop1[g_pos] = s_id
    #     nptable[g_pos] = 2

    # print(f"            Matched {n_matched:,} boundary nodes ({time.time()-t0:.1f}s)")

    # unmatched = len(global_bnodes) - n_matched
    # if unmatched > 0:
    #     print(f"            WARNING: {unmatched} boundary nodes failed to match\n")
    # else:
    #     print()


    # def mkeline(num_nodes, num_elements, nodes_per_element, netable):
    # Search for elements around a node
    # icheck = 0
    # numean = np.zeros(num_nodes, dtype=int)
    
    # nean_max = 10
    
    # while True:
    #     nean = np.zeros((num_nodes, nean_max), dtype=int)
    #     numean[:] = 0
    #     restart = False
    #     for m in range(num_elements):
    #         for i in range(3):
    #             n = nodes_per_element[m,i]
    #             numean[n] += 1
    #             if numean[n] > nean_max:
    #                 icheck += 1
    #                 if icheck == 1:
    #                     print('Salvere000: System require much more nean_max')
    #                 nean_max += 5
    #                 print(f'System try to recalculate with nean_max = {nean_max}')
    #                 restart = True
    #                 break
    #             nean[n, numean[n]-1] = m
    #         if restart:
    #             break
    #     if not restart:
    #         break
    # if icheck != 0:
    #     print('SUccess!!!!!!!!!!')

    # # Search the border of global grid and subgrid
    # nbl = 0
    # for m in range(num_elements):
    #     if (netable[m] == 1):
    #         for i in range(3):
    #             n1, n2 = eline(i, nodes_per_element, m, num_elements)
    #             icheck = 0
    #             for i1 in range(numean[n1]):
    #                 m1 = nean[n1, i1]
    #                 for i2 in range(numean[n2]):
    #                     m2 = nean[n2, i2]
    #                     if m1 == m2:
    #                         if m1 != m:
    #                             icheck += 1
    #                             break
    #             if (icheck != 0):
    #                 break
    #         if (icheck == 0):
    #             nblnc[0, nbl] = n1
    #             nblnc[1, nbl] = n2
    #             nbl += 1
    # return nbl, nblnc
# def mkeline(num_nodes, num_elements, nodes_per_element, netable):
    # """Find boundary edges using dictionary-based adjacency."""
    # from collections import defaultdict
    # import time

    # t_start = time.time()

    # # Convert DataFrame to numpy array ONCE for fast indexing
    # if hasattr(nodes_per_element, 'values'):
    #     elem_array = nodes_per_element.values.astype(np.int32)
    # else:
    #     elem_array = nodes_per_element

    # # Count how many elements are active
    # n_active = int(np.sum(netable))
    # print(f"              mkeline: Processing {n_active:,}/{num_elements:,} elements")

    # # Build node-to-elements mapping
    # t0 = time.time()
    # node_to_elems = defaultdict(list)
    # for elem_idx in range(num_elements):
    #     for i in range(3):
    #         node_id = elem_array[elem_idx, i]
    #         node_to_elems[node_id].append(elem_idx)
    # print(f"              mkeline: Built node-element map ({time.time()-t0:.1f}s)")

    # # Find boundary edges
    # t0 = time.time()
    # boundary_edges = []
    # edge_patterns = [(0, 1), (1, 2), (2, 0)]

    # # Progress tracking
    # progress_interval = max(1, num_elements // 20)  # Report every 5%
    # checked = 0

    # for elem_idx in range(num_elements):
    #     if elem_idx % progress_interval == 0 and elem_idx > 0:
    #         pct = 100 * elem_idx / num_elements
    #         print(f"              mkeline: ... {pct:.0f}% ({elem_idx:,}/{num_elements:,}, found {len(boundary_edges):,} edges)")

    #     if netable[elem_idx] != 1:
    #         continue

    #     checked += 1
    #     elem_nodes = elem_array[elem_idx]
    
    #     for i, j in edge_patterns:
    #         n1, n2 = elem_nodes[i], elem_nodes[j]

    #         # Find elements containing both nodes
    #         elems_n1 = set(node_to_elems[n1])
    #         elems_n2 = set(node_to_elems[n2])
    #         shared_elems = elems_n1 & elems_n2
    #         shared_elems.discard(elem_idx)
    #         # ✅ FILTER to only subgrid elements
    #         shared_elems = {e for e in shared_elems if netable[e] == 1}
    #         # Check if neighbor is in region
    #         has_neighbor = len(shared_elems) > 0

    #         if not has_neighbor:
    #             boundary_edges.append((n1, n2))
    #             break

    # print(f"              mkeline: Found {len(boundary_edges):,} boundary edges from {checked:,} elements ({time.time()-t_start:.1f}s total)")

    # # ADD THIS CHECK:
    # node_3_edges = [edge for edge in boundary_edges if 3 in edge]
    # if node_3_edges:
    #     print(f"              WARNING: Found {len(node_3_edges)} edges containing node 3: {node_3_edges}")
    
    # nbl = len(boundary_edges)
    # nblnc = np.array(boundary_edges, dtype=np.int32).T if boundary_edges else np.zeros((2, 0), dtype=np.int32)

    # return nbl, nblnc
# def mkeline(num_nodes, num_elements, nodes_per_element, netable):
#     """Ultra-fast boundary edge detection using pandas."""
#     import time
#     import pandas as pd
    
#     t_start = time.time()
    
#     # Convert to numpy arrays
#     if hasattr(nodes_per_element, 'values'):
#         elem_array = nodes_per_element.values.astype(np.int32)
#     else:
#         elem_array = np.asarray(nodes_per_element, dtype=np.int32)
    
#     netable = np.asarray(netable, dtype=np.int8)
    
#     # Filter to active elements
#     active_mask = (netable == 1)
#     active_elems = elem_array[active_mask]
#     n_active = len(active_elems)
    
#     print(f"              mkeline: Processing {n_active:,}/{num_elements:,} elements")
    
#     # Create all edges
#     t0 = time.time()
#     edges = np.empty((n_active * 3, 2), dtype=np.int32)
#     edges[0::3] = active_elems[:, [0, 1]]
#     edges[1::3] = active_elems[:, [1, 2]]
#     edges[2::3] = active_elems[:, [2, 0]]
    
#     # Sort each edge
#     edges.sort(axis=1)
    
#     print(f"              mkeline: Created edges ({time.time()-t0:.1f}s)")
    
#     # Use pandas for fast groupby
#     t0 = time.time()
#     df = pd.DataFrame(edges, columns=['n1', 'n2'])
#     edge_counts = df.groupby(['n1', 'n2']).size()
    
#     # Boundary edges appear once
#     boundary_edges_df = edge_counts[edge_counts == 1].reset_index()[['n1', 'n2']]
#     boundary_edges = boundary_edges_df.values
    
#     print(f"              mkeline: Found {len(boundary_edges):,} boundary edges ({time.time()-t_start:.1f}s total)")
    
#     nbl = len(boundary_edges)
#     nblnc = boundary_edges.T if len(boundary_edges) > 0 else np.zeros((2, 0), dtype=np.int32)
    
#     return nbl, nblnc

# def eline(i, nm, m):
    # if i == 2:
    #     n1 = nm[m, 0]
    #     n2 = nm[m, 1]
    #     return n1, n2
    # elif i == 0:
    #     n1 = nm[m,1]
    #     n2 = nm[m,2]
    #     return n1, n2
    # elif i == 1:
    #     n1 = nm[m, 2]
    #     n2 = nm[m, 0]
    #     return n1, n2
