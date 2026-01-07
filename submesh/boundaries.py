"""ADCIRC boundary data structures.

This module provides dataclass-based representations for ADCIRC mesh boundaries,
replacing the previous dictionary-based structure with parallel lists.
"""

from dataclasses import dataclass
from typing import Optional, List, Literal

# ADCIRC boundary type constants
IBTYPE_SIMPLE = frozenset({0, 1, 2, 10, 11, 12, 20, 21, 22, 30, 52})
IBTYPE_BARRIER = frozenset({3, 13, 23})
IBTYPE_CONNECTED = frozenset({4, 24, 64})

BoundaryCategory = Literal['open', 'land']


@dataclass(frozen=True)
class BoundaryNode:
    """A single node in a boundary segment with all its attributes.

    Attributes:
        node_id: The mesh node ID
        connected_to: For IBTYPE_CONNECTED boundaries (4, 24, 64)
        barrier_heights: For IBTYPE_BARRIER (3,13,23) as 2-tuple,
                        or IBTYPE_CONNECTED (4,24,64) as 3-tuple
    """
    node_id: int
    connected_to: Optional[int] = None
    barrier_heights: Optional[tuple[float, ...]] = None

    def __post_init__(self):
        """Validate barrier_heights tuple length."""
        if self.barrier_heights is not None:
            length = len(self.barrier_heights)
            if self.connected_to is not None:
                # Connected boundary: expect 3 values
                if length != 3:
                    raise ValueError(
                        f"Connected boundary node {self.node_id} requires "
                        f"3 barrier values, got {length}"
                    )
            else:
                # Barrier boundary: expect 2 values
                if length != 2:
                    raise ValueError(
                        f"Barrier boundary node {self.node_id} requires "
                        f"2 barrier values, got {length}"
                    )

    def is_simple(self) -> bool:
        """Check if this is a simple node (no extra attributes)."""
        return self.connected_to is None and self.barrier_heights is None

    def is_barrier(self) -> bool:
        """Check if this has barrier heights but no connection."""
        return self.connected_to is None and self.barrier_heights is not None

    def is_connected(self) -> bool:
        """Check if this is a connected boundary node."""
        return self.connected_to is not None


@dataclass(frozen=True)
class BoundarySegment:
    """A contiguous boundary segment consisting of ordered nodes.

    Attributes:
        nodes: Ordered tuple of BoundaryNode objects
        category: 'open' or 'land'
        ibtype: ADCIRC boundary type code (None for open boundaries)
    """
    nodes: tuple[BoundaryNode, ...]
    category: BoundaryCategory
    ibtype: Optional[int] = None

    def __post_init__(self):
        """Validate segment consistency."""
        if not self.nodes:
            raise ValueError("BoundarySegment must have at least one node")

        if self.category == 'land' and self.ibtype is None:
            raise ValueError("Land boundaries must have an ibtype")

        if self.category == 'open' and self.ibtype is not None:
            raise ValueError("Open boundaries should not have an ibtype")

        # Validate all nodes are consistent with ibtype
        if self.ibtype is not None:
            self._validate_nodes_match_ibtype()

    def _validate_nodes_match_ibtype(self):
        """Ensure all nodes match the declared ibtype."""
        if self.ibtype in IBTYPE_SIMPLE:
            if not all(n.is_simple() for n in self.nodes):
                raise ValueError(
                    f"ibtype {self.ibtype} (simple) but some nodes have attributes"
                )
        elif self.ibtype in IBTYPE_BARRIER:
            if not all(n.is_barrier() for n in self.nodes):
                raise ValueError(
                    f"ibtype {self.ibtype} (barrier) but some nodes lack barriers"
                )
        elif self.ibtype in IBTYPE_CONNECTED:
            if not all(n.is_connected() for n in self.nodes):
                raise ValueError(
                    f"ibtype {self.ibtype} (connected) but some nodes not connected"
                )

    def node_ids(self) -> List[int]:
        """Extract just the node IDs in order."""
        return [n.node_id for n in self.nodes]

    def remap(self, node_id_map: dict[int, int]) -> Optional['BoundarySegment']:
        """Remap node IDs using the provided mapping.

        Args:
            node_id_map: Dictionary mapping old node IDs to new node IDs

        Returns:
            New BoundarySegment with remapped nodes, or None if no nodes remain
            after remapping.
        """
        remapped_nodes = []

        for node in self.nodes:
            new_node_id = node_id_map.get(node.node_id)
            if new_node_id is None:
                continue  # Node not in submesh, skip

            # Handle connected_to remapping
            new_connected_to = None
            if node.connected_to is not None:
                new_connected_to = node_id_map.get(node.connected_to)
                # If connected node not in map, this boundary becomes invalid
                if new_connected_to is None:
                    continue

            remapped_nodes.append(BoundaryNode(
                node_id=new_node_id,
                connected_to=new_connected_to,
                barrier_heights=node.barrier_heights
            ))

        if not remapped_nodes:
            return None

        return BoundarySegment(
            nodes=tuple(remapped_nodes),
            category=self.category,
            ibtype=self.ibtype
        )


# Type alias for the top-level boundary structure
BoundaryDict = dict[BoundaryCategory, List[BoundarySegment]]
