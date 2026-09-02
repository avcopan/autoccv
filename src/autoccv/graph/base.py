"""Generic graph engine.

Uses NetworkX for graph representation, with Pydantic node and edge data
validation. Node keys are arbitrary hashables (typically ints); node and edge
data are validated against user-defined :class:`Node` and :class:`Edge` models.

This module is deliberately self-contained: it depends only on ``networkx`` and
``pydantic`` and has no intra-package imports, so it can be copied into a new
project as a starting point for graph algorithms.
"""

import copy
import enum
from collections.abc import Collection, Iterator
from typing import Any, TypeVar

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from pydantic import BaseModel

EdgeKey = tuple[int, int]


class FieldEnumMismatchError(TypeError):
    """A model's nested ``Field`` enum does not match its declared fields."""


class FieldEnumModel(BaseModel):
    """Base model whose field names can be mirrored by a nested ``Field`` enum.

    A subclass that needs its field names as string constants (e.g. as NetworkX
    attribute keys) declares a nested ``class Field(enum.StrEnum)``. Its members
    are checked for exact parity with the model fields at class-creation time, so
    the enum can never silently drift from the schema. Because ``StrEnum``
    members are real ``str`` instances, ``Atom.Field.symbol`` can be used
    directly as a key -- no ``.value`` needed.
    """

    @classmethod
    def __pydantic_init_subclass__(cls, **kwargs: object) -> None:
        """Check any nested ``Field`` enum for parity with the model fields."""
        super().__pydantic_init_subclass__(**kwargs)

        field_enum = getattr(cls, "Field", None)
        if not (isinstance(field_enum, type) and issubclass(field_enum, enum.Enum)):
            return

        enum_names = {member.value for member in field_enum}
        model_names = set(cls.model_fields)
        if enum_names != model_names:
            missing = sorted(model_names - enum_names)
            extra = sorted(enum_names - model_names)
            parts = []
            if missing:
                parts.append(f"missing from Field: {missing}")
            if extra:
                parts.append(f"not a model field: {extra}")
            detail = "; ".join(parts)
            msg = f"{cls.__name__}.Field is out of sync with its fields ({detail})"
            raise FieldEnumMismatchError(msg)


class Node(FieldEnumModel):
    """Base class for node data. Subclass to add validated fields."""


class Edge(FieldEnumModel):
    """Base class for edge data. Subclass to add validated fields."""


NodeT = TypeVar("NodeT", bound=Node)
EdgeT = TypeVar("EdgeT", bound=Edge)


class Graph[NodeT: Node, EdgeT: Edge](nx.Graph):
    """Generic graph with typed node and edge data."""

    node_type: type[NodeT]
    edge_type: type[EdgeT]

    def validate(self) -> None:
        """Validate node and edge data against the node and edge classes."""
        for *_, data in self.nodes(data=True):
            self.node_type.model_validate(data)

        for *_, data in self.edges(data=True):
            self.edge_type.model_validate(data)


# Properties
def node_keys[NodeT: Node, EdgeT: Edge](gra: Graph[NodeT, EdgeT]) -> list[int]:
    """Get the list of node keys.

    Args:
        gra: A graph.

    Returns:
        The node keys.
    """
    return list(gra.nodes())


# Transformations
def remove_edges[NodeT: Node, EdgeT: Edge](
    gra: Graph[NodeT, EdgeT],
    edges: Collection[EdgeKey],
    *,
    in_place: bool = False,
) -> Graph[NodeT, EdgeT]:
    """Remove edges from a graph.

    Args:
        gra: A graph.
        edges: The edges to remove.
        in_place: Whether to modify the graph in place instead of returning a copy.

    Returns:
        The graph with the edges removed.
    """
    gra = gra if in_place else copy.deepcopy(gra)
    gra.remove_edges_from(edges)
    return gra


# Algorithms
def graph_matcher[NodeT: Node, EdgeT: Edge](
    gra1: Graph[NodeT, EdgeT], gra2: Graph[NodeT, EdgeT]
) -> GraphMatcher:
    """Get a NetworkX graph matcher comparing node and edge data fields.

    Args:
        gra1: The first graph.
        gra2: The second graph.

    Returns:
        A graph matcher that matches nodes and edges on all data fields.
    """
    node_fields = gra1.node_type.model_fields.keys()
    edge_fields = gra1.edge_type.model_fields.keys()

    def node_match(n1: dict[str, Any], n2: dict[str, Any]) -> bool:
        return all(n1[field] == n2[field] for field in node_fields)

    def edge_match(e1: dict[str, Any], e2: dict[str, Any]) -> bool:
        return all(e1[field] == e2[field] for field in edge_fields)

    return GraphMatcher(gra1, gra2, node_match=node_match, edge_match=edge_match)


def isomorphisms[NodeT: Node, EdgeT: Edge](
    gra1: Graph[NodeT, EdgeT], gra2: Graph[NodeT, EdgeT]
) -> Iterator[dict[int, int]]:
    """Iterate over isomorphisms between two graphs.

    Args:
        gra1: The first graph.
        gra2: The second graph.

    Yields:
        Node mappings from the first graph to the second.
    """
    return graph_matcher(gra1, gra2).isomorphisms_iter()


def isomorphism[NodeT: Node, EdgeT: Edge](
    gra1: Graph[NodeT, EdgeT], gra2: Graph[NodeT, EdgeT]
) -> dict[int, int] | None:
    """Get an isomorphism between two graphs, if there is one.

    Args:
        gra1: The first graph.
        gra2: The second graph.

    Returns:
        A node mapping from the first graph to the second, or ``None``.
    """
    return next(isomorphisms(gra1, gra2), None)


# Comparisons
def is_isomorphic[NodeT: Node, EdgeT: Edge](
    gra1: Graph[NodeT, EdgeT], gra2: Graph[NodeT, EdgeT]
) -> bool:
    """Determine whether two graphs are isomorphic.

    Args:
        gra1: The first graph.
        gra2: The second graph.

    Returns:
        Whether the graphs are isomorphic.
    """
    return graph_matcher(gra1, gra2).is_isomorphic()
