"""Transition-state graphs.

A transition-state graph is a molecular graph whose edges carry an optional
:class:`Change` marking bonds that form or break over the course of a reaction.
"""

from enum import StrEnum

import networkx as nx
from rdkit.Chem import rdchem

from .base import EdgeKey, Graph
from .mol import Atom, Bond, MolGraph


class Change(StrEnum):
    """A bond change over the course of a reaction."""

    FORMED = "formed"
    BROKEN = "broken"
    FLEETING = "fleeting"


class TransBond(Bond):
    """Represents a bond in a transition-state graph."""

    class Field(StrEnum):
        """Field names of :class:`TransBond`, for use as graph attribute keys."""

        change = "change"

    change: Change | None

    def to_rdkit_bond_type(self) -> rdchem.BondType:
        """Convert to an RDKit Bond Type.

        Returns:
            The RDKit bond type.
        """
        if self.change is not None:
            return rdchem.BondType.HYDROGEN
        return rdchem.BondType.SINGLE


class TransGraph(Graph[Atom, TransBond]):
    """Transition-state graph."""

    node_type = Atom
    edge_type = TransBond


FORMED_BOND = TransBond(change=Change.FORMED)
BROKEN_BOND = TransBond(change=Change.BROKEN)


def from_bond_changes(gra: MolGraph, bond_changes: dict[EdgeKey, Change]) -> TransGraph:
    """Construct a transition-state graph from a graph and its bond changes.

    Args:
        gra: A molecular graph.
        bond_changes: The formed and broken bonds, keyed by bond.

    Returns:
        The transition-state graph.
    """
    ts_gra = TransGraph()
    ts_gra.add_nodes_from(gra.nodes(data=True))
    ts_gra.add_edges_from(gra.edges(), change=None)
    formed_bonds = {k for k, c in bond_changes.items() if c == Change.FORMED}
    broken_bonds = {k for k, c in bond_changes.items() if c == Change.BROKEN}
    ts_gra.add_edges_from(formed_bonds, change=Change.FORMED)
    ts_gra.add_edges_from(broken_bonds, change=Change.BROKEN)
    ts_gra.validate()
    return ts_gra


def bond_changes(gra: TransGraph) -> dict[EdgeKey, Change]:
    """Extract the formed and broken bonds from a transition-state graph.

    Args:
        gra: A transition-state graph.

    Returns:
        The bond changes, keyed by bond.
    """
    change = nx.get_edge_attributes(gra, TransBond.Field.change)
    return {k: v for k, v in change.items() if v is not None}


def formed_bonds(gra: TransGraph) -> set[EdgeKey]:
    """Extract the formed bonds from a transition-state graph.

    Args:
        gra: A transition-state graph.

    Returns:
        The formed bonds.
    """
    changes = bond_changes(gra)
    return {k for k, v in changes.items() if v == Change.FORMED}


def broken_bonds(gra: TransGraph) -> set[EdgeKey]:
    """Extract the broken bonds from a transition-state graph.

    Args:
        gra: A transition-state graph.

    Returns:
        The broken bonds.
    """
    changes = bond_changes(gra)
    return {k for k, v in changes.items() if v == Change.BROKEN}


def reverse(gra: TransGraph) -> TransGraph:
    """Reverse the direction of a transition-state graph.

    Args:
        gra: A transition-state graph.

    Returns:
        The reversed transition-state graph.
    """
    changes = bond_changes(gra)
    changes = {
        k: Change.FORMED if v == Change.BROKEN else Change.BROKEN
        for k, v in changes.items()
    }
    return from_bond_changes(products_graph(gra), changes)


def reactants_graph(gra: TransGraph) -> MolGraph:
    """Extract the reactant graph from a transition-state graph.

    Args:
        gra: A transition-state graph.

    Returns:
        The reactant graph.
    """
    rct_gra = MolGraph()
    rct_gra.add_nodes_from(gra.nodes(data=True))
    rct_gra.add_edges_from(gra.edges(data=True))
    rct_gra.remove_edges_from(formed_bonds(gra))
    return rct_gra


def products_graph(gra: TransGraph) -> MolGraph:
    """Extract the product graph from a transition-state graph.

    Args:
        gra: A transition-state graph.

    Returns:
        The product graph.
    """
    prd_gra = MolGraph()
    prd_gra.add_nodes_from(gra.nodes(data=True))
    prd_gra.add_edges_from(gra.edges(data=True))
    prd_gra.remove_edges_from(broken_bonds(gra))
    return prd_gra
