"""Molecular graphs.

Builds a chemistry layer on top of the generic :mod:`autoccv.graph.base` engine:
atoms as nodes, bonds as edges, and conversions to and from RDKit / SMILES /
InChI. Bond order information is excluded by design.
"""

from collections.abc import Mapping, Sequence

from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem.rdchem import Mol, RWMol

from ..element import bonding_capacity
from .base import Edge, Graph, Node, node_keys


# RDKit helpers (inlined to avoid an extra dependency)
def _mol_from_smiles(smi: str) -> Mol:
    """Get an RDKit molecule with explicit hydrogens from a SMILES string."""
    return Chem.AddHs(Chem.MolFromSmiles(smi))


def _mol_from_inchi(chi: str) -> Mol:
    """Get an RDKit molecule with explicit hydrogens from an InChI string."""
    return Chem.AddHs(Chem.MolFromInchi(chi, sanitize=False, removeHs=False))


def _mol_to_inchi(mol: Mol) -> str:
    """Get a standard InChI string from an RDKit molecule."""
    return Chem.inchi.MolBlockToInchi(Chem.rdmolfiles.MolToMolBlock(mol))


def _add_atom_numbers(mol: Mol, to_number: Mapping[int, int]) -> Mol:
    """Return a copy of an RDKit molecule with ``atomLabel`` props set."""
    mol = Mol(mol)
    for atom in mol.GetAtoms():
        number = to_number[atom.GetIdx()]
        atom.SetProp("atomLabel", f"{atom.GetSymbol()}{number}")
    return mol


class Atom(Node):
    """Represents an atom in a molecule."""

    symbol: str

    def to_rdkit_atom(self) -> rdchem.Atom:
        """Convert to an RDKit Atom.

        Returns:
            The RDKit atom.
        """
        rd_atom = rdchem.Atom(self.symbol)
        rd_atom.SetNoImplicit(True)  # noqa: FBT003
        return rd_atom


class Bond(Edge):
    """Represents a bond between two atoms in a molecule."""

    def to_rdkit_bond_type(self) -> rdchem.BondType:
        """Convert to an RDKit Bond Type.

        Returns:
            The RDKit bond type.
        """
        return rdchem.BondType.SINGLE


class MolGraph(Graph[Atom, Bond]):
    """Molecular graph."""

    node_type = Atom
    edge_type = Bond


# Properties
def symbols(gra: MolGraph, keys: Sequence[int] | None = None) -> list[str]:
    """Get the atomic symbols of atoms.

    Args:
        gra: A molecular graph.
        keys: The atom keys to get symbols for; if ``None``, get all of them.

    Returns:
        The atomic symbols.
    """
    keys = node_keys(gra) if keys is None else keys
    return [gra.nodes[key][Atom.symbol] for key in keys]


def element_bonding_capacities(
    gra: MolGraph, keys: Sequence[int] | None = None
) -> list[int]:
    """Get the element bonding capacities of atoms.

    Args:
        gra: A molecular graph.
        keys: The atom keys to consider; if ``None``, consider all of them.

    Returns:
        The element bonding capacities.
    """
    return [bonding_capacity(symb) for symb in symbols(gra, keys)]


def degrees(gra: MolGraph, keys: Sequence[int] | None = None) -> list[int]:
    """Get the degrees of atoms.

    Args:
        gra: A molecular graph.
        keys: The atom keys to consider; if ``None``, consider all of them.

    Returns:
        The atom degrees.
    """
    keys = node_keys(gra) if keys is None else keys
    return [gra.degree[key] for key in keys]


def open_valences(gra: MolGraph, keys: Sequence[int] | None = None) -> list[int]:
    """Get the open valences of atoms.

    This is the bonding capacity minus the degree, i.e. the number of additional
    bonds that could be formed.

    Args:
        gra: A molecular graph.
        keys: The atom keys to consider; if ``None``, consider all of them.

    Returns:
        The open valences.
    """
    caps = element_bonding_capacities(gra, keys)
    degs = degrees(gra, keys)
    return [cap - deg for cap, deg in zip(caps, degs, strict=True)]


# Conversions from other types
def from_smiles(smi: str) -> MolGraph:
    """Instantiate a molecular graph from a SMILES string.

    Args:
        smi: A SMILES string.

    Returns:
        The molecular graph.
    """
    return from_rdkit_mol(_mol_from_smiles(smi))


def from_inchi(chi: str) -> MolGraph:
    """Instantiate a molecular graph from an InChI string.

    Args:
        chi: An InChI string.

    Returns:
        The molecular graph.
    """
    return from_rdkit_mol(_mol_from_inchi(chi))


def from_rdkit_mol(mol: Mol) -> MolGraph:
    """Instantiate a molecular graph from an RDKit molecule.

    Args:
        mol: An RDKit molecule.

    Returns:
        The molecular graph.
    """
    gra = MolGraph()

    for rd_atom in mol.GetAtoms():
        atom = Atom(symbol=rd_atom.GetSymbol())
        gra.add_node(rd_atom.GetIdx(), **atom.model_dump())

    for rd_bond in mol.GetBonds():
        bond = Bond()
        gra.add_edge(
            rd_bond.GetBeginAtomIdx(), rd_bond.GetEndAtomIdx(), **bond.model_dump()
        )
    gra.validate()
    return gra


# Conversions to other types
def inchi(gra: MolGraph) -> str:
    """Get an InChI string from a molecular graph.

    Args:
        gra: A molecular graph.

    Returns:
        The InChI string.
    """
    return _mol_to_inchi(rdkit_mol(gra))


def rdkit_mol[NodeT: Atom, EdgeT: Bond](
    gra: Graph[NodeT, EdgeT],
    *,
    label: bool = False,
) -> Mol:
    """Convert a molecular graph to an RDKit molecule.

    Args:
        gra: A molecular graph.
        label: Whether to label atoms with their graph keys.

    Returns:
        The RDKit molecule.
    """
    mol, to_key = rdkit_mol_with_index_map(gra)
    if label:
        mol = _add_atom_numbers(mol, to_key)
    return mol


def rdkit_mol_with_index_map[NodeT: Atom, EdgeT: Bond](
    gra: Graph[NodeT, EdgeT],
) -> tuple[Mol, dict[int, int]]:
    """Convert a molecular graph to an RDKit molecule with an index map.

    Args:
        gra: A molecular graph.

    Returns:
        The RDKit molecule and a mapping from RDKit atom index to graph key.
    """
    rw_mol = RWMol()
    to_idx: dict[int, int] = {}

    for key in sorted(gra.nodes()):
        atom = gra.node_type.model_validate(gra.nodes[key])
        idx = rw_mol.AddAtom(atom.to_rdkit_atom())
        to_idx[key] = idx

    for key1, key2 in gra.edges():
        bond = gra.edge_type.model_validate(gra.edges[key1, key2])
        rw_mol.AddBond(to_idx[key1], to_idx[key2], order=bond.to_rdkit_bond_type())

    to_key = dict(map(reversed, to_idx.items()))
    return rw_mol.GetMol(), to_key
