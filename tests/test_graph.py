"""Graph tests."""

import enum

import pytest

from autoccv import graph
from autoccv.graph.base import FieldEnumMismatchError, FieldEnumModel


def test__smiles() -> None:
    """Test graph smiles."""
    water_smiles = "O"
    water_inchi = "InChI=1S/H2O/h1H2"
    water = graph.from_smiles(water_smiles)
    assert graph.inchi(water) == water_inchi


def test__inchi() -> None:
    """Test graph inchi."""
    water_inchi = "InChI=1S/H2O/h1H2"
    water = graph.from_inchi(water_inchi)
    assert graph.inchi(water) == water_inchi


def test__remove_edges() -> None:
    """Test graph remove edges."""
    water_smiles = "O"
    oh_h_smiles = "[OH].[H]"
    water = graph.from_smiles(water_smiles)
    oh_h_ref = graph.from_smiles(oh_h_smiles)
    oh_h = graph.remove_edges(water, [(0, 1)])
    assert graph.is_isomorphic(oh_h, oh_h_ref)


def test__symbols() -> None:
    """Test graph symbols."""
    water_smiles = "O"
    water = graph.from_smiles(water_smiles)
    assert graph.symbols(water) == ["O", "H", "H"]


def test__field_enum_matches_fields() -> None:
    """A nested Field enum in parity with the model fields is accepted."""

    class Good(FieldEnumModel):
        class Field(enum.StrEnum):
            a = "a"
            b = "b"

        a: int
        b: str

    assert Good.Field.a == "a"
    assert set(Good.Field) == set(Good.model_fields)


def test__field_enum_missing_member() -> None:
    """A field with no Field member fails at class-creation time."""
    with pytest.raises(FieldEnumMismatchError, match="missing from Field"):

        class MissingMember(FieldEnumModel):
            class Field(enum.StrEnum):
                a = "a"

            a: int
            b: str


def test__field_enum_extra_member() -> None:
    """A Field member that is not a field fails at class-creation time."""
    with pytest.raises(FieldEnumMismatchError, match="not a model field"):

        class ExtraMember(FieldEnumModel):
            class Field(enum.StrEnum):
                a = "a"
                c = "c"

            a: int


@pytest.mark.parametrize(
    ("smi", "ref"),
    [
        ("[H]", [1]),
        ("[He]", [0]),
        ("O", [0, 0, 0]),
        ("[OH]", [1, 0]),
        ("C=C", [1, 1, 0, 0, 0, 0]),
        ("C#C", [2, 2, 0, 0]),
        ("O[O]", [0, 1, 0]),
    ],
)
def test__open_valences(smi: str, ref: list[int]) -> None:
    """Test graph open valences."""
    gra = graph.from_smiles(smi)
    assert graph.open_valences(gra) == ref
