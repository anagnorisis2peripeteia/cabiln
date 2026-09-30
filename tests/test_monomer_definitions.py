"""Shared definitions remain current and detached from individual sequences."""

import os

import pytest
from fastapi.testclient import TestClient
from rdkit import Chem

from pyPept import monomer_store, recognition
from pyPept.molecule import Molecule
from pyPept.sequence import Sequence, get_monomer_info
from pyPept.web.app import create_app


@pytest.fixture
def definitions(tmp_path, monkeypatch):
    _, records = monomer_store._load_sdf()
    path = tmp_path / "monomers.sdf"
    with Chem.SDWriter(str(path)) as writer:
        for name in ("A", "G"):
            writer.write(records[name])
    monkeypatch.setenv("CABILN_MONOMER_LIBRARY", str(path))
    monomer_store._invalidate_sdf()
    yield path
    monomer_store._invalidate_sdf()


def canonical(source):
    return Chem.MolToSmiles(Molecule(Sequence(source), depiction=None).mol)


def test_sequences_cannot_mutate_shared_definitions(definitions):
    expected = canonical("A-G")
    sequence = Sequence("A-G")
    sequence.monomer_df.at["A", "m_romol"].GetAtomWithIdx(1).SetFormalCharge(1)
    sequence.monomer_df.at["A", "m_Rgroups"][0] = "[F]"
    sequence.monomer_df.attrs["_synonyms"]["CorruptAlias"] = "A"
    sequence.monomer_df.drop(index="G", inplace=True)

    assert canonical("A-G") == expected
    with pytest.raises(ValueError, match="not in the monomer library"):
        Sequence("CorruptAlias")
    synthetic = Sequence("G-<N[C@@H](CCC(F)(F)F)C(=O)O>-A")
    assert any(name.startswith("__SYN_") for name in synthetic.monomer_df.index)
    assert set(get_monomer_info(str(definitions)).index) == {"A", "G"}


def test_resolved_assembly_keeps_bound_templates_and_occurrence_ownership(definitions):
    from pyPept.peptide import Peptide

    sequence = Sequence("A-G")
    nitrogen = next(
        atom for atom in sequence.s_monomers[0]["m_romol"].GetAtoms()
        if atom.GetAtomicNum() == 7
    )
    nitrogen.SetMonomerInfo(
        Chem.AtomPDBResidueInfo(atomName=" N  ", residueName="ALA", residueNumber=7)
    )
    peptide = Peptide.from_sequence(sequence, (17, 4))
    peptide = Peptide(tuple(reversed(peptide.occurrences)), peptide.connections)
    glycine = sequence.s_monomers[1]
    glycine["m_Rgroups"][1] = "[18O-]"
    carbonyl = next(
        atom.GetNeighbors()[0] for atom in glycine["m_romol"].GetAtoms()
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == 2
    )
    carbonyl.SetIsotope(13)

    assembled = Molecule(peptide, depiction=None)
    expected = Chem.MolFromSmiles("N[C@@H](C)C(=O)NCC(=O)O")
    assert Chem.MolToSmiles(assembled.mol) == Chem.MolToSmiles(expected)
    changed = Chem.MolFromSmiles("N[C@@H](C)C(=O)NC[13C](=O)[18O-]")
    assert Chem.MolToSmiles(Molecule(sequence, depiction=None).mol) == (
        Chem.MolToSmiles(changed)
    )
    owners = assembled.get_residue_atom_map()
    assert {identity: len(atoms) for identity, atoms in owners.items()} == {17: 5, 4: 5}
    assert sorted(atom for atoms in owners.values() for atom in atoms) == list(
        range(assembled.mol.GetNumAtoms())
    )
    pdb_atoms = [
        assembled.mol.GetAtomWithIdx(index).GetPDBResidueInfo()
        for index in owners[17]
    ]
    assert [
        (info.GetName(), info.GetResidueNumber()) for info in pdb_atoms if info
    ] == [(" N  ", 7)]


def test_synonym_changes_invalidate_parser_and_rendered_results(definitions):
    aliases = definitions.with_name("monomers.csv")
    aliases.write_text("token,synonyms\nA,CustomAlias\n")
    with TestClient(create_app()) as client:
        first = client.post("/render", json={"cabiln": "CustomAlias"})
        assert first.status_code == 200, first.text
        first_mol = Chem.MolFromMolBlock(first.json()["mol_block"])
        assert Chem.MolToSmiles(first_mol) == canonical("A")

        replacement = aliases.with_suffix(".new")
        replacement.write_text("token,synonyms\nG,CustomAlias\n")
        replacement.replace(aliases)
        second = client.post("/render", json={"cabiln": "CustomAlias"})
        assert second.status_code == 200, second.text
        second_mol = Chem.MolFromMolBlock(second.json()["mol_block"])
        assert Chem.MolToSmiles(second_mol) == canonical("G")

        aliases.unlink()
        assert client.post("/render", json={"cabiln": "CustomAlias"}).status_code == 400


def replace_preserving_metadata(path, before, after):
    """Emulate an external atomic update that preserves size and mtime."""
    metadata = path.stat()
    content = path.read_bytes()
    changed = content.replace(before, after)
    assert changed != content and len(changed) == len(content)
    replacement = path.with_suffix(".new")
    replacement.write_bytes(changed)
    os.utime(replacement, ns=(metadata.st_atime_ns, metadata.st_mtime_ns))
    replacement.replace(path)
    assert path.stat().st_size == metadata.st_size
    assert path.stat().st_mtime_ns == metadata.st_mtime_ns


def test_palette_version_header_tracks_unchanged_metadata_and_alias_updates(definitions):
    with TestClient(create_app()) as client:
        initial = client.get("/monomers")
        initial_version = initial.headers.get("x-library-version")
        assert initial_version
        assert isinstance(initial.json(), list)
        cached = client.get("/monomers")
        assert cached.headers["x-library-version"] == initial_version
        assert cached.json() == initial.json()

        metadata = definitions.stat()
        records = list(Chem.SDMolSupplier(str(definitions), removeHs=False))
        alanine = records[0]
        alpha = next(
            atom for atom in alanine.GetAtoms()
            if atom.GetAtomicNum() == 6
            and any(neighbor.GetAtomicNum() == 7 for neighbor in atom.GetNeighbors())
        )
        alpha.SetAtomicNum(14)
        Chem.SanitizeMol(alanine)
        replacement = definitions.with_suffix(".new")
        with Chem.SDWriter(str(replacement)) as writer:
            for record in records:
                writer.write(record)
        os.utime(replacement, ns=(metadata.st_atime_ns, metadata.st_mtime_ns))
        replacement.replace(definitions)
        assert definitions.stat().st_size == metadata.st_size
        assert definitions.stat().st_mtime_ns == metadata.st_mtime_ns

        changed = client.get("/monomers")
        assert changed.status_code == 200
        # Labels/sites stay the same; chemistry evidence must reflect the change.
        def without_quality(rows):
            return [
                {key: value for key, value in row.items() if key != "quality"}
                for row in rows
            ]
        assert without_quality(changed.json()) == without_quality(initial.json())
        old_a = next(row for row in initial.json() if row["abbr"] == "A")
        new_a = next(row for row in changed.json() if row["abbr"] == "A")
        assert new_a["quality"]["status"] == "unreviewed"
        assert new_a["quality"]["definition_hash"] != old_a["quality"]["definition_hash"]
        assert changed.headers["x-library-version"] != initial_version

        aliases = definitions.with_name("monomers.csv")
        aliases.write_text("token,synonyms\nA,CustomAlias\n")
        renamed = client.get("/monomers")
        assert renamed.json() == changed.json()
        assert renamed.headers["x-library-version"] != changed.headers["x-library-version"]


def test_palette_does_not_certify_a_revision_that_changes_during_loading(
    definitions, monkeypatch
):
    from pyPept.web import cache, monomers

    attachment_sites = monomers.attachment_sites
    changed = False

    def inspect_and_change_library(*args):
        nonlocal changed
        if not changed:
            changed = True
            definitions.with_name("monomers.csv").write_text(
                "token,synonyms\nA,CustomAlias\n"
            )
        return attachment_sites(*args)

    with TestClient(create_app()) as client:
        # Trigger a request-time cache miss after startup has prepared the palette.
        cache.clear_render_cache()
        monkeypatch.setattr(monomers, "attachment_sites", inspect_and_change_library)
        interrupted = client.get("/monomers")
        assert interrupted.status_code == 200
        assert "x-library-version" not in interrupted.headers
        stable = client.get("/monomers")
        assert stable.status_code == 200
        assert stable.headers["x-library-version"]
        assert stable.json() == interrupted.json()


def test_preserved_mtime_replacement_refreshes_parser_and_recognition(definitions):
    source = Chem.MolFromSmiles("N[C@@H](C)C(=O)NCC(=O)O")
    expected = Chem.MolToSmiles(source)
    assert canonical("A-G") == expected
    assert {
        candidate.symbol for candidate in recognition.recognize(source).covers[0]
    } == {"A", "G"}
    previous = monomer_store.library_version(include_aliases=False)

    replace_preserving_metadata(definitions, b"\nA\n", b"\nB\n")

    assert monomer_store.library_version(include_aliases=False) != previous
    assert set(monomer_store._load_sdf()[1]) == {"B", "G"}
    assert set(monomer_store.monomer_table().index) == {"B", "G"}
    with pytest.raises(ValueError, match="not in the monomer library"):
        Sequence("A-G")
    assert canonical("B-G") == expected
    assert {
        candidate.symbol for candidate in recognition.recognize(source).covers[0]
    } == {"B", "G"}


def test_alias_updates_refresh_definitions_without_recompiling_structures(
    definitions, monkeypatch
):
    aliases = definitions.with_name("monomers.csv")
    aliases.write_text("token,synonyms\nA,CustomAlias\n")
    source = Chem.MolFromSmiles("N[C@@H](C)C(=O)NCC(=O)O")
    expected = recognition.recognize(source)
    assert canonical("CustomAlias") == canonical("A")
    structure_version = monomer_store.library_version(include_aliases=False)
    complete_version = monomer_store.library_version()

    def unexpected_compile(*_args, **_kwargs):
        pytest.fail("CSV aliases must not invalidate SDF recognition patterns")

    monkeypatch.setattr(recognition, "_compile_patterns", unexpected_compile)
    replace_preserving_metadata(aliases, b"A,CustomAlias", b"G,CustomAlias")

    assert monomer_store.library_version() != complete_version
    assert monomer_store.library_version(include_aliases=False) == structure_version
    assert canonical("CustomAlias") == canonical("G")
    assert recognition.recognize(source) == expected


@pytest.mark.parametrize("loader", ["raw", "table"])
def test_replacement_during_load_is_retried_with_matching_snapshot(
    definitions, monkeypatch, loader
):
    calls = []
    if loader == "raw":
        read = Chem.SDMolSupplier

        def interrupted_read(*args, **kwargs):
            result = list(read(*args, **kwargs))
            calls.append(None)
            if len(calls) == 1:
                replace_preserving_metadata(definitions, b"\nA\n", b"\nB\n")
            return result

        monkeypatch.setattr(Chem, "SDMolSupplier", interrupted_read)
        stamp, _, records = monomer_store._load_sdf_snapshot()
        assert stamp == monomer_store.library_version(include_aliases=False)
        assert set(records) == {"B", "G"}
    else:
        read = monomer_store._read_monomer_table

        def interrupted_read(path):
            result = read(path)
            calls.append(None)
            if len(calls) == 1:
                replace_preserving_metadata(definitions, b"\nA\n", b"\nB\n")
            return result

        monkeypatch.setattr(monomer_store, "_read_monomer_table", interrupted_read)
        assert set(monomer_store.monomer_table().index) == {"B", "G"}
    assert len(calls) == 2
    assert set(monomer_store._load_sdf()[1]) == {"B", "G"}
    assert set(monomer_store.monomer_table().index) == {"B", "G"}


@pytest.mark.parametrize("loader", ["raw", "table"])
def test_continually_replaced_library_fails_instead_of_caching_mixed_data(
    definitions, monkeypatch, loader
):
    calls = []

    def replace_after_read():
        calls.append(None)
        old, new = (b"A", b"B") if len(calls) % 2 else (b"B", b"A")
        replace_preserving_metadata(
            definitions, b"\n" + old + b"\n", b"\n" + new + b"\n"
        )

    if loader == "raw":
        read = Chem.SDMolSupplier

        def interrupted_read(*args, **kwargs):
            result = list(read(*args, **kwargs))
            replace_after_read()
            return result

        monkeypatch.setattr(Chem, "SDMolSupplier", interrupted_read)
        load = monomer_store._load_sdf
    else:
        read = monomer_store._read_monomer_table

        def interrupted_read(path):
            result = read(path)
            replace_after_read()
            return result

        monkeypatch.setattr(monomer_store, "_read_monomer_table", interrupted_read)
        load = monomer_store.monomer_table

    with pytest.raises(ValueError, match="changed repeatedly"):
        load()
    assert len(calls) == 3


def test_raw_structure_snapshot_preserves_explicit_hydrogens(definitions):
    molecules = [Chem.AddHs(mol) for mol in monomer_store._load_sdf()[0]]
    with Chem.SDWriter(str(definitions)) as writer:
        for mol in molecules:
            writer.write(mol)
    stamp, loaded, _ = monomer_store._load_sdf_snapshot()
    assert stamp == monomer_store.library_version(include_aliases=False)
    assert all(
        any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()) for mol in loaded
    )
    assert [mol.GetNumAtoms() for mol in loaded] == [
        mol.GetNumAtoms() for mol in molecules
    ]
