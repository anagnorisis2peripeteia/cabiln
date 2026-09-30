"""Small generators shared by the chemistry property tests."""

from contextlib import contextmanager
from pathlib import Path
import random
from tempfile import TemporaryDirectory

import pytest
from rdkit import Chem

from pyPept import monomer_store


def reordered(molecule, seed):
    order = list(range(molecule.GetNumAtoms()))
    random.Random(seed).shuffle(order)
    result = Chem.RenumberAtoms(molecule, order)
    # RDKit renumbering drops molecule-level SDF fields; keep the definition.
    for name in molecule.GetPropNames(includePrivate=True, includeComputed=False):
        result.SetProp(name, molecule.GetProp(name))
    assert Chem.MolToSmiles(result) == Chem.MolToSmiles(molecule)
    return result


def alternate_smiles(source, seed):
    molecule = Chem.MolFromSmiles(source)
    assert molecule is not None, source
    return Chem.MolToSmiles(reordered(molecule, seed), canonical=False)


def write_library(path, molecules):
    with Chem.SDWriter(str(path)) as writer:
        for molecule in molecules:
            writer.write(molecule)


def _clear_library_caches():
    from pyPept import recognition
    from pyPept.web.cache import clear_render_cache

    monomer_store._invalidate_sdf()
    monomer_store._binding_cache.clear()
    recognition._pattern_cache.clear()
    recognition._core_cache.clear()
    clear_render_cache()


@contextmanager
def isolated_library(symbols):
    """A fresh library and environment for each Hypothesis example, not fixture."""
    bundled = Path(monomer_store.__file__).with_name("data") / "monomers.sdf"
    selected = {
        molecule.GetProp("m_abbr"): molecule
        for molecule in Chem.SDMolSupplier(str(bundled), removeHs=False)
        if molecule is not None and molecule.GetProp("m_abbr") in symbols
    }
    assert set(selected) == set(symbols)
    with TemporaryDirectory(prefix="cabiln-chemistry-fuzz-") as temporary:
        path = Path(temporary) / "monomers.sdf"
        write_library(path, (selected[symbol] for symbol in symbols))
        with pytest.MonkeyPatch.context() as patch:
            patch.setenv("CABILN_MONOMER_LIBRARY", str(path))
            _clear_library_caches()
            try:
                yield path
            finally:
                _clear_library_caches()
        _clear_library_caches()
