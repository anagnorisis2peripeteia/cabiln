"""Monomer CLI assembly, registration, output, and argument handling."""

import os
import pathlib
import sys
import tempfile
import unittest as _unittest

from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


class TestMonomerBuilderCLI(_unittest.TestCase):
    """Tests for cli_monomer.smiles_from_cabiln, register_monomer, and main()."""

    # ------------------------------------------------------------------ #
    # smiles_from_cabiln
    # ------------------------------------------------------------------ #

    def test_smiles_from_cabiln_single_residue(self):
        """smiles_from_cabiln('A') returns a parseable SMILES for Alanine."""
        from pyPept.interfaces.cli_monomer import smiles_from_cabiln
        smi = smiles_from_cabiln('A')
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"smiles_from_cabiln('A') gave invalid SMILES: {smi}"

    def test_smiles_from_cabiln_single_residue_correct_atoms(self):
        """Single-residue Gly assembly gives a SMILES with N and C=O but no dummy atoms."""
        from pyPept.interfaces.cli_monomer import smiles_from_cabiln
        smi = smiles_from_cabiln('G')
        assert '[*]' not in smi and '[0*]' not in smi, \
            f"Dummy atoms unexpectedly present in assembled SMILES: {smi}"
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None

    def test_smiles_from_cabiln_capped_residue(self):
        """smiles_from_cabiln('fmoc-K-am') returns parseable SMILES (capped Lys)."""
        from pyPept.interfaces.cli_monomer import smiles_from_cabiln
        smi = smiles_from_cabiln('fmoc-K-am')
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None, f"SMILES invalid: {smi}"

    def test_smiles_from_cabiln_inline_cap(self):
        """smiles_from_cabiln with inline cap notation assembles correctly."""
        from pyPept.interfaces.cli_monomer import smiles_from_cabiln
        smi = smiles_from_cabiln('fmoc-C.trt(4,2)-am')
        mol = Chem.MolFromSmiles(smi)
        assert mol is not None

    # ------------------------------------------------------------------ #
    # register_monomer (uses a temp SDF so the library is never mutated)
    # ------------------------------------------------------------------ #

    def _make_temp_sdf(self):
        """Return path to an empty temp SDF that is cleaned up after the test."""
        tmp = tempfile.NamedTemporaryFile(suffix='.sdf', delete=False, mode='w')
        tmp.close()
        self.addCleanup(lambda: pathlib.Path(tmp.name).unlink(missing_ok=True))
        return tmp.name

    def test_register_monomer_plain_smiles(self):
        """register_monomer with plain SMILES appends a record to the SDF."""
        from pyPept.interfaces.cli_monomer import register_monomer
        sdf = self._make_temp_sdf()
        # Alanine SMILES — simple and well-understood
        result = register_monomer('N[C@@H](C)C(=O)O', symbol='TestAla',
                                  sdf_path=sdf)
        assert result.chuckles, "Expected non-empty CHUCKLES"
        content = pathlib.Path(sdf).read_text(encoding='utf-8')
        assert '$$$$' in content, "SDF record terminator missing"
        assert 'TestAla' in content, "Symbol not written to SDF"

    def test_register_monomer_returns_activation_result(self):
        """register_monomer returns an ActivationResult with leaving and chem_types."""
        from pyPept.interfaces.cli_monomer import register_monomer
        from pyPept.interfaces.monomer_pipeline import ActivationResult
        sdf = self._make_temp_sdf()
        result = register_monomer('N[C@@H](CS)C(=O)O', symbol='TestCys', sdf_path=sdf)
        assert isinstance(result, ActivationResult)
        assert result.leaving
        assert result.chem_types

    def test_register_monomer_multiple_append(self):
        """Two calls to register_monomer both appear in the SDF (no overwrite)."""
        from pyPept.interfaces.cli_monomer import register_monomer
        sdf = self._make_temp_sdf()
        register_monomer('N[C@@H](C)C(=O)O', symbol='Mon1', sdf_path=sdf)
        register_monomer('NCC(=O)O', symbol='Mon2', sdf_path=sdf)
        content = pathlib.Path(sdf).read_text(encoding='utf-8')
        assert content.count('$$$$') == 2, \
            f"Expected 2 records, got {content.count('$$$$')}: {content[:400]}"
        assert 'Mon1' in content and 'Mon2' in content

    def test_register_monomer_sets_type_and_subtype(self):
        """Custom m_type and m_subtype are written to the SDF properties."""
        from pyPept.interfaces.cli_monomer import register_monomer
        sdf = self._make_temp_sdf()
        register_monomer('NCC(=O)O', symbol='GlyCap',
                         m_type='cap', m_subtype='cap', sdf_path=sdf)
        content = pathlib.Path(sdf).read_text(encoding='utf-8')
        assert 'm_type' in content
        assert 'cap' in content

    def test_register_monomer_rgroups_written(self):
        """m_Rgroups property is present in the SDF record."""
        from pyPept.interfaces.cli_monomer import register_monomer
        sdf = self._make_temp_sdf()
        register_monomer('N[C@@H](C)C(=O)O', symbol='TestAla2', sdf_path=sdf)
        content = pathlib.Path(sdf).read_text(encoding='utf-8')
        assert 'm_Rgroups' in content

    # ------------------------------------------------------------------ #
    # main() via argparse (CLI integration)
    # ------------------------------------------------------------------ #

    def _run_main(self, argv):
        """Run cli_monomer.main() with the given argv list; return (stdout, stderr, rc)."""
        from pyPept.interfaces import cli_monomer
        import io
        from unittest.mock import patch
        out_buf = io.StringIO()
        err_buf = io.StringIO()
        rc = 0
        with patch('sys.argv', ['pyPept-monomer-add'] + argv), \
             patch('sys.stdout', out_buf), \
             patch('sys.stderr', err_buf):
            try:
                cli_monomer.main()
            except SystemExit as exc:
                rc = int(exc.code) if exc.code is not None else 0
        return out_buf.getvalue(), err_buf.getvalue(), rc

    def test_main_smiles_path(self):
        """CLI --smiles path exits 0 and prints 'Registered'."""
        sdf = self._make_temp_sdf()
        out, err, rc = self._run_main([
            '--smiles', 'N[C@@H](C)C(=O)O',
            '--symbol', 'CLIAla',
            '--sdf', sdf,
        ])
        assert rc == 0, f"Non-zero exit code; stderr: {err}"
        assert 'Registered' in out, f"Expected 'Registered' in stdout: {out}"
        assert 'CLIAla' in out

    def test_main_from_cabiln_path(self):
        """CLI --from-cabiln path assembles and registers correctly."""
        sdf = self._make_temp_sdf()
        out, err, rc = self._run_main([
            '--from-cabiln', 'G',
            '--symbol', 'CLIGly',
            '--sdf', sdf,
        ])
        assert rc == 0, f"Non-zero exit code; stderr: {err}"
        assert 'Assembled SMILES' in out, f"Assembly message missing: {out}"
        assert 'Registered' in out

    def test_main_sdf_is_written(self):
        """CLI writes a valid SDF record to the given --sdf path."""
        sdf = self._make_temp_sdf()
        self._run_main([
            '--smiles', 'NCC(=O)O',
            '--symbol', 'CLIGly2',
            '--sdf', sdf,
        ])
        content = pathlib.Path(sdf).read_text(encoding='utf-8')
        assert '$$$$' in content
        assert 'CLIGly2' in content

    def test_main_bad_smiles_exits_nonzero(self):
        """CLI exits non-zero when --smiles cannot be pre-activated."""
        sdf = self._make_temp_sdf()
        out, err, rc = self._run_main([
            '--smiles', 'NOTVALIDSMILES???',
            '--symbol', 'Bad',
            '--sdf', sdf,
        ])
        assert rc != 0, f"Expected non-zero exit for invalid SMILES; got rc={rc}"

    def test_main_missing_symbol_exits_nonzero(self):
        """CLI exits non-zero when --symbol is omitted."""
        sdf = self._make_temp_sdf()
        out, err, rc = self._run_main([
            '--smiles', 'NCC(=O)O',
            '--sdf', sdf,
        ])
        assert rc != 0

    def test_main_no_source_exits_nonzero(self):
        """CLI exits non-zero when neither --smiles nor --from-cabiln is given."""
        sdf = self._make_temp_sdf()
        out, err, rc = self._run_main([
            '--symbol', 'NoSrc',
            '--sdf', sdf,
        ])
        assert rc != 0
