"""Resolve peptide notation into monomer occurrences and numbered connections.

From pyPept: a python library to generate atomistic 2D and 3D representations
of peptides (Journal of Cheminformatics, 2023).
"""

__credits__ = ["Rodrigo Ochoa", "J.B. Brown", "Thomas Fox"]
__license__ = "MIT"

import copy
import re
import warnings
from dataclasses import dataclass

import numpy as np
from rdkit import Chem

from pyPept.source import SourceText, join as _source_join, sub as _source_sub
from pyPept.notation import MAX_NOTATION_CHARACTERS, split_outside
# Compatibility imports preserve the existing library entry points.
from pyPept.notation_conversion import (
    _OLD_BILN_RE, biln_to_cabiln, cabiln_to_branch, cabiln_to_bracket,
    colorize_cabiln,
)
from pyPept.notation_lowering import _expand_inline_caps, _preprocess_cabiln
from pyPept.pdb_names import correct_pdb_atoms, greekify, get_atom_by_name, get_monomer_codes
from pyPept.monomer_store import library_resource
from pyPept.synthetic import expand_synthetic_tokens, add_synthetic_monomers


def _attachment_idx(mol, slot):
    """Return attachment atom index for R-group slot (1-based), or None.

    Locates the dummy atom whose isotope equals slot (CHUCKLES convention:
    slot 1 → [1*], slot 2 → [2*], …) and returns its neighbour's index —
    i.e. the heavy atom where the inter-monomer bond will form.
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == slot:
            nb = atom.GetNeighbors()
            return nb[0].GetIdx() if nb else None
    return None


def _rgroup_atom_idx(mol, slot):
    """Return the dummy atom index for R-group slot (1-based), or None."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == slot:
            return atom.GetIdx()
    return None


def _slot_for_attachment(mol, atom_idx):
    """Return R-group slot (1-based) whose dummy neighbours atom_idx, or None."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            for nb in atom.GetNeighbors():
                if nb.GetIdx() == atom_idx:
                    return atom.GetIsotope()
    return None


def _emit_warning(message, category=UserWarning, stacklevel=1, *, warning_sink=None):
    """Keep request diagnostics local while retaining the library warning API."""
    if warning_sink is None:
        warnings.warn(message, category, stacklevel=stacklevel + 1)
    else:
        warning_sink(str(message))


def _bond_chemistry_diagnostic(mol1, at1, mol2, at2, bond_label='', warning_sink=None):
    """Emit established exotic-bond warnings and return a legacy rejection.

    Atom pairs describe these diagnostics, not numbered-site eligibility. The
    latter is determined by effective attachment types and the reaction index.
    """
    def _is_carbonyl_carbon(mol, atom):
        return (atom.GetAtomicNum() == 6
                and any(nb.GetAtomicNum() == 8
                        and mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                               .GetBondTypeAsDouble() == 2.0
                        for nb in atom.GetNeighbors()))

    a1 = mol1.GetAtomWithIdx(at1)
    a2 = mol2.GetAtomWithIdx(at2)
    sym1, sym2 = a1.GetAtomicNum(), a2.GetAtomicNum()
    pair = frozenset([sym1, sym2])

    # N(7)–C(6): amide/isopeptide — C must be carbonyl
    if pair == frozenset([7, 6]):
        c_mol, c_atom = (mol1, a1) if sym1 == 6 else (mol2, a2)
        if _is_carbonyl_carbon(c_mol, c_atom):
            return
        _emit_warning(
            f"Bond {bond_label}: N–C join where C is not a carbonyl carbon. "
            "This forms a C–N bond without amide character (e.g. reductive amination "
            "product). Intentional?",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # S(16)–S(16): disulfide
    if pair == frozenset([16]):
        return

    # O(8)–C(6): ester — C must be carbonyl
    if pair == frozenset([8, 6]):
        c_mol, c_atom = (mol1, a1) if sym1 == 6 else (mol2, a2)
        if _is_carbonyl_carbon(c_mol, c_atom):
            return
        _emit_warning(
            f"Bond {bond_label}: O–C join where C is not a carbonyl carbon. "
            "This forms an ether, not an ester. Intentional?",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # Se(34)–Se(34): diselenide
    if pair == frozenset([34]):
        return

    # N(7)–N(7): hydrazide / hydrazone
    if pair == frozenset([7]):
        _emit_warning(
            f"Bond {bond_label}: N–N join (hydrazide/hydrazone). "
            "Unusual in peptide chemistry — check R-group assignments.",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # S–C(=O): thioester; S–C(aliphatic): thioether
    if pair == frozenset([16, 6]):
        c_mol, c_atom = (mol1, a1) if sym1 == 6 else (mol2, a2)
        if _is_carbonyl_carbon(c_mol, c_atom):
            _emit_warning(
                f"Bond {bond_label}: S–C(=O) thioester bond. "
                "Valid for native chemical ligation but unusual for standard assembly. "
                "Intentional?",
                UserWarning, stacklevel=4, warning_sink=warning_sink,
            )
        else:
            _emit_warning(
                f"Bond {bond_label}: S–C(aliphatic) thioether bond. "
                "Valid for thioether-linked protecting groups (e.g. trt, acm) "
                "and related ligation chemistry. Intentional?",
                UserWarning, stacklevel=4, warning_sink=warning_sink,
            )
        return

    # S–N: sulfenamide
    if pair == frozenset([16, 7]):
        _emit_warning(
            f"Bond {bond_label}: S–N sulfenamide bond. "
            "Valid for Cys PTMs, oxidative macrolactamisation, bioconjugation, "
            "and NCL variant chemistry, but unusual in standard assembly. "
            "Intentional?",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # C(6)–C(6): RCM olefin staple (alkene–alkene metathesis) or other C–C bond
    if pair == frozenset([6]):
        def _is_vinyl(mol, atom):
            return any(
                mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2.0
                and nb.GetAtomicNum() == 6
                for nb in atom.GetNeighbors()
            )
        if _is_vinyl(mol1, a1) and _is_vinyl(mol2, a2):
            return  # all-hydrocarbon alkene staple (RCM)
        _emit_warning(
            f"Bond {bond_label}: C–C inter-monomer bond (non-vinyl). "
            "Unusual — check R-group assignments unless this is intentional "
            "bioconjugation chemistry.",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # O(8)–N(7): hydroxamic acid / hydroxylamine / isoxazole linkage
    if pair == frozenset([8, 7]):
        _emit_warning(
            f"Bond {bond_label}: O–N hydroxylamine/hydroxamic acid bond. "
            "Valid for O-amino acid caps (OBn_, OMe_) and hydroxamate "
            "bioconjugation, but unusual in standard peptide assembly. "
            "Intentional?",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # Retained for callers of the bare-molecule compatibility helper. A
    # registered numbered-site reaction can support additional element pairs.
    sym_names = {6: 'C', 7: 'N', 8: 'O', 16: 'S', 34: 'Se'}
    s1 = sym_names.get(sym1, str(sym1))
    s2 = sym_names.get(sym2, str(sym2))
    return (
        f"Bond {bond_label}: {s1}–{s2} inter-monomer bond has no recognised "
        "peptide chemistry context. Check that the correct R-group numbers "
        "were specified for both monomers."
    )


def _check_bond_chemistry(mol1, at1, mol2, at2, bond_label='', warning_sink=None):
    """Compatibility diagnostics for bare molecules without numbered sites."""
    diagnostic = _bond_chemistry_diagnostic(
        mol1, at1, mol2, at2, bond_label, warning_sink)
    if diagnostic:
        raise ValueError(diagnostic)


def _check_parsed_connection(mol1, at1, slot1, mol2, at2, slot2, *,
                                leaving_groups1, leaving_groups2,
                                bond_label='', warning_sink=None):
    """Use registered reactions before legacy atom-pair diagnostics.

    Sequence historically also accepts some graphs without an assembly reaction.
    The legacy diagnostic preserves that parse-only contract; reaction support
    belongs to resolve_connection, as used by assembly and the builder.
    """
    from pyPept.attachments import resolve_connection

    connection = resolve_connection(
        mol1, slot1, mol2, slot2,
        leaving_groups1=leaving_groups1, leaving_groups2=leaving_groups2)
    if connection.reaction is not None:
        return
    diagnostic = _bond_chemistry_diagnostic(
        mol1, at1, mol2, at2, bond_label, warning_sink)
    if diagnostic:
        raise ValueError(diagnostic)


class SequenceConstants:
    """
    A class to hold defaults values related to pyPept.Sequence objects.
    """
    def_path = "pyPept.data"
    def_lib_filename = "monomers.sdf"
    monomer_join = "-"
    chain_separator = "."
    csv_separator = ","
    helm_polymer = '|'
    max_rgroups = 4


@dataclass
class ValidationReport:
    """Result of Sequence.validate(). Check .ok before calling .build()."""
    biln: str
    ok: bool
    bonds: list    # (m1_idx, slot1, m2_idx, slot2) tuples
    warnings: list # chemistry warning strings
    errors: list   # parse/validation error strings

    def __post_init__(self):
        self._seq = None  # set by Sequence.validate()

    def build(self):
        """Return the constructed Sequence, or raise ValueError if validation failed."""
        if not self.ok:
            msg = self.errors[0] if self.errors else "validation failed"
            raise ValueError(f"Cannot build invalid sequence: {msg}")
        return self._seq


class Sequence:
    """
    This class holds the information for a given sequence in BILN format,
    parsing monomer and bond information.
    """

    def __init__(self, input_biln, path=SequenceConstants.def_path,
                 monomer_lib=SequenceConstants.def_lib_filename, fmt=None, *, warning_sink=None, track_source=False):
        """
        Instantiate a pyPept.Sequence object with the input BILN/CABILN sequence.

        :param input_biln: CABILN sequence string (default), or old BILN when fmt='biln'.
        :type input_biln: str
        :param path: an override option to specify a new location of a monomer
                     library.
        :type path: str
        :param monomer_lib: an override option to specific a different monomer
                           library file name.
        :type monomer_lib: str
        :param fmt: pass ``'biln'`` to accept old BILN crosslink notation
            (``Token(bondID,Rgroup)``), which is auto-converted to CABILN.
            Without this flag, old BILN notation raises ``ValueError``.
        :type fmt: str or None
        :param track_source: retain original token locations in ``s_sources``.
        """
        if not isinstance(input_biln, str) or not input_biln.strip():
            raise ValueError("CABILN must be a non-empty string.")
        if len(input_biln) > MAX_NOTATION_CHARACTERS:
            raise ValueError(
                f'Notation exceeds {MAX_NOTATION_CHARACTERS:,} characters; '
                'split it into smaller documents.')
        if fmt not in (None, 'cabiln', 'biln'):
            raise ValueError(f"Unknown sequence format: {fmt!r}")

        self.s_inputbiln = input_biln
        self._warning_sink = warning_sink
        self.s_mid = -1
        self.s_bonds = []
        self.s_nbonds = 0
        self._used_slots = set()
        self.s_monomers = []
        self.s_nmonomers = 0
        self.__is_valid = True

        if track_source and fmt == 'biln':
            raise ValueError('Source tracking requires CABILN; convert old BILN first.')

        # Old BILN bare-integer crosslink notation requires explicit opt-in via fmt='biln'.
        m = _OLD_BILN_RE.search(input_biln)
        if m:
            tok, bid, rg = m.group(1), m.group(2), m.group(3)
            if fmt != 'biln':
                raise ValueError(
                    f"Old BILN crosslink notation detected: {tok}({bid},{rg}). "
                    f"Pass fmt='biln' to auto-convert, or call biln_to_cabiln() explicitly."
                )
            input_biln = biln_to_cabiln(input_biln)

        # Protect synthetic SMILES before interpreting CABILN separators.
        expanded = SourceText.original(input_biln) if track_source else input_biln
        self.s_sources = []

        expanded, self._synthetic_aa_smiles, self._synthetic_caps = expand_synthetic_tokens(expanded)
        expanded, _branch_rgroup = _expand_inline_caps(expanded)

        seq = split_outside(expanded,
                            by_element=SequenceConstants.monomer_join,
                            outside='[]')
        if any(not token.strip() for token in seq):
            raise ValueError("Empty monomer between backbone separators.")
        seq = _source_join(SequenceConstants.monomer_join, seq)
        self.s_biln = _source_sub(r'[-]*\.[-]*',
                             SequenceConstants.chain_separator, seq)

        monomer_df_filepath = library_resource(path, monomer_lib)
        self.monomer_df = get_monomer_info(str(monomer_df_filepath))

        add_synthetic_monomers(self.monomer_df, self._synthetic_aa_smiles, self._synthetic_caps)

        self.__parse_biln()
        self.__parse_biln_bonds()

    def __parse_biln(self):
        """Split BILN string into individual monomers, do some basic checks,
        then construct the monomers and store the chemical representations
        of the monomers in the Sequence object.
        """

        biln = self.s_biln
        chains = split_outside(biln,
                               by_element=SequenceConstants.chain_separator,
                               outside='[]')

        self.s_chains = {}
        self.s_chains["s_nChains"] = len(chains)
        self.s_chains["s_cType"] = [None] * len(chains)
        self.s_chains["s_monomerIDs"] = [None] * len(chains)

        m_idx = 0
        for num_chain, chain in enumerate(chains):
            residues = split_outside(chain,
                                     by_element=SequenceConstants.monomer_join,
                                     outside='[]')
            monomer_types = set()
            monomer_ids = []
            for _res_pos, res in enumerate(residues):
                # Strip bond annotations: (n,m) and (!n,m) crosslink IDs
                resname = _source_sub(r'\((?:!\w+|\d+),\d+\)', '', res)
                if isinstance(resname, SourceText):
                    self.s_sources.append(resname.tracker.occurrence(resname))
                resname = re.sub(r'^[\[{](.*?)[\]}]$', '\\1', resname)


                if resname not in self.monomer_df.index:
                    if resname in self.monomer_df.attrs.get('_ambiguous_aliases', {}):
                        raise ValueError(
                            f"Ambiguous monomer alias '{resname}' names multiple library entries; "
                            "use a canonical symbol instead.")
                    # 1) Check stored abbreviations and unambiguous CSV aliases.
                    synonyms = self.monomer_df.attrs.get('_synonyms', {})
                    if resname in synonyms:
                        resname = synonyms[resname]
                    # 2) Check degenerate cap aliases (R1/R2 structural pairs)
                    elif resname in self.monomer_df.attrs.get('_degen_aliases', {}):
                        degen = self.monomer_df.attrs['_degen_aliases']
                        variants = degen[resname]
                        res_idx = _res_pos
                        if res_idx == 0:
                            chosen = [v for v in variants
                                      if v.endswith('_') and not v.startswith('_')]
                            if not chosen:
                                chosen = [v for v in variants if v.endswith('_')]
                        elif res_idx == len(residues) - 1:
                            chosen = [v for v in variants
                                      if v.startswith('_') and not v.endswith('_')]
                            if not chosen:
                                chosen = [v for v in variants if v.startswith('_')]
                        else:
                            raise ValueError(
                                f"Degenerate cap '{resname}' used mid-chain "
                                f"(position {res_idx}); cannot resolve variant.")
                        if chosen:
                            resname = chosen[0]
                        else:
                            raise ValueError(
                                f"Monomer {resname!r} has no terminal variant for this position.")
                    else:
                        raise ValueError(
                            f"Monomer {resname!r} is not in the monomer library.")

                mm_info = self.monomer_df.loc[resname, :].to_dict()
                mm_value = {'m_name': resname,
                            'm_name_in_biln': res,
                            'm_chainID': num_chain}
                mm_value.update(mm_info)

                self.__append_monomer(mm_value)

                monomer_ids.append(m_idx)
                monomer_types.add(mm_info['m_type'])
                m_idx += 1

            if len(monomer_types) == 1:
                if monomer_types == {'aa'}:
                    self.s_chains["s_cType"][num_chain] = 'peptide'
                else:
                    self.s_chains["s_cType"][num_chain] = 'chem'
            else:
                self.s_chains["s_cType"][num_chain] = 'mixed'

            self.s_chains["s_monomerIDs"][num_chain] = monomer_ids


    def __append_monomer(self, monomer):
        """Place a new monomer at the end of the sequence, i.e.
        append the current monomer to the monomer list.

        :param monomer: monomer dictionary of added monomer
        :type monomer: dictionary
        """
        self.s_mid += 1

        keys_needed = ['m_name', 'm_abbr', 'm_name_in_biln', 'm_type',
                       'm_subtype', 'm_chainID', 'm_Rgroups', 'm_romol',
                       'm_chem_types']

        for key in keys_needed:
            try:
                monomer[key]
            except KeyError:
                raise ValueError(f'Key {key} missing in monomer description')

        m_dict = {}
        m_dict.update({key: monomer[key] for key in keys_needed})

        self.s_nmonomers += 1
        self.s_monomers.append(m_dict)

    def __parse_biln_bonds(self):
        """
        Function to parse bond information from the BILN string
        """

        bond_info = []
        bond_info_helm = []

        residues = split_outside(self.s_biln,
                                 SequenceConstants.chain_separator +
                                 SequenceConstants.monomer_join, '[]')
        nres = len(residues)
        for num_res, res in enumerate(residues):

            if num_res < nres - 1:
                if self.__get_monomer_prop('m_chainID', num_res) == \
                        self.__get_monomer_prop('m_chainID', num_res + 1):
                    # Read attachment points directly from the mol via isotope labels
                    # (CHUCKLES convention) rather than the sidecar index property.
                    mol1 = self.__get_monomer_prop('m_romol', num_res)
                    mol2 = self.__get_monomer_prop('m_romol', num_res + 1)
                    r1_attach = _attachment_idx(mol1, 2)  # slot 2 = R2 = C-term exit
                    r2_attach = _attachment_idx(mol2, 1)  # slot 1 = R1 = N-term entry
                    if r1_attach is None or r2_attach is None:
                        raise ValueError(
                            f"Cannot form backbone R2→R1 bond between residues "
                            f"{num_res} and {num_res + 1}: one or both monomers "
                            f"lack the required attachment point. Check that "
                            f"terminal-only monomers (e.g. ac, am, fmoc) are "
                            f"not placed mid-chain.")
                    else:
                        _check_parsed_connection(
                            mol1, r1_attach, 2, mol2, r2_attach, 1,
                            leaving_groups1=self.__get_monomer_prop('m_Rgroups', num_res),
                            leaving_groups2=self.__get_monomer_prop('m_Rgroups', num_res + 1),
                            bond_label=f'backbone {num_res}→{num_res+1}',
                            warning_sink=self._warning_sink)
                        self.__add_bond(
                            num_res, r1_attach, num_res + 1, r2_attach, slot1=2, slot2=1)

            # Find connections with other chains.
            # Bond IDs may be plain integers ("1") or !n markers ("!1").
            match = re.findall(r"\((!\w+|\d+),(\d+)\)", res)
            for i, m_value in enumerate(match):
                b_idx, rgrp_idx = m_value[0], int(m_value[1])
                mol = self.__get_monomer_prop('m_romol', num_res)
                attach = _attachment_idx(mol, rgrp_idx)
                if attach is None:
                    raise ValueError(
                        f"Residue {num_res} has no R{rgrp_idx} attachment "
                        f"point but BILN specifies a bond using R{rgrp_idx}.")
                bond_info.append([b_idx, num_res, rgrp_idx, attach])
                bond_info_helm.append([b_idx, num_res, rgrp_idx])

            # For multiple branching (two crosslinks on the same residue)
            match = re.findall(r"\((!\w+|\d+),(\d+);(!\w+|\d+),(\d+)\)", res)
            for i, m_value in enumerate(match):
                b1_idx, rgrp1_idx = m_value[0], int(m_value[1])
                b2_idx, rgrp2_idx = m_value[2], int(m_value[3])
                mol = self.__get_monomer_prop('m_romol', num_res)
                attach1 = _attachment_idx(mol, rgrp1_idx)
                attach2 = _attachment_idx(mol, rgrp2_idx)
                if attach1 is None:
                    raise ValueError(
                        f"Residue {num_res} has no R{rgrp1_idx} attachment "
                        f"point but BILN specifies a bond using R{rgrp1_idx}.")
                if attach2 is None:
                    raise ValueError(
                        f"Residue {num_res} has no R{rgrp2_idx} attachment "
                        f"point but BILN specifies a bond using R{rgrp2_idx}.")
                bond_info.append([b1_idx, num_res, rgrp1_idx, int(attach1)])
                bond_info.append([b2_idx, num_res, rgrp2_idx, int(attach2)])
                bond_info_helm.append([b1_idx, num_res, rgrp1_idx])
                bond_info_helm.append([b2_idx, num_res, rgrp2_idx])

        nbonds = int(len(bond_info))

        if (nbonds % 2) == 1:
            raise ValueError("Must be an even number of extra bonds.")

        if nbonds > 0:
            # Group bond_info entries by bond identifier (str — may be "1" or "!1").
            grouped: dict = {}
            for entry in bond_info:
                grouped.setdefault(entry[0], []).append(entry)

            bond = []
            for bid, entries in grouped.items():
                if len(entries) != 2:
                    raise ValueError(
                        f"Bond identifier {bid!r} has {len(entries)} endpoint(s); "
                        "expected exactly 2 (one on each partner residue).")
                bond.append(entries)

            for bondx in bond:

                if bondx[0][1] > bondx[1][1]:
                    bondx[0], bondx[1] = bondx[1], bondx[0]

                m1, at1, m2, at2 = bondx[0][1], bondx[0][3], bondx[1][1], bondx[1][3]
                slot1 = bondx[0][2]  # 1-based slot for m1's attachment
                slot2 = bondx[1][2]  # 1-based slot for m2's attachment
                mol1 = self.__get_monomer_prop('m_romol', m1)
                mol2 = self.__get_monomer_prop('m_romol', m2)
                _check_parsed_connection(
                    mol1, at1, slot1, mol2, at2, slot2,
                    leaving_groups1=self.__get_monomer_prop('m_Rgroups', m1),
                    leaving_groups2=self.__get_monomer_prop('m_Rgroups', m2),
                    bond_label=f'crosslink bond-id {bondx[0][0]!r} '
                               f'(R{bondx[0][2]} of residue {m1} <-> '
                               f'R{bondx[1][2]} of residue {m2})',
                    warning_sink=self._warning_sink)
                self.__add_bond(m1, at1, m2, at2, slot1=slot1, slot2=slot2)

        self.__only_unique_bonds()


    def __add_bond(self, m_id1, atom1, m_id2, atom2, slot1=None, slot2=None):
        """
        Add an entry in the bond list

        :param m_id1: index of monomer 1
        :param atom1: index of bond atom in monomer 1
        :param m_id2: index of monomer 2
        :param atom2: index of bond atom in monomer 2
        :param slot1: 1-based R-group slot for atom1 (optional; avoids ambiguity
                      when an atom neighbours multiple dummies, e.g. backbone N
                      which has both [1*] and [3*]).
        :param slot2: 1-based R-group slot for atom2 (optional; same rationale).
        """
        slot1 = slot1 if slot1 is not None else _slot_for_attachment(
            self.s_monomers[m_id1]['m_romol'], atom1)
        slot2 = slot2 if slot2 is not None else _slot_for_attachment(
            self.s_monomers[m_id2]['m_romol'], atom2)
        endpoints = ((m_id1, slot1), (m_id2, slot2))
        if endpoints[0] == endpoints[1]:
            raise ValueError("A bond cannot join the same attachment slot to itself.")
        for monomer_id, slot in endpoints:
            if (monomer_id, slot) in self._used_slots:
                raise ValueError(
                    f"Residue {monomer_id} R{slot} is already used by another bond.")
        self._used_slots.update(endpoints)
        self.s_nbonds += 1
        entry = [m_id1, atom1, m_id2, atom2]
        if slot1 is not None:
            entry.extend([slot1, slot2])
        self.s_bonds.append(entry)

    def get_monomer(self, idx):
        """
        Return the dictionary of ith monomer by index

        :param idx: index of monomer
        :type idx: int
        """

        mon = None
        try:
            mon = self.s_monomers[idx]
        except IndexError:
            warnings.warn(f'Cannot access monomer index {idx} in sequence')

        return mon

    def __get_monomer_prop(self, property_val, idx):
        """
        Return a property of the i-th monomer

        :param property_val: which property to fetch
        :type property_val: str
        :param idx: index of monomer
        :type idx: int
        :return: desired property
        """

        m_prop = None
        try:
            m_dict = self.get_monomer(idx)
        except ValueError:
            m_prop = None

        try:
            m_prop = m_dict[property_val]
        except ValueError:
            warnings.warn(
                f'Cannot access property {property_val} in monomer {idx}')

        return m_prop

    def length(self):
        """
        Return the length of the sequence
        """
        return len(self.s_monomers)

    def __len__(self):
        """
        Return the length of the sequence.
        """
        return self.length()

    def __only_unique_bonds(self):
        """
        In __parseBILNBonds it might happen that a bond is added twice
        (once from the amide section, once if it is)
        in addition set explicitly by (1,2) sets.
        Here, we filter these bonds out

        :return: filtered bond list with unique bonds only
        """
        seen = set()
        unique = []
        for b in self.s_bonds:
            key = tuple(b)
            if key not in seen:
                seen.add(key)
                unique.append(b)
        self.s_bonds = unique

    def is_valid(self):
        """Flag if the initialized pyPept.Sequence is valid.
        This includes check of R-groups present in the monomer dictionary,
        and R-group connectivity when explicitly given.

        :return: bool
        """
        return self.__is_valid

    @classmethod
    def validate(cls, biln, path=SequenceConstants.def_path,
                 monomer_lib=SequenceConstants.def_lib_filename, fmt=None):
        """Dry-run parse and chemistry-check a BILN/CABILN string.

        Constructs the Sequence, captures all warnings and errors, and returns
        a ValidationReport without raising.  Call report.build() to get the
        Sequence object if report.ok is True.

        :param biln: BILN or CABILN sequence string
        :param path: monomer library package path (default: built-in)
        :param monomer_lib: monomer SDF filename (default: monomers.sdf)
        :param fmt: pass ``'biln'`` to accept old BILN crosslink notation.
        :returns: ValidationReport
        """
        caught_warns = []
        errors = []
        seq = None

        try:
            seq = cls(biln, path=path, monomer_lib=monomer_lib, fmt=fmt,
                      warning_sink=caught_warns.append)
        except ValueError as exc:
            errors.append(str(exc))

        bonds = []
        if seq is not None:
            for entry in seq.s_bonds:
                m1, m2 = entry[0], entry[2]
                slot1 = entry[4] if len(entry) > 4 else None
                slot2 = entry[5] if len(entry) > 5 else None
                bonds.append((m1, slot1, m2, slot2))

        ok = seq is not None and not errors
        report = ValidationReport(
            biln=biln, ok=ok, bonds=bonds,
            warnings=caught_warns, errors=errors,
        )
        report._seq = seq
        return report

    def __bool__(self):
        """
        Returns if a pyPept.Sequence object has been constructed with a valid
        BILN.

        :seealso: is_valid

        >>> s = Sequence('P-E-P')
        >>> bool(s)
        True
        >>> s = Sequence('P-nonsense-P')
        >>> bool(s)
        False
        """
        return self.is_valid()


def get_monomer_info(path):
    """Return a detached monomer table from the shared, versioned library."""
    from pyPept.monomer_store import monomer_table

    return monomer_table(path)
