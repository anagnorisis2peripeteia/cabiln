"""Expand inline SMILES into local monomer definitions without changing the library.

Token expansion protects SMILES punctuation before notation lowering. Definition
construction retains the existing synthetic-site and terminal-group conventions.
"""

import hashlib
import re

from pyPept.source import sub as _source_sub, synthetic as _source_synthetic
from pyPept.structure import parse_template_smiles

_SYN_NCAP_RE = re.compile(r'<([^<>]*)>_')
_SYN_CCAP_RE = re.compile(r'_<([^<>]*)>')
_SYN_AA_RE = re.compile(r'<([^<>]*)>')


def _build_synthetic_aa(sidechain_smi):
    """Build an alpha-AA ROMol template from a sidechain SMILES (or a full
    residue template).

    Two input forms accepted:
      1. Sidechain-only (no `[1*]` / `[2*]`): wrapped as
         `[1*]NC(<sc>)C(=O)[2*]`.  Cα stereo is unset.
      2. Full template (contains both `[1*]` and `[2*]`): used as-is.
         This lets the walker preserve Cα stereo by emitting e.g.
         `<[1*]N[C@@H](C)C([2*])=O>` for L-alanine-like residues.

    Layout convention: R1 (isotope 1) on backbone N, R2 (isotope 2) on
    carbonyl C, matching the SDF library.  Returns None on parse error.
    Empty / whitespace `sidechain_smi` produces a glycine equivalent.
    """
    if sidechain_smi is None:
        sidechain_smi = ''
    sc = sidechain_smi.strip()
    if not sc:
        smi = '[1*]NCC(=O)[2*]'
    elif re.search(r'\[1\*[:\]]', sc) and re.search(r'\[2\*[:\]]', sc):
        # Full template form: accept [1*], [1*:n], [2*], [2*:n]. Atom-map digits
        # after the isotope are in-band flags (e.g. [2*:99] = aldehyde terminus).
        smi = sc
    else:
        smi = f'[1*]NC({sc})C(=O)[2*]'
    try:
        return parse_template_smiles(smi)
    except Exception:
        return None


def _build_synthetic_cap(cap_smi, side):
    """Build a synthetic cap ROMol from a SMILES.

    side='n':  N-cap — the cap's C(=O) bonds to the next residue's R1
               via standard backbone_amide.  Output has R2 only.
    side='c':  C-cap — the cap's N bonds to the previous residue's R2.
               Output has R1 only.

    Two input forms accepted (same as _build_synthetic_aa):
      - Bare SMILES (no `[1*]` / `[2*]`): walker is responsible for
        placing the dummy atom; we wrap minimally.  For N-cap, append
        `[2*]` to the LAST atom; for C-cap, prepend `[1*]` to the first.
      - Full template containing the appropriate `[N*]`: used as-is.
    """
    if cap_smi is None:
        cap_smi = ''
    s = cap_smi.strip()
    if not s:
        return None
    if side == 'n':
        if '[2*]' in s:
            smi = s
        else:
            smi = f'{s}[2*]'
    elif side == 'c':
        if '[1*]' in s:
            smi = s
        else:
            smi = f'[1*]{s}'
    else:
        return None
    try:
        return parse_template_smiles(smi)
    except Exception:
        return None


def _synthetic_aa_symbol(sidechain_smi):
    """Stable synthetic-AA monomer symbol from sidechain SMILES."""
    h = hashlib.md5((sidechain_smi or '').encode()).hexdigest()[:10]
    return f'__SYN_{h}'


def _synthetic_cap_symbol(cap_smi, side):
    """Stable synthetic cap monomer symbol from cap SMILES + side ('n' / 'c')."""
    h = hashlib.md5(f'{side}:{cap_smi or ""}'.encode()).hexdigest()[:10]
    return f'__SYN{side.upper()}CAP_{h}'


def _infer_synth_chem_type(mol, slot):
    """Return (chem_type, leaving_group) for an R-group dummy in a synthetic
    monomer template.  Slot 1 → backbone_n, slot 2 → backbone_c, slot 3+ →
    inferred from the dummy's neighbour atom (sidechain anchor).

    The default backbone-C leaving group is [OH] (carboxylic acid terminus).
    To emit a peptide-aldehyde terminus instead (C(=O)H), the caller writes the
    slot-2 dummy as [2*:99] in SMILES — atom-map 99 is an in-band marker that
    overrides _lg to [H] without changing the CABILN grammar."""
    if slot == 1:
        return ('backbone_n', '[H]')
    if slot == 2:
        for atom in mol.GetAtoms():
            if (atom.GetAtomicNum() == 0 and atom.GetIsotope() == 2
                    and atom.GetAtomMapNum() == 99):
                return ('backbone_c', '[H]')
        return ('backbone_c', '[OH]')
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == slot:
            nbs = atom.GetNeighbors()
            if not nbs:
                return (None, None)
            anchor = nbs[0]
            sym = anchor.GetSymbol()
            if sym == 'S':
                return ('thiol', '[H]')
            if sym == 'Se':
                return ('selenol', '[H]')
            if sym == 'N':
                # Could be amine or guanidinium etc. Default to amine_primary.
                return ('amine_primary', '[H]')
            if sym == 'O':
                # Phenolic O if anchor is bonded to aromatic C; else aliphatic OH
                for nb2 in anchor.GetNeighbors():
                    if nb2.GetIsAromatic():
                        return ('aryl_phenol_o', '[H]')
                return ('hydroxyl', '[H]')
            if sym == 'C':
                for nb2 in anchor.GetNeighbors():
                    if nb2.GetAtomicNum() == 8:
                        b = mol.GetBondBetweenAtoms(anchor.GetIdx(), nb2.GetIdx())
                        if b and b.GetBondTypeAsDouble() == 2.0:
                            return ('carboxyl', '[OH]')
                for nb2 in anchor.GetNeighbors():
                    if nb2.GetSymbol() in ('Cl', 'Br', 'I'):
                        return ('alkyl_halide_c', None)
                return ('carbon', '[H]')
            return (None, None)
    return (None, None)


def expand_synthetic_tokens(source):
    """Return protected source and the AA/cap definitions used by its occurrences."""
    amino_acids, caps = {}, {}

    def amino_acid(match):
        smiles = match.group(1)
        symbol = _synthetic_aa_symbol(smiles)
        amino_acids[symbol] = smiles
        return _source_synthetic(match, symbol)

    def cap(match, side):
        smiles = match.group(1)
        symbol = _synthetic_cap_symbol(smiles, side)
        symbol = symbol + '_' if side == 'n' else '_' + symbol
        caps[symbol] = (smiles, side)
        return _source_synthetic(match, symbol)

    # Caps must consume their underscore before the bare AA pattern can match.
    source = _source_sub(_SYN_NCAP_RE, lambda match: cap(match, 'n'), source)
    source = _source_sub(_SYN_CCAP_RE, lambda match: cap(match, 'c'), source)
    return _source_sub(_SYN_AA_RE, amino_acid, source), amino_acids, caps


def add_synthetic_monomers(table, amino_acids, caps):
    """Add local definitions to a detached library table, preserving existing rows."""
    def add(symbol, mol, name, kind, groups, chemistry):
        mol.SetProp('m_chem_types', chemistry)
        table.loc[symbol] = {
            'm_Rgroups': groups, 'm_abbr': symbol, 'm_name': name,
            'm_type': kind, 'm_subtype': 'synthetic',
            'm_chem_types': chemistry, 'ID': symbol, 'm_romol': mol,
        }

    for symbol, smiles in amino_acids.items():
        if symbol in table.index:
            continue
        mol = _build_synthetic_aa(smiles)
        if mol is None:
            raise ValueError(f"Synthetic AA token '<{smiles}>' contains invalid SMILES.")
        slots = sorted({atom.GetIsotope() for atom in mol.GetAtoms()
                        if atom.GetAtomicNum() == 0 and atom.GetIsotope() >= 1})
        groups = [''] * max(slots, default=2)
        chemistry = []
        for slot in slots:
            chem_type, leaving = _infer_synth_chem_type(mol, slot)
            if chem_type is not None:
                groups[slot - 1] = leaving if leaving is not None else ''
                chemistry.append(f'{slot}:{chem_type}')
        add(symbol, mol, f'synthetic_aa<{smiles}>', 'aa', groups, ','.join(chemistry))

    for symbol, (smiles, side) in caps.items():
        if symbol in table.index:
            continue
        mol = _build_synthetic_cap(smiles, side)
        if mol is None:
            raise ValueError(f"Synthetic {side}-cap token contains invalid SMILES: {smiles!r}")
        groups, chemistry = (['', '[OH]'], '2:backbone_c') if side == 'n' else (['[H]', ''], '1:backbone_n')
        add(symbol, mol, f'synthetic_{side}cap<{smiles}>', 'cap', groups, chemistry)
