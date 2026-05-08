"""
Monomer pre-activation pipeline for pyPept.

Converts full monomer SMILES (in any notation) to pre-activated CHUCKLES
fragments with isotope-labelled dummy atoms and SMILES leaving group metadata.

Supports auto-detection of:
  - Backbone R1 (N-terminal amine) and R2 (C-terminal carboxyl) via a
    graph-topology rule: the N and carboxyl C pair with the shortest bond
    path form the backbone.  Works for α-, β-, γ-amino acids and any chain
    length without hard-coded stereo or H-count assumptions.
  - Sidechain R3+ wherever chemistry makes sense (thiol, amine, carboxyl,
    hydroxyl, aromatic NH) — following the principle that any chemically
    viable bond site should be available; whether it is *used* is a
    sequence-level decision.

Each slot carries a chem_type string (e.g. 'backbone_n', 'thiol',
'amine_primary') that the assembly layer uses to look up the correct
SMIRKS reaction entry in reactions.yaml without any if-branches.

Leaving group derivation is auto-inferred via SMARTS and stored as SMILES
fragments (e.g. '[H]', '[OH]') used by Sequence.__remove_rgroups() to
restore free termini when a slot is left unbonded.

Usage:
    from pyPept.interfaces.monomer_pipeline import pre_activate

    result = pre_activate("N[C@@H](CS)C(=O)O")
    # result.chuckles   = '[1*]N([3*])[C@@H](CS[4*])C([2*])=O'
    # result.leaving    = {1: '[H]', 2: '[OH]', 3: '[H]', 4: '[H]'}
    # result.chem_types = {1: 'backbone_n', 2: 'backbone_c', 3: 'backbone_n_mod', 4: 'thiol'}

    # Legacy tuple unpacking still works:
    chuckles, leaving, chem_types, err = pre_activate("N[C@@H](CS)C(=O)O")
    # err is always None — failures now raise ActivationError

From publication: pyPept: a python library to generate atomistic 2D and 3D
representations of peptides, Journal of Cheminformatics, 2023.
Updated 2025.
"""

__credits__ = ["J.B. Brown", "Thomas Fox", "Cameron Beesley"]
__license__ = "MIT"

import csv
import re
import warnings
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import PandasTools, SDWriter, rdDepictor

class ActivationError(ValueError):
    """Raised when pre_activate cannot generate a valid CHUCKLES fragment."""


@dataclass
class ActivationResult:
    """Result of pre_activate(). Iterable for backward-compatible tuple unpacking."""
    chuckles: str
    leaving: dict
    chem_types: dict

    def __iter__(self):
        yield self.chuckles
        yield self.leaving
        yield self.chem_types
        yield None  # err slot — always None; failures raise ActivationError

    def __getitem__(self, idx):
        return list(self)[idx]


# Suppress phantom RDKit H-removal warnings that fire when dummy atoms are added
# adjacent to atoms whose implicit-H count is then recalculated by SanitizeMol.
RDLogger.DisableLog('rdApp.warning')


# ── Standard amino acid lookup (FASTA / BILN token → SMILES) ─────────────
# Used by normalize_input() to convert single-letter or token-name inputs.
_AA_SMILES = {
    # FASTA single-letter codes
    'G': 'NCC(=O)O',
    'A': 'N[C@@H](C)C(=O)O',
    'V': 'N[C@@H](C(C)C)C(=O)O',
    'L': 'N[C@@H](CC(C)C)C(=O)O',
    'I': 'N[C@@H]([C@@H](C)CC)C(=O)O',
    'P': 'OC(=O)[C@@H]1CCCN1',
    'F': 'N[C@@H](Cc1ccccc1)C(=O)O',
    'W': 'N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O',
    'M': 'N[C@@H](CCSC)C(=O)O',
    'S': 'N[C@@H](CO)C(=O)O',
    'T': 'N[C@@H]([C@@H](O)C)C(=O)O',
    'C': 'N[C@@H](CS)C(=O)O',
    'Y': 'N[C@@H](Cc1ccc(O)cc1)C(=O)O',
    'H': 'N[C@@H](Cc1cnc[nH]1)C(=O)O',
    'D': 'N[C@@H](CC(=O)O)C(=O)O',
    'E': 'N[C@@H](CCC(=O)O)C(=O)O',
    'N': 'N[C@@H](CC(=O)N)C(=O)O',
    'Q': 'N[C@@H](CCC(=O)N)C(=O)O',
    'K': 'N[C@@H](CCCCN)C(=O)O',
    'R': 'N[C@@H](CCCNC(=N)N)C(=O)O',
    # Common BILN / token names for NCAAs and D-AAs
    'DAla': 'N[C@H](C)C(=O)O',
    'DVal': 'N[C@H](C(C)C)C(=O)O',
    'DLeu': 'N[C@H](CC(C)C)C(=O)O',
    'DIle': 'N[C@H]([C@H](C)CC)C(=O)O',
    'DPhe': 'N[C@H](Cc1ccccc1)C(=O)O',
    'DTrp': 'N[C@H](Cc1c[nH]c2ccccc12)C(=O)O',
    'DSer': 'N[C@H](CO)C(=O)O',
    'DThr': 'N[C@H]([C@H](O)C)C(=O)O',
    'DCys': 'N[C@H](CS)C(=O)O',
    'DTyr': 'N[C@H](Cc1ccc(O)cc1)C(=O)O',
    'DHis': 'N[C@H](Cc1cnc[nH]1)C(=O)O',
    'DAsp': 'N[C@H](CC(=O)O)C(=O)O',
    'DGlu': 'N[C@H](CCC(=O)O)C(=O)O',
    'DAsn': 'N[C@H](CC(=O)N)C(=O)O',
    'DGln': 'N[C@H](CCC(=O)N)C(=O)O',
    'DLys': 'N[C@H](CCCCN)C(=O)O',
    'DArg': 'N[C@H](CCCNC(=N)N)C(=O)O',
    'Aib':  'CC(N)(C)C(=O)O',
    'Orn':  'N[C@@H](CCCN)C(=O)O',
    'Hyp':  'OC(=O)[C@@H]1C[C@@H](O)CN1',
    'Nle':  'N[C@@H](CCCC)C(=O)O',
    'Nva':  'N[C@@H](CCC)C(=O)O',
    'Sec':  'N[C@@H](C[SeH])C(=O)O',
}

# CHUCKLES pattern: contains [n*] dummy atoms (isotope-labelled wildcards)
_CHUCKLES_RE = re.compile(r'\[\d+\*\]')

# SMILES heuristic: contains chemistry characters not found in token names
_SMILES_CHARS = set('()[]=#@+\\/')


def normalize_input(raw):
    """
    Convert any supported monomer input notation to SMILES.

    Accepted formats:
      - SMILES: any string containing SMILES chemistry characters
      - FASTA: single uppercase letter (A, G, C, ...)
      - BILN token / NCAA name: e.g. DAla, Hyp, Aib
      - CHUCKLES: pre-activated fragment — returned as-is (caller must detect)

    :param raw: raw input string from the CSV `input` column.
    :returns: (smiles_or_chuckles, is_chuckles, error) where error is None on success.
    """
    raw = raw.strip()
    if not raw:
        return None, False, 'Empty input'

    # CHUCKLES: already pre-activated (has [n*] dummies) — signal to caller
    if _CHUCKLES_RE.search(raw):
        return raw, True, None

    # Explicit SMILES: contains SMILES-specific characters
    if _SMILES_CHARS & set(raw):
        mol = Chem.MolFromSmiles(raw)
        if mol is None:
            return None, False, f"Cannot parse SMILES: {raw!r}"
        return raw, False, None

    # FASTA single-letter or BILN token name — look up in table
    candidate = _AA_SMILES.get(raw)
    if candidate:
        return candidate, False, None

    # Last attempt: try RDKit anyway (catches unusual SMILES without special chars)
    mol = Chem.MolFromSmiles(raw)
    if mol is not None:
        return raw, False, None

    return None, False, (
        f"Unrecognised input {raw!r}. Provide SMILES, a FASTA letter, "
        f"a known BILN token (DAla, Hyp, ...), or pre-filled CHUCKLES. "
        f"HELM is not supported as a per-monomer input — use SMILES instead."
    )


# ── Leaving group inference rules ─────────────────────────────────────────
# Each entry: (SMARTS with attachment atom mapped :1, leaving_group_smiles)
# First match wins; more specific rules must come first.
LEAVING_GROUP_RULES = [
    # Halide LGs (reagent form) — checked first so acid chlorides beat free acids
    ('[CX3:1](=O)[ClX1]',        '[Cl]'),  # acid chloride
    ('[CX3:1](=O)[BrX1]',        '[Br]'),  # acid bromide
    ('[SX4:1](=O)(=O)[ClX1]',    '[Cl]'),  # sulfonyl chloride (TsCl, MsCl, Pbf-Cl)
    ('[CX4:1][IX1]',              '[I]'),   # alkyl iodide (MeI)
    ('[CX4:1][BrX1]',            '[Br]'),  # alkyl bromide (EtBr, BnBr)
    ('[CX4:1][ClX1]',            '[Cl]'),  # alkyl chloride (TrtCl, AcmCl)
    # Free-form LGs (condensation / HATU coupling)
    ('[CX3H1:1](=O)',      '[H]'),   # aldehyde C → H leaves
    ('[CX3:1](=O)[OX2H1]', '[OH]'),  # carboxylic acid C → OH leaves
    ('[SX4:1](=O)(=O)[OX2H1]', '[OH]'),  # sulfonyl-OH S → OH leaves
    ('[NX3;H3:1]',         '[H]'),   # unsubstituted amine (NH3 cap) → H leaves
    ('[NX3;H2:1]',         '[H]'),   # primary amine N → H leaves
    ('[NX3;H1;R:1]',       '[H]'),   # proline-type ring N → H leaves
    ('[NX3;H1:1]',         '[H]'),   # secondary amine N (N-methyl backbone) → H leaves
    ('[SX2H1:1]',          '[H]'),   # thiol S → H leaves
    ('[SeX2H1:1]',         '[H]'),   # selenol Se → H leaves
    ('[OX2H1:1]',          '[H]'),   # hydroxyl O → H leaves
    # Carbon nucleophile (Grignard) — last resort
    ('[cH1:1]',            '[H]'),   # aromatic C-H
    ('[CX4;H1,H2,H3:1]',  '[H]'),   # sp3 C-H
]

# Compiled once at import time.
_LG_PATTERNS = [(Chem.MolFromSmarts(s), lg) for s, lg in LEAVING_GROUP_RULES]


# ── Backbone detection — topology patterns ────────────────────────────────
# Atom-type patterns for the graph-distance backbone search.
# N: trivalent N with at least one H (covers primary amine, Pro ring N,
#    secondary N, but excludes tertiary amines like NMe2 which can't donate H).
# C: carboxyl carbonyl carbon (C(=O)OH).
_BB_N_PAT        = Chem.MolFromSmarts('[NX3;H1,H2,H3]')
_BB_COOH_PAT     = Chem.MolFromSmarts('[CX3:1](=O)[OX2H1]')
_BB_ALDEHYDE_PAT = Chem.MolFromSmarts('[CX3H1:1](=O)')
_BB_ALCOHOL_PAT  = Chem.MolFromSmarts('[OX2H1:1][CX4]')
_BB_LACTONE_PAT  = Chem.MolFromSmarts('[CX3:1](=O)[OX2;R]')

# Sidechain rules derived from the unified registry in reaction_library.
# Uses pre_smarts column (H≥1 required — there must be an H to replace with dummy).
# All entries including label_only are included to reserve slots in pre_activate.
from pyPept.interfaces.reaction_library import _CHEM_TYPE_REGISTRY
_SIDECHAIN_RULES = [(pre_smarts, lg, ct, lo) for ct, pre_smarts, lg, _infer, lo in _CHEM_TYPE_REGISTRY]

_SC_PATTERNS  = [(Chem.MolFromSmarts(s), lg, ct, lo) for s, lg, ct, lo in _SIDECHAIN_RULES]
_SECOND_H_PAT = Chem.MolFromSmarts('[NX3;H2:1]')


# ── Public API ────────────────────────────────────────────────────────────

def find_backbone_slots(mol):
    """
    Identify R1 (N-terminal) and R2 (C-terminal) attachment atoms via a
    graph-topology rule: the (N, carboxyl-C) pair with the shortest bond
    path is the backbone pair.

    This replaces the earlier 4-SMARTS approach and handles α-, β-, γ-amino
    acids, D-amino acids, and any chain length without stereo or H-count
    constraints.

    :param mol: RDKit mol with explicit H (call AddHs first).
    :returns: {0: N_atom_idx, 1: carbonyl_C_atom_idx} or None if not found.
    """
    n_idxs = [m[0] for m in mol.GetSubstructMatches(_BB_N_PAT)]
    c_idxs = [m[0] for m in mol.GetSubstructMatches(_BB_COOH_PAT)]

    # Fallback: if no COOH, try aldehyde → alcohol → lactone as C-terminal analogue
    if not c_idxs:
        c_idxs = [m[0] for m in mol.GetSubstructMatches(_BB_ALDEHYDE_PAT)]
    if not c_idxs:
        c_idxs = [m[0] for m in mol.GetSubstructMatches(_BB_ALCOHOL_PAT)]
    if not c_idxs:
        c_idxs = [m[0] for m in mol.GetSubstructMatches(_BB_LACTONE_PAT)]

    if not n_idxs or not c_idxs:
        return None

    best_n, best_c, best_dist = None, None, float('inf')
    for n_idx in n_idxs:
        for c_idx in c_idxs:
            path = Chem.GetShortestPath(mol, n_idx, c_idx)
            if path and len(path) - 1 < best_dist:
                best_dist = len(path) - 1
                best_n, best_c = n_idx, c_idx

    return {1: best_n, 2: best_c} if best_n is not None else None


_CAP_SULFONYL_PAT = Chem.MolFromSmarts('[SX4:1](=O)(=O)[OX2H1]')
_CAP_ALCOHOL_PAT  = Chem.MolFromSmarts('[OX2H1:1][CX4]')
# Reagent-form halide patterns (electrophilic caps entered with their LG)
_CAP_ACYL_HALIDE_PAT    = Chem.MolFromSmarts('[CX3:1](=O)[Cl,Br]')
_CAP_SULFONYL_CL_PAT    = Chem.MolFromSmarts('[SX4:1](=O)(=O)[ClX1]')
_CAP_ALKYL_HALIDE_PAT   = Chem.MolFromSmarts('[CX4:1][Cl,Br,I]')
# Carbon nucleophile caps (Grignard / organometallic — last resort)
_CAP_BENZYLIC_CH_PAT    = Chem.MolFromSmarts('[CX4;H1,H2,H3:1][c]')
_CAP_TERTIARY_CH_PAT    = Chem.MolFromSmarts('[CX4;H1:1]([CX4])([CX4])[CX4]')
_CAP_AROMATIC_CH_PAT    = Chem.MolFromSmarts('[cH1:1]')


def _pick_alpha_n(mol, n_idxs):
    """Among multiple amine N atoms, pick the one closest to a C=O or C-OH."""
    targets = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 8:
                    targets.append(atom.GetIdx())
                    break
    if not targets:
        return n_idxs[0]
    best_n, best_dist = n_idxs[0], float('inf')
    for n_idx in n_idxs:
        for c_idx in targets:
            path = Chem.GetShortestPath(mol, n_idx, c_idx)
            if path and len(path) - 1 < best_dist:
                best_dist = len(path) - 1
                best_n = n_idx
    return best_n


def find_cap_slots(mol):
    """
    Fallback for single-ended cap monomers that lack the full N+COOH backbone.

    Handles both free-form (COOH, sulfonyl-OH, amine, alcohol) and
    reagent-form (acid chloride, sulfonyl chloride, alkyl halide) caps.
    All electrophilic caps → R2 (slot 2). Nucleophilic caps → R1 (slot 1).

    :param mol: RDKit mol with explicit H.
    :returns: (slots_dict, chem_types_dict) or (None, None) if unresolvable.
    """
    c_idxs = [m[0] for m in mol.GetSubstructMatches(_BB_COOH_PAT)]
    n_idxs = [m[0] for m in mol.GetSubstructMatches(_BB_N_PAT)]
    s_idxs = [m[0] for m in mol.GetSubstructMatches(_CAP_SULFONYL_PAT)]

    # Multi-arm crosslinker: 3+ alkyl halides → defer to sidechain-only fallback
    if not n_idxs and not c_idxs and not s_idxs:
        alk_multi = [m[0] for m in mol.GetSubstructMatches(_CAP_ALKYL_HALIDE_PAT)]
        if len(alk_multi) >= 3:
            return None, None

    # Free-form electrophilic caps
    if c_idxs and not n_idxs:
        return {2: c_idxs[0]}, {2: 'backbone_c'}
    if s_idxs and not n_idxs and not c_idxs:
        return {2: s_idxs[0]}, {2: 'backbone_c'}

    # Reagent-form electrophilic caps (halide LG still attached)
    if not n_idxs and not c_idxs and not s_idxs:
        acyl_h = [m[0] for m in mol.GetSubstructMatches(_CAP_ACYL_HALIDE_PAT)]
        if acyl_h:
            return {2: acyl_h[0]}, {2: 'backbone_c'}
        sul_cl = [m[0] for m in mol.GetSubstructMatches(_CAP_SULFONYL_CL_PAT)]
        if sul_cl:
            return {2: sul_cl[0]}, {2: 'backbone_c'}
        alk_h = [m[0] for m in mol.GetSubstructMatches(_CAP_ALKYL_HALIDE_PAT)]
        if alk_h and len(alk_h) < 3:
            return {2: alk_h[0]}, {2: 'backbone_c'}

    # Nucleophilic caps — but prefer electrophilic alkyl halide if present
    # (e.g. acm: Cl-CH2-NH-COCH3 has both amine and alkyl chloride)
    if n_idxs and not c_idxs:
        alk_h = [m[0] for m in mol.GetSubstructMatches(_CAP_ALKYL_HALIDE_PAT)]
        if alk_h and len(alk_h) < 3:
            return {2: alk_h[0]}, {2: 'backbone_c'}
        best = _pick_alpha_n(mol, n_idxs) if len(n_idxs) > 1 else n_idxs[0]
        return {1: best}, {1: 'backbone_n'}
    if not n_idxs and not c_idxs and not s_idxs:
        o_idxs = [m[0] for m in mol.GetSubstructMatches(_CAP_ALCOHOL_PAT)]
        if o_idxs:
            return {1: o_idxs[0]}, {1: 'backbone_n'}

    # Carbon nucleophile caps (Grignard / organometallic — last resort)
    for pat in (_CAP_BENZYLIC_CH_PAT, _CAP_TERTIARY_CH_PAT, _CAP_AROMATIC_CH_PAT):
        hits = [m[0] for m in mol.GetSubstructMatches(pat)]
        if hits:
            return {1: hits[0]}, {1: 'carbon'}

    return None, None


def find_sidechain_slots(mol, assigned_atoms, start_slot=4):
    """
    Auto-detect sidechain attachment atoms from atoms not already assigned.

    Assigns sequential slots with no gaps starting from start_slot.
    backbone_n_mod is handled by pre_activate before this runs.

    :param mol: RDKit mol with explicit H.
    :param assigned_atoms: atom indices already consumed by backbone slots.
    :param start_slot: first slot number to assign; caller sets this to
                       len(backbone) so numbering is always gap-free.
    :returns: {slot: (attachment_idx, leaving_smiles, chem_type)}.
    """
    slots = {}
    next_slot = start_slot
    seen = set(assigned_atoms)
    pre_scan = frozenset(seen)  # backbone atoms — excluded from second-H pass
    protected = set()  # non-attachment atoms in multi-atom SMARTS — off-limits
    first_ct = {}      # atom_idx -> first assigned chem_type

    for patt, leaving, chem_type, label_only in _SC_PATTERNS:
        matches = sorted(mol.GetSubstructMatches(patt), key=lambda m: m[0])
        for match in matches:
            idx = match[0]
            if idx not in seen and idx not in protected:
                lg = leaving if leaving is not None else infer_leaving_group(mol, idx)
                slots[next_slot] = (idx, lg, chem_type)
                seen.add(idx)
                first_ct[idx] = chem_type
                next_slot += 1
                if not label_only:
                    for other_idx in match[1:]:
                        protected.add(other_idx)

    # Second H on sidechain primary amines only (not backbone atoms).
    for match in mol.GetSubstructMatches(_SECOND_H_PAT):
        idx = match[0]
        if (idx in seen and idx not in pre_scan
                and first_ct.get(idx) == 'amine_primary'):
            slots[next_slot] = (idx, '[H]', 'amine_secondary')
            next_slot += 1

    return slots


def infer_leaving_group(mol, attachment_idx):
    """
    Infer leaving group SMILES for an attachment atom via SMARTS rules.

    :param mol: RDKit mol with explicit H.
    :param attachment_idx: index of the attachment atom.
    :returns: leaving group SMILES string, or None if unrecognised.
    """
    for patt, leaving in _LG_PATTERNS:
        for match in mol.GetSubstructMatches(patt):
            if match[0] == attachment_idx:
                return leaving
    return None


def validate_leaving_group(mol, attachment_idx, leaving_smiles):
    """
    Verify that leaving_smiles is structurally present on the attachment atom.

    :param mol: RDKit mol with explicit H.
    :param attachment_idx: index of the attachment atom.
    :param leaving_smiles: SMILES of the proposed leaving group.
    :returns: (ok: bool, message: str)
    """
    lg_mol = Chem.MolFromSmiles(leaving_smiles)
    if lg_mol is None:
        return False, f"Cannot parse leaving group SMILES '{leaving_smiles}'"

    attach_atom = mol.GetAtomWithIdx(attachment_idx)
    composite = f'[#{attach_atom.GetAtomicNum()}:1]{leaving_smiles}'
    patt = Chem.MolFromSmarts(composite)
    if patt is None:
        return False, f"Cannot build validation SMARTS from '{leaving_smiles}'"

    for match in mol.GetSubstructMatches(patt):
        if match[0] == attachment_idx:
            return True, 'ok'

    neighbors = [f"{nb.GetSymbol()}(H{nb.GetTotalNumHs()})"
                 for nb in attach_atom.GetNeighbors()]
    return False, (
        f"'{leaving_smiles}' not found on "
        f"{attach_atom.GetSymbol()}[{attachment_idx}]. "
        f"Actual neighbors: {neighbors}"
    )


def find_all_backbone_slots(mol):
    """Return ALL (N, COOH) backbone pairings sorted by bond-path length.

    Unlike find_backbone_slots which returns only the shortest-path (alpha)
    assignment, this returns every distinct COOH paired with the nearest N,
    ordered by increasing path length.  Use this to generate β, γ, δ backbone
    orientation registrations for multi-COOH monomers like Asp, Gla.

    :param mol: RDKit mol with explicit H (call AddHs first).
    :returns: list of ({1: n_idx, 2: c_idx}, path_len) sorted by path_len.
              First entry matches find_backbone_slots (shortest = alpha).
    """
    n_idxs = [m[0] for m in mol.GetSubstructMatches(_BB_N_PAT)]
    c_idxs = [m[0] for m in mol.GetSubstructMatches(_BB_COOH_PAT)]
    if not n_idxs or not c_idxs:
        return []

    triples = []
    for n_idx in n_idxs:
        for c_idx in c_idxs:
            path = Chem.GetShortestPath(mol, n_idx, c_idx)
            if path:
                triples.append((len(path) - 1, n_idx, c_idx))
    triples.sort()

    seen_c = set()
    result = []
    for dist, n_idx, c_idx in triples:
        if c_idx not in seen_c:
            seen_c.add(c_idx)
            result.append(({1: n_idx, 2: c_idx}, dist))
    return result


_DIST_SUFFIXES = {2: '', 3: '_b', 4: '_g', 5: '_d'}


def pre_activate_all(smiles):
    """Generate CHUCKLES registrations for every backbone orientation.

    Calls pre_activate once per distinct (N, COOH) pairing found by
    find_all_backbone_slots.  Monomers with only one COOH return a single-
    element list identical to pre_activate.  Monomers with two COOHs (e.g.
    Asp, Glu) return two entries; three COOHs return three, and so on.

    :param smiles: full monomer SMILES.
    :returns: list of (suffix, ActivationResult) sorted by path length.
              Suffix is '' (alpha), '_b' (beta, path 3), '_g' (gamma, path 4),
              '_d' (delta, path 5), or '_xN' for longer paths.
    :raises ActivationError: if no backbone pairings found.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ActivationError(f"Invalid SMILES: {smiles}")
    mol_h = Chem.AddHs(mol)
    pairings = find_all_backbone_slots(mol_h)
    if not pairings:
        raise ActivationError(f"No backbone pairings found in: {smiles}")
    results = []
    for backbone_dict, dist in pairings:
        suffix = _DIST_SUFFIXES.get(dist, f'_x{dist}')
        result = pre_activate(smiles, backbone_indices=backbone_dict)
        results.append((suffix, result))
    return results


def pre_activate(smiles, slot_overrides=None, leaving_overrides=None,
                 backbone_indices=None):
    """
    Convert a full monomer SMILES to a pre-activated CHUCKLES fragment.

    Backbone R1/R2 are always auto-detected via SMARTS backbone patterns.
    Sidechain R3+ are auto-detected wherever a chemically viable bond site
    exists; the sequence author decides whether to use them.

    :param smiles: full monomer SMILES (any stereo/notation, all atoms present).
    :param slot_overrides: {slot: attachment_atom_idx} to add/override slots.
    :param leaving_overrides: {slot: leaving_smiles} to override auto-derived
           leaving groups.  Each override is validated against the mol.
    :param backbone_indices: {1: n_idx, 2: c_idx} to force specific backbone
           atoms, bypassing find_backbone_slots.  Atom indices must correspond
           to the mol AFTER AddHs.  Used by pre_activate_all for β/γ variants.
    :returns: ActivationResult with .chuckles, .leaving, .chem_types.
              Raises ActivationError on failure.
              Also supports legacy tuple unpacking: chuckles, leaving, chem_types, err = pre_activate(...)
              where err is always None (failures raise instead).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ActivationError(f"Invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)

    _is_cap = False
    _sidechain_only = False
    if backbone_indices is not None:
        backbone = dict(backbone_indices)
    else:
        backbone = find_backbone_slots(mol)
    if backbone is None:
        backbone, backbone_chem_types = find_cap_slots(mol)
        _is_cap = True
        if backbone is None:
            # Last resort: sidechain-only molecule (multi-arm crosslinkers)
            _sc_all = find_sidechain_slots(mol, assigned_atoms=set(), start_slot=4)
            if _sc_all:
                backbone = {s: info[0] for s, info in _sc_all.items()}
                backbone_chem_types = {s: info[2] for s, info in _sc_all.items()}
                _sidechain_only = True
            else:
                raise ActivationError(
                    "Cannot identify backbone (N + COOH) or a single cap terminus. "
                    "Provide pre-filled CHUCKLES in the input column."
                )
    else:
        r2_atom = mol.GetAtomWithIdx(backbone[2])
        r2_element = r2_atom.GetAtomicNum()
        if r2_element == 8:
            r2_ct = 'backbone_o'
        elif any(
            nb.GetAtomicNum() == 8 and nb.IsInRing()
            for nb in r2_atom.GetNeighbors()
            if mol.GetBondBetweenAtoms(r2_atom.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 1.0
        ):
            r2_ct = 'lactone_c'
        else:
            r2_ct = 'backbone_c'
        backbone_chem_types = {1: 'backbone_n', 2: r2_ct}
        # backbone_n_mod: second H on backbone N (N-methylation slot).
        # Pro's ring N has only 1 H and is skipped naturally.
        n_atom = mol.GetAtomWithIdx(backbone[1])
        h_count = sum(1 for nb in n_atom.GetNeighbors() if nb.GetAtomicNum() == 1)
        if h_count >= 2:
            slot = max(backbone.keys()) + 1  # 3 for standard AA
            backbone[slot] = backbone[1]
            backbone_chem_types[slot] = 'backbone_n_mod'

    # Exclude backbone atoms from the sidechain scan. The backbone COOH's
    # hydroxyl O is also excluded so downstream patterns (hydroxyl, etc.)
    # don't spuriously fire on it.
    _bb_excluded = set(backbone.values())
    if 2 in backbone:
        _c = mol.GetAtomWithIdx(backbone[2])
        for _nb in _c.GetNeighbors():
            if _nb.GetAtomicNum() == 8:
                _bond = mol.GetBondBetweenAtoms(backbone[2], _nb.GetIdx())
                if _bond.GetBondTypeAsDouble() == 1.0:
                    _bb_excluded.add(_nb.GetIdx())

    if _sidechain_only:
        # All slots already determined — skip further sidechain detection
        sidechain = {}
        sidechain_leaving = {s: _sc_all[s][1] for s in _sc_all}
    else:
        _sc_start = max(backbone.keys()) + 1
        if _is_cap:
            _sc_start = max(_sc_start, 4)
        sidechain = find_sidechain_slots(
            mol,
            assigned_atoms=_bb_excluded,
            start_slot=_sc_start,
        )
    slots = {**backbone,
             **{s: info[0] for s, info in sidechain.items()}}
    sidechain_leaving = {**sidechain_leaving} if _sidechain_only else {s: info[1] for s, info in sidechain.items()}

    chem_types = {
        **backbone_chem_types,
        **{s: info[2] for s, info in sidechain.items()},
    }

    if slot_overrides:
        slots.update(slot_overrides)

    # Resolve leaving groups
    leaving = {}
    for slot, attach_idx in slots.items():
        if leaving_overrides and slot in leaving_overrides:
            lg = leaving_overrides[slot]
            ok, msg = validate_leaving_group(mol, attach_idx, lg)
            if not ok:
                raise ActivationError(f"Slot R{slot} leaving group invalid: {msg}")
        else:
            lg = sidechain_leaving.get(slot) or infer_leaving_group(mol, attach_idx)
            if lg is None:
                nb_desc = [f"{nb.GetSymbol()}(H{nb.GetTotalNumHs()})"
                           for nb in mol.GetAtomWithIdx(attach_idx).GetNeighbors()]
                raise ActivationError(
                    f"Cannot infer leaving group for slot R{slot} "
                    f"({mol.GetAtomWithIdx(attach_idx).GetSymbol()}[{attach_idx}]). "
                    f"Neighbors: {nb_desc}. Provide via leaving_overrides."
                )
        leaving[slot] = lg

    # Build pre-activated fragment.
    # Process slots in order so that when the same atom appears in two slots
    # (e.g. backbone N in R1 and backbone-N modification slot), the second call
    # picks the OTHER H — not the one already claimed by the first slot.
    emol = Chem.RWMol(mol)
    atoms_to_remove = []
    _already_claimed = set()
    for slot, attach_idx in sorted(slots.items()):
        indices = _find_leaving_atoms(emol, attach_idx, leaving[slot],
                                      exclude=_already_claimed)
        atoms_to_remove.extend(indices)
        _already_claimed.update(indices)

    # Add dummies at high indices (before removal to avoid index shifts)
    for slot, attach_idx in sorted(slots.items()):
        dummy = Chem.Atom(0)
        dummy.SetIsotope(slot)
        new_idx = emol.AddAtom(dummy)
        emol.AddBond(attach_idx, new_idx, Chem.BondType.SINGLE)

    for idx in sorted(set(atoms_to_remove), reverse=True):
        emol.RemoveAtom(idx)

    mol_final = Chem.RemoveHs(emol.GetMol())
    try:
        Chem.SanitizeMol(mol_final)
    except Exception as e:
        raise ActivationError(f"Sanitization failed: {e}") from e

    rdDepictor.Compute2DCoords(mol_final)

    return ActivationResult(
        chuckles=Chem.MolToSmiles(mol_final),
        leaving=leaving,
        chem_types=chem_types,
    )


# ── CSV authoring pipeline ────────────────────────────────────────────────

_CSV_COLUMNS = [
    'token', 'input', 'name', 'type', 'synonyms',
    'chuckles',
    'r1_leaving', 'r2_leaving', 'r3_leaving', 'r4_leaving',
    'r5_leaving', 'r6_leaving',
]

_LG_COLS = {
    1: 'r1_leaving', 2: 'r2_leaving', 3: 'r3_leaving',
    4: 'r4_leaving', 5: 'r5_leaving', 6: 'r6_leaving',
}


def derive_monomers(csv_path, rebuild=False):
    """
    Derive CHUCKLES and leaving groups for all monomers in a CSV.

    Pure transform — reads the CSV, runs pre_activate() for any row that
    lacks a valid CHUCKLES (or all 'aa' rows when rebuild=True), and returns
    the updated rows as a list of dicts.  Nothing is written to disk.

    :param csv_path: path to the monomer CSV.
    :param rebuild: if True, re-derive all 'aa' type rows even if CHUCKLES present.
    :returns: (rows, errors) where rows is list[dict] and errors is list[str].
    """
    csv_path = Path(csv_path)
    with open(csv_path, newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))

    errors = []
    for row in rows:
        token = row.get('token', '').strip()
        inp = row.get('input', '').strip()
        existing_chuckles = row.get('chuckles', '').strip()
        row_type = row.get('type', '').strip()
        chem_types_derived = {}

        chuckles_valid = (
            bool(existing_chuckles and Chem.MolFromSmiles(existing_chuckles))
            and not (rebuild and row_type == 'aa')
        )

        if chuckles_valid:
            continue  # keep existing columns as-is

        smiles_or_chuck, is_chuckles, norm_err = normalize_input(inp)
        if norm_err:
            errors.append(f"{token}: {norm_err}")
            warnings.warn(f"Skipping monomer '{token}': {norm_err}")
            continue

        if is_chuckles:
            row['chuckles'] = smiles_or_chuck
            continue

        leaving_overrides = None if (rebuild and row_type == 'aa') else (
            {slot: val for slot, col in _LG_COLS.items()
             if (val := (row.get(col) or '').strip()) and val != 'None'} or None
        )
        try:
            result = pre_activate(smiles_or_chuck, leaving_overrides=leaving_overrides)
        except ActivationError as e:
            errors.append(f"{token}: {e}")
            warnings.warn(f"Skipping monomer '{token}': {e}")
            continue

        row['chuckles'] = result.chuckles
        chem_types_derived = result.chem_types
        for slot, col in _LG_COLS.items():
            row[col] = result.leaving.get(slot, '') or ''
        row['_chem_types'] = chem_types_derived  # transient; not written to CSV

    if errors:
        warnings.warn(f"{len(errors)} monomer(s) failed:\n" + "\n".join(errors))

    return rows, errors


def write_sdf(rows, output_sdf):
    """
    Write a list of monomer row dicts to an SDF file.

    :param rows: list of dicts as returned by derive_monomers().
    :param output_sdf: path for the output SDF.
    :returns: (ok_count, error_list)
    """
    writer = SDWriter(str(output_sdf))
    errors = []

    for row in rows:
        token = row.get('token', '').strip()
        chuckles = row.get('chuckles', '').strip()
        if not chuckles:
            continue

        leaving = {
            slot: (row.get(col) or '').strip() or None
            for slot, col in _LG_COLS.items()
        }
        chem_types_derived = row.get('_chem_types', {})

        mol = Chem.MolFromSmiles(chuckles)
        if mol is None:
            errors.append(f"{token}: cannot parse CHUCKLES '{chuckles}'")
            continue

        rdDepictor.Compute2DCoords(mol)
        mol.SetProp('symbol', token)
        mol.SetProp('m_abbr', token)
        mol.SetProp('m_name', row.get('name', token))
        mol.SetProp('m_type', row.get('type', 'aa'))
        mol.SetProp('m_subtype', 'natural' if row.get('type', '') == 'aa' else 'cap')

        max_slot = max((s for s, v in leaving.items() if v), default=3)
        lg_vals = [leaving.get(s, None) for s in range(1, max_slot + 1)]
        mol.SetProp('m_Rgroups', ','.join(str(v) if v else 'None' for v in lg_vals))

        if chem_types_derived:
            mol.SetProp('m_chem_types',
                        ','.join(f"{s}:{ct}" for s, ct in sorted(chem_types_derived.items())))

        writer.write(mol)

    writer.close()
    ok = sum(1 for r in rows if r.get('chuckles')) - len(errors)
    return ok, errors


def build_library_from_csv(csv_path, output_sdf, rebuild=False):
    """
    Derive CHUCKLES for all monomers and write to SDF + update CSV in place.

    Thin wrapper around derive_monomers() + write_sdf() for convenience.

    :param csv_path: path to the monomer CSV file.
    :param output_sdf: path for the output SDF.
    :param rebuild: if True, re-derive all 'aa' rows even if CHUCKLES present.
    :returns: (ok_count, errors)
    """
    csv_path = Path(csv_path)
    rows, errors = derive_monomers(csv_path, rebuild=rebuild)
    write_sdf(rows, output_sdf)

    # Write updated CSV back (strip transient _chem_types key)
    for row in rows:
        row.pop('_chem_types', None)
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer_csv = csv.DictWriter(f, fieldnames=_CSV_COLUMNS, extrasaction='ignore')
        writer_csv.writeheader()
        writer_csv.writerows(rows)

    return len(rows) - len(errors), errors


def import_helm_sdf(helm_sdf_path, csv_out_path, peptide_only=True):
    """
    Convert a PistoiaHELM-format monomer SDF into a monomer CSV for this pipeline.

    Reads the HELM SDF and writes a CSV with columns matching _CSV_COLUMNS.
    The `input` column is populated with canonical SMILES extracted from each
    mol block (R-group dummy atoms stripped); `chuckles` and r*_leaving are left
    blank so build_library_from_csv() will derive them via SMARTS pre-activation.

    :param helm_sdf_path: path to the HELM monomer SDF.
    :param csv_out_path: path for the output CSV.
    :param peptide_only: if True (default), skip non-PEPTIDE polymer types.
    :returns: (ok_count, skipped_count)
    """
    helm_sdf_path = Path(helm_sdf_path)
    df = PandasTools.LoadSDF(str(helm_sdf_path))

    rows = []
    skipped = 0
    for _, row in df.iterrows():
        p_type = str(row.get('polymerType', '')).strip()
        if peptide_only and p_type != 'PEPTIDE':
            skipped += 1
            continue

        mol = row.get('ROMol')
        if mol is None:
            skipped += 1
            continue

        # Strip HELM R-group dummy atoms (symbol starts with 'R' + digit)
        emol = Chem.RWMol(mol)
        r_indices = sorted(
            [a.GetIdx() for a in emol.GetAtoms()
             if a.GetSymbol().startswith('R') and len(a.GetSymbol()) > 1],
            reverse=True,
        )
        for idx in r_indices:
            emol.RemoveAtom(idx)
        try:
            Chem.SanitizeMol(emol)
            smiles = Chem.MolToSmiles(emol)
        except Exception:
            skipped += 1
            continue

        symbol   = str(row.get('symbol', '')).strip().replace('-', '_')
        name     = str(row.get('name', symbol)).strip()
        m_type   = str(row.get('monomerType', '')).strip()
        nat      = str(row.get('naturalAnalog', '')).strip()
        csv_type = 'aa' if m_type == 'Backbone' else 'cap'
        synonyms = nat if nat and nat != symbol else ''

        rows.append({
            'token':      symbol,
            'input':      smiles,
            'name':       name,
            'type':       csv_type,
            'synonyms':   synonyms,
            'chuckles':   '',
            'r1_leaving': '',
            'r2_leaving': '',
            'r3_leaving': '',
            'r4_leaving': '',
        })

    csv_out_path = Path(csv_out_path)
    with open(csv_out_path, 'w', newline='', encoding='utf-8') as f:
        writer_csv = csv.DictWriter(f, fieldnames=_CSV_COLUMNS, extrasaction='ignore')
        writer_csv.writeheader()
        writer_csv.writerows(rows)

    return len(rows), skipped


# ── Private helpers ───────────────────────────────────────────────────────

def _find_leaving_atoms(mol, attach_idx, leaving_smiles, exclude=None):
    """Return atom indices to remove from mol for this leaving group.

    :param exclude: set of atom indices already claimed by an earlier slot on the
        same or a different attachment atom.  Used when the same heavy atom appears
        in two slots (e.g. backbone N in R1 and backbone-N mod slot) so each slot
        claims a DIFFERENT H neighbour.
    """
    exclude = exclude or set()
    lg_mol = Chem.MolFromSmiles(leaving_smiles)
    lg_root_num = lg_mol.GetAtomWithIdx(0).GetAtomicNum()
    attach_atom = mol.GetAtomWithIdx(attach_idx)

    for nb in attach_atom.GetNeighbors():
        if nb.GetAtomicNum() != lg_root_num:
            continue
        if nb.GetIdx() in exclude:
            continue  # already claimed by an earlier slot
        if leaving_smiles == '[OH]':
            h_count = sum(1 for h in nb.GetNeighbors() if h.GetAtomicNum() == 1)
            if h_count == 0 and nb.GetTotalNumHs() == 0:
                continue
            to_remove = [nb.GetIdx()]
            to_remove += [h.GetIdx() for h in nb.GetNeighbors()
                          if h.GetAtomicNum() == 1]
            return to_remove
        else:
            return [nb.GetIdx()]
    return []
