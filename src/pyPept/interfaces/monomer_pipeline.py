"""Convert complete monomers to explicit, numbered attachment templates.

Perception resolves chemical handles before canonical numbering. Backbone or
cap choices that are chemically distinct require author intent. Replacing the
numbered dummies with their leaving groups must recover the exact input graph.
Existing CSV templates retain their sites unless a rebuild is requested.
"""

__credits__ = ["J.B. Brown", "Thomas Fox", "Cameron Beesley"]
__license__ = "MIT"

import csv
import re
import warnings
from collections import deque
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import PandasTools, SDWriter

from pyPept.structure import parse_template_smiles, require_supported_stereo
from pyPept.monomer_store import format_chem_types, monomer_record, parse_chem_types
from pyPept.site_chemistry import Perception, canonical_atom_order, nitrogen_chem_type

ACTIVATION_POLICY = "canonical-sites-v1"

class ActivationError(ValueError):
    """Raised when pre_activate cannot generate a valid CHUCKLES fragment."""


@dataclass
class ActivationResult:
    """Result of pre_activate(). Iterable for backward-compatible tuple unpacking."""
    chuckles: str
    leaving: dict
    chem_types: dict
    policy: str = ACTIVATION_POLICY

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


# Standard amino acid lookup (FASTA / BILN token → SMILES)
# Used by normalize_input() to convert single-letter or token-name inputs.
_AA_SMILES = {
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
    'T': 'N[C@@H]([C@H](O)C)C(=O)O',
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
    'DThr': 'N[C@H]([C@@H](O)C)C(=O)O',
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


# Backbone detection — topology patterns
# Atom-type patterns for the graph-distance backbone search.
# N: trivalent N with at least one H (covers primary amine, Pro ring N,
#    secondary N, but excludes tertiary amines like NMe2 which can't donate H).
# C: carboxyl carbonyl carbon (C(=O)OH).
_BB_N_PAT        = Chem.MolFromSmarts('[NX3;H1,H2,H3;!$(N-C=[O,S,N])]')
# A cyclic lactam N next to the alpha carbon retains its R1 backbone role,
# while its actual functionality remains amide_nh (e.g. pyroglutamate).
_BB_LACTAM_N_PAT = Chem.MolFromSmarts(
    '[NX3;H1;R:1]([CX3]=O)[CX4][CX3](=O)[OX2H1]')
_BB_COOH_PAT     = Chem.MolFromSmarts('[CX3:1](=O)[OX2H1]')
_BB_ALDEHYDE_PAT = Chem.MolFromSmarts('[CX3H1:1](=O)')
_BB_ALCOHOL_PAT  = Chem.MolFromSmarts('[OX2H1:1][CX4]')
_BB_LACTONE_PAT  = Chem.MolFromSmarts('[CX3:1](=O)[OX2;R]')

class BackboneAmbiguity(ActivationError):
    """Several chemically distinct backbone/cap choices need author intent."""

    def __init__(self, choices):
        self.choices = tuple(dict(choice) for choice in choices)
        super().__init__(
            "Ambiguous backbone or cap orientation: choose the intended attachment "
            "path in the preview, or supply backbone_indices."
        )


def _distinct_orientations(mol, choices):
    """Collapse symmetry-equivalent paths while retaining original atom indices."""
    if len(choices) <= 1:
        return choices
    unique = {}
    for choice in choices:
        marked = Chem.Mol(mol)
        for atom in marked.GetAtoms():
            atom.SetAtomMapNum(0)
        for slot, index in choice.items():
            marked.GetAtomWithIdx(index).SetAtomMapNum(slot)
        key = Chem.MolToSmiles(Chem.RemoveHs(marked))
        unique.setdefault(key, choice)
    return [unique[key] for key in sorted(unique)]


def _choose_orientation(mol, choices):
    choices = _distinct_orientations(mol, choices)
    if len(choices) > 1:
        raise BackboneAmbiguity(choices)
    return choices[0] if choices else None


def _backbone_pairs(mol):
    nitrogen = _backbone_n_indices(mol)
    carbon = [match[0] for match in mol.GetSubstructMatches(_BB_COOH_PAT)]
    if not nitrogen and carbon:
        acid_oxygen = {neighbor.GetIdx() for index in carbon
                       for neighbor in mol.GetAtomWithIdx(index).GetNeighbors()
                       if neighbor.GetAtomicNum() == 8}
        nitrogen = [match[0] for match in mol.GetSubstructMatches(_BB_ALCOHOL_PAT)
                    if match[0] not in acid_oxygen]
    if not carbon:
        for pattern in (_BB_ALDEHYDE_PAT, _BB_ALCOHOL_PAT, _BB_LACTONE_PAT):
            carbon = [match[0] for match in mol.GetSubstructMatches(pattern)]
            if carbon:
                break
    pairs = []
    for first in nitrogen:
        for second in carbon:
            if first == second:
                continue
            path = Chem.GetShortestPath(mol, first, second)
            if path:
                pairs.append((len(path) - 1, {1: first, 2: second}))
    return pairs


def _backbone_n_indices(mol):
    perception = Perception(mol)
    return sorted({match[0] for pattern in (_BB_N_PAT, _BB_LACTAM_N_PAT)
                   for match in mol.GetSubstructMatches(pattern)
                   if perception.nitrogen(match[0]) != 'protected_amine'})


def find_backbone_slots(mol):
    """Choose the unique shortest N/O-to-C backbone, modulo graph symmetry."""
    pairs = _backbone_pairs(mol)
    if not pairs:
        return None
    distance = min(item[0] for item in pairs)
    return _choose_orientation(mol, [choice for length, choice in pairs if length == distance])


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
    targets = [atom.GetIdx() for atom in mol.GetAtoms()
               if atom.GetAtomicNum() == 6 and
               any(nb.GetAtomicNum() == 8 for nb in atom.GetNeighbors())]
    choices = []
    for index in n_idxs:
        distances = [len(Chem.GetShortestPath(mol, index, target)) - 1 for target in targets]
        choices.append((min((value for value in distances if value >= 0), default=0), index))
    nearest = min(distance for distance, _ in choices)
    return _choose_orientation(mol, [{1: index} for distance, index in choices if distance == nearest])[1]


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
    n_idxs = _backbone_n_indices(mol)
    s_idxs = [m[0] for m in mol.GetSubstructMatches(_CAP_SULFONYL_PAT)]

    # Multi-arm crosslinker: 3+ alkyl halides → defer to sidechain-only fallback
    if not n_idxs and not c_idxs and not s_idxs:
        alk_multi = [m[0] for m in mol.GetSubstructMatches(_CAP_ALKYL_HALIDE_PAT)]
        if len(alk_multi) >= 3:
            return None, None

    # Free-form electrophilic caps
    if c_idxs and not n_idxs:
        return _choose_orientation(mol, [{2: index} for index in c_idxs]), {2: 'backbone_c'}
    if s_idxs and not n_idxs and not c_idxs:
        return _choose_orientation(mol, [{2: index} for index in s_idxs]), {2: 'element_16'}

    # Reagent-form electrophilic caps (halide LG still attached)
    if not n_idxs and not c_idxs and not s_idxs:
        acyl_h = [m[0] for m in mol.GetSubstructMatches(_CAP_ACYL_HALIDE_PAT)]
        if acyl_h:
            return _choose_orientation(mol, [{2: index} for index in acyl_h]), {2: 'backbone_c'}
        sul_cl = [m[0] for m in mol.GetSubstructMatches(_CAP_SULFONYL_CL_PAT)]
        if sul_cl:
            return _choose_orientation(mol, [{2: index} for index in sul_cl]), {2: 'element_16'}
        alk_h = [m[0] for m in mol.GetSubstructMatches(_CAP_ALKYL_HALIDE_PAT)]
        if alk_h and len(alk_h) < 3:
            return _choose_orientation(mol, [{2: index} for index in alk_h]), {2: 'alkyl_halide_c'}

    # Nucleophilic caps — but prefer electrophilic alkyl halide if present
    # (e.g. acm: Cl-CH2-NH-COCH3 has both amine and alkyl chloride)
    if n_idxs and not c_idxs:
        alk_h = [m[0] for m in mol.GetSubstructMatches(_CAP_ALKYL_HALIDE_PAT)]
        if alk_h and len(alk_h) < 3:
            return _choose_orientation(mol, [{2: index} for index in alk_h]), {2: 'alkyl_halide_c'}
        best = _pick_alpha_n(mol, n_idxs) if len(n_idxs) > 1 else n_idxs[0]
        return {1: best}, {1: 'backbone_n'}
    if not n_idxs and not c_idxs and not s_idxs:
        o_idxs = [m[0] for m in mol.GetSubstructMatches(_CAP_ALCOHOL_PAT)]
        if o_idxs:
            return _choose_orientation(mol, [{1: index} for index in o_idxs]), {1: 'backbone_o'}

    # Carbon nucleophile caps (Grignard / organometallic — last resort)
    for pat in (_CAP_BENZYLIC_CH_PAT, _CAP_TERTIARY_CH_PAT, _CAP_AROMATIC_CH_PAT):
        hits = [m[0] for m in mol.GetSubstructMatches(pat)]
        if hits:
            choice = _choose_orientation(mol, [{1: index} for index in hits])
            kind = 'aryl_c_anchor' if mol.GetAtomWithIdx(choice[1]).GetIsAromatic() else 'carbon'
            return choice, {1: kind}

    return None, None


def find_sidechain_slots(mol, assigned_atoms, start_slot=4, *, atom_order=None):
    """Resolve chemical handles, then number their replaceable groups canonically."""
    sites = Perception(mol).raw_sites(
        assigned_atoms, canonical_atom_order(mol) if atom_order is None else atom_order)
    slots = {}
    for site in sites:
        leaving = site.leaving or infer_leaving_group(mol, site.atom)
        for _ in range(site.capacity):
            slots[start_slot + len(slots)] = (site.atom, leaving, site.functionality)
    return slots


def infer_leaving_group(mol, attachment_idx):
    """Read replaceable groups from bonds after a site has been identified.

    Reagent halides precede acid OH, which precedes H. Formic acid consequently
    loses OH at its acyl port, rather than being mistaken for an aldehyde.
    """
    atom = mol.GetAtomWithIdx(attachment_idx)
    neighbors = list(atom.GetNeighbors())
    if atom.GetAtomicNum() in (6, 16):
        for element, group in ((53, '[I]'), (35, '[Br]'), (17, '[Cl]')):
            if any(nb.GetAtomicNum() == element for nb in neighbors):
                return group
        if any(nb.GetAtomicNum() == 8 and
               mol.GetBondBetweenAtoms(attachment_idx, nb.GetIdx()).GetBondTypeAsDouble() == 2
               for nb in neighbors):
            if any(nb.GetAtomicNum() == 8 and
                   mol.GetBondBetweenAtoms(attachment_idx, nb.GetIdx()).GetBondTypeAsDouble() == 1
                   and (nb.GetTotalNumHs() or any(h.GetAtomicNum() == 1 for h in nb.GetNeighbors()))
                   for nb in neighbors):
                return '[OH]'
    if atom.GetTotalNumHs() or any(nb.GetAtomicNum() == 1 for nb in neighbors):
        return '[H]'
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
    """Every distinct backbone orientation, ordered by distance and graph identity."""
    pairs = _backbone_pairs(mol)
    results = []
    for distance in sorted({length for length, _ in pairs}):
        choices = _distinct_orientations(mol, [choice for length, choice in pairs if length == distance])
        results.extend((choice, distance) for choice in choices)
    return results


_DIST_SUFFIXES = {2: '', 3: '_b', 4: '_g', 5: '_d'}


def pre_activate_all(smiles):
    """Generate CHUCKLES registrations for every backbone orientation.

    Calls pre_activate once per chemically distinct endpoint pairing found by
    find_all_backbone_slots, including longer alternatives such as beta-Asp.
    Symmetry-equivalent orientations produce one entry. Equal-length distinct
    orientations receive numbered suffixes in canonical graph order.

    :param smiles: full monomer SMILES.
    :returns: list of (suffix, ActivationResult) sorted by path length.
              Suffix is '' (alpha), '_b' (beta, path 3), '_g' (gamma, path 4),
              '_d' (delta, path 5), or '_xN' for longer paths.
    :raises ActivationError: if no backbone pairings found.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ActivationError(f"Invalid SMILES: {smiles}")
    try:
        require_supported_stereo(mol)
    except ValueError as exc:
        raise ActivationError(str(exc)) from exc
    mol_h = Chem.AddHs(mol)
    pairings = find_all_backbone_slots(mol_h)
    if not pairings:
        raise ActivationError(f"No backbone pairings found in: {smiles}")
    results, counts = [], {}
    for backbone_dict, dist in pairings:
        suffix = _DIST_SUFFIXES.get(dist, f'_x{dist}')
        counts[dist] = counts.get(dist, 0) + 1
        if counts[dist] > 1:
            suffix += f'_{counts[dist]}'
        result = pre_activate(smiles, backbone_indices=backbone_dict)
        results.append((suffix, result))
    return results


def pre_activate(smiles, slot_overrides=None, leaving_overrides=None,
                 backbone_indices=None):
    """
    Convert a full monomer SMILES to a pre-activated CHUCKLES fragment.

    The shortest supported backbone pair determines R1/R2 unless explicitly
    selected. R3 is reserved for the second backbone N hydrogen. Detected
    sidechain handles receive R4+ in canonical atom order, independently of
    rule execution order. Distinct tied backbones/caps raise BackboneAmbiguity.

    :param smiles: full monomer SMILES (any stereo/notation, all atoms present).
    :param slot_overrides: {slot: attachment_atom_idx} to add/override slots.
    :param leaving_overrides: {slot: leaving_smiles} to override auto-derived
           leaving groups.  Each override is validated against the mol.
    :param backbone_indices: {1: n_idx, 2: c_idx} to select a supported backbone,
           or a single R1/R2 entry for a cap. Indices refer to the parsed input
           SMILES before hydrogens are added. Used by pre_activate_all for β/γ
           variants and to resolve explicit orientation choices.
    :returns: ActivationResult with .chuckles, .leaving, .chem_types, .policy.
              Raises ActivationError on failure.
              Also supports legacy tuple unpacking: chuckles, leaving, chem_types, err = pre_activate(...)
              where err is always None (failures raise instead).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ActivationError(f"Invalid SMILES: {smiles}")
    try:
        require_supported_stereo(mol)
    except ValueError as exc:
        raise ActivationError(str(exc)) from exc
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        raise ActivationError('Pre-activation requires a full monomer without numbered dummies')
    atom_order = canonical_atom_order(mol)
    mol = Chem.AddHs(mol)

    _sidechain_only = False
    if backbone_indices is not None:
        backbone = dict(backbone_indices)
        if (not backbone or set(backbone) - {1, 2}
                or any(not isinstance(index, int) or index < 0 or index >= mol.GetNumAtoms()
                       or mol.GetAtomWithIdx(index).GetAtomicNum() < 2
                       for index in backbone.values())
                or len(set(backbone.values())) != len(backbone)):
            raise ActivationError('Backbone indices must identify distinct heavy atoms for R1/R2')
        if set(backbone) == {1, 2}:
            choices = [choice for _, choice in _backbone_pairs(mol)]
        else:
            try:
                cap, _ = find_cap_slots(mol)
                choices = [cap] if cap else []
            except BackboneAmbiguity as error:
                choices = error.choices
        if backbone not in choices:
            raise ActivationError('Backbone indices do not select a supported backbone or cap orientation')
    else:
        backbone = find_backbone_slots(mol)
    if backbone is None:
        backbone, backbone_chem_types = find_cap_slots(mol)
        if backbone is None:
            # Last resort: sidechain-only molecule (multi-arm crosslinkers)
            _sc_all = find_sidechain_slots(mol, assigned_atoms=set(), start_slot=4, atom_order=atom_order)
            if _sc_all:
                backbone = {s: info[0] for s, info in _sc_all.items()}
                backbone_chem_types = {s: info[2] for s, info in _sc_all.items()}
                _sidechain_only = True
            else:
                raise ActivationError(
                    "Cannot identify backbone (N + COOH) or a single cap terminus. "
                    "Provide pre-filled CHUCKLES in the input column."
                )
    elif set(backbone) != {1, 2}:
        slot, index = next(iter(backbone.items()))
        atom = mol.GetAtomWithIdx(index)
        if slot == 1:
            kind = {7: 'backbone_n', 8: 'backbone_o'}.get(atom.GetAtomicNum(), 'carbon')
            if atom.GetIsAromatic():
                kind = 'aryl_c_anchor'
        elif atom.GetAtomicNum() == 16:
            kind = 'element_16'
        else:
            kind = 'backbone_c' if atom.GetHybridization() == Chem.HybridizationType.SP2 else 'alkyl_halide_c'
        backbone_chem_types = {slot: kind}
    else:
        r2_atom = mol.GetAtomWithIdx(backbone[2])
        r2_element = r2_atom.GetAtomicNum()
        if r2_element == 8:
            r2_ct = 'hydroxyl'
        elif any(
            nb.GetAtomicNum() == 8 and nb.IsInRing()
            for nb in r2_atom.GetNeighbors()
            if mol.GetBondBetweenAtoms(r2_atom.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 1.0
        ):
            r2_ct = 'lactone_c'
        else:
            r2_ct = 'backbone_c'
        r1_atom = mol.GetAtomWithIdx(backbone[1])
        if r1_atom.GetAtomicNum() == 8:
            r1_ct = 'backbone_o'
        else:
            n_kind = nitrogen_chem_type(mol, backbone[1])
            r1_ct = ('backbone_n' if n_kind in ('amine_primary', 'amine_secondary')
                     else n_kind)
        backbone_chem_types = {1: r1_ct, 2: r2_ct}
        # backbone_n_mod: second H on backbone N (N-methylation slot).
        # Pro's ring N has only 1 H and is skipped naturally.
        # Depsipeptide O atoms don't get backbone_n_mod.
        n_atom = mol.GetAtomWithIdx(backbone[1])
        if n_atom.GetAtomicNum() == 7 and sum(1 for nb in n_atom.GetNeighbors() if nb.GetAtomicNum() == 1) >= 2:
            slot = max(backbone.keys()) + 1  # 3 for standard AA
            backbone[slot] = backbone[1]
            backbone_chem_types[slot] = 'backbone_n_mod'

    _bb_excluded = backbone_exclusions(mol, backbone)
    atom_order = canonical_atom_order(mol, {index: slot for slot, index in backbone.items() if slot in (1, 2)})

    if _sidechain_only:
        sidechain = {}
        sidechain_leaving = {s: _sc_all[s][1] for s in _sc_all}
    else:
        # R3 is reserved even when a secondary backbone N has no second H.
        _sc_start = max(max(backbone.keys()) + 1, 4)
        sidechain = find_sidechain_slots(
            mol,
            assigned_atoms=_bb_excluded,
            start_slot=_sc_start,
            atom_order=atom_order,
        )
    slots = {**backbone,
             **{s: info[0] for s, info in sidechain.items()}}
    sidechain_leaving = {**sidechain_leaving} if _sidechain_only else {s: info[1] for s, info in sidechain.items()}

    chem_types = {
        **backbone_chem_types,
        **{s: info[2] for s, info in sidechain.items()},
    }

    if slot_overrides:
        if any(not isinstance(slot, int) or slot <= 0 or not isinstance(index, int)
               or index < 0 or index >= mol.GetNumAtoms()
               or mol.GetAtomWithIdx(index).GetAtomicNum() < 2
               for slot, index in slot_overrides.items()):
            raise ActivationError('Attachment overrides need positive slots and valid heavy atoms')
        slots.update(slot_overrides)
        for slot in slot_overrides:
            sidechain_leaving.pop(slot, None)

    leaving = {}
    if set(leaving_overrides or {}) - set(slots):
        raise ActivationError('Leaving-group overrides refer to missing attachment slots')
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

    return activate_sites(mol, slots, leaving, chem_types,
                          exact_leaving=set(leaving_overrides or {}))


def backbone_exclusions(mol, backbone):
    """Exclude all equally short backbone paths, never one atom-order tie."""
    excluded = set(backbone.values())
    if 1 in backbone and 2 in backbone:
        def distances(start):
            found, pending = {start: 0}, deque([start])
            while pending:
                current = pending.popleft()
                for atom in mol.GetAtomWithIdx(current).GetNeighbors():
                    index = atom.GetIdx()
                    if index not in found:
                        found[index] = found[current] + 1
                        pending.append(index)
            return found

        first, second = distances(backbone[1]), distances(backbone[2])
        length = first.get(backbone[2])
        excluded.update(index for index in first.keys() & second.keys()
                        if first[index] + second[index] == length)
    if 2 in backbone:
        for bond in mol.GetAtomWithIdx(backbone[2]).GetBonds():
            other = bond.GetOtherAtomIdx(backbone[2])
            if mol.GetAtomWithIdx(other).GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 1:
                excluded.add(other)
    return excluded


def activate_sites(mol, slots, leaving, chem_types, *, exact_leaving=()):
    """Build a template from resolved sites and prove exact source restoration.

    The input has explicit H. Slots and leaving groups are resolved before any
    bond changes. This same assembly boundary serves raw and authored inputs.
    """
    from pyPept.leaving_groups import restore_leaving_groups

    source_identity = Chem.MolToSmiles(Chem.RemoveHs(restore_leaving_groups(mol, leaving)))
    leaving = dict(leaving)
    existing = {atom.GetIsotope() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0}
    # Process slots in order so a shared attachment atom uses a different H
    # for each slot, including backbone N with both R1 and a modification site.
    emol = Chem.RWMol(mol)
    atoms_to_remove = []
    _already_claimed = set()
    for slot, attach_idx in sorted(slots.items()):
        if slot in existing:
            continue
        explicit = slot in exact_leaving
        indices = _find_leaving_atoms(emol, attach_idx, leaving[slot],
                                      exclude=_already_claimed, allow_isotopes=not explicit)
        if not indices:
            raise ActivationError(f'R{slot} has no unclaimed {leaving[slot]} leaving group')
        group_atom = emol.GetAtomWithIdx(indices[0])
        if not explicit and group_atom.GetIsotope():
            group = 'OH' if group_atom.GetAtomicNum() == 8 else group_atom.GetSymbol()
            leaving[slot] = f'[{group_atom.GetIsotope()}{group}]'
        atoms_to_remove.extend(indices)
        _already_claimed.update(indices)

    # Add dummies at high indices (before removal to avoid index shifts)
    for slot, attach_idx in sorted(slots.items()):
        if slot in existing:
            continue
        dummy = Chem.Atom(0)
        dummy.SetIsotope(slot)
        new_idx = emol.AddAtom(dummy)
        emol.AddBond(attach_idx, new_idx, Chem.BondType.SINGLE)

    for idx in sorted(set(atoms_to_remove), reverse=True):
        emol.RemoveAtom(idx)

    try:
        mol_final = Chem.RemoveHs(emol.GetMol())
    except (ValueError, RuntimeError) as e:
        raise ActivationError(f"Sanitization failed: {e}") from e

    perception = Perception(mol_final)
    chem_types = {
        atom.GetIsotope(): perception.classify(
            atom.GetNeighbors()[0].GetIdx(), atom.GetIsotope(),
            leaving[atom.GetIsotope()], chem_types.get(atom.GetIsotope())
        ).reaction_type
        for atom in mol_final.GetAtoms() if atom.GetAtomicNum() == 0
    }

    restored = restore_leaving_groups(mol_final, leaving)
    if Chem.MolToSmiles(restored) != source_identity:
        raise ActivationError('Attachment detection changed the source structure; supply explicit sites and leaving groups')

    return ActivationResult(
        chuckles=Chem.MolToSmiles(mol_final),
        leaving=leaving,
        chem_types=chem_types,
    )


# CSV authoring pipeline

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


def _row_leaving(row):
    """Read legacy columns and additional explicit slots without renumbering."""
    return {
        int(match.group(1)): value.strip()
        for column, value in row.items()
        if (match := re.fullmatch(r"r([1-9][0-9]*)_leaving", column))
        and value
        and value.strip() != "None"
    }


def _row_chemistry(row, molecule):
    from pyPept.attachments import attachment_sites

    declared = parse_chem_types(row.get("chem_types", ""))
    if declared:
        return declared
    leaving = _row_leaving(row)
    groups = [leaving.get(slot) for slot in range(1, max(leaving, default=0) + 1)]
    return {
        site["slot"]: site["chem_type"] for site in attachment_sites(molecule, groups)
    }


def derive_monomers(csv_path, rebuild=False):
    """Read CSV and derive activated structures and persistent slot metadata.

    Existing CSV files need no schema migration. Chemistry declarations are
    added when absent; explicit declarations and additional leaving-group
    columns survive unchanged builds. No files are written by this transform.
    """
    csv_path = Path(csv_path)
    with open(csv_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    errors = []
    for row in rows:
        token = row.get("token", "").strip()
        existing_chuckles = row.get("chuckles", "").strip()
        row_type = row.get("type", "").strip()
        try:
            existing = parse_template_smiles(existing_chuckles) if existing_chuckles else None
            if existing is not None:
                chem_types = _row_chemistry(row, existing)
                if rebuild:
                    from pyPept.monomer_migration import reprocess_monomer

                    authored = monomer_record(
                        existing, token, _row_leaving(row), chem_types,
                        name=row.get("name", token), m_type=row_type or "aa",
                        m_subtype=row.get("subtype") or ("natural" if row_type == "aa" else "cap"),
                        strict_metadata=False,
                    )
                    processed = reprocess_monomer(authored).molecule
                    row["chuckles"] = Chem.MolToSmiles(processed)
                    row["activation_policy"] = processed.GetProp("m_activation_policy")
                    groups = processed.GetProp("m_Rgroups").split(",")
                    for slot in set(_row_leaving(row)) | set(range(1, len(groups) + 1)):
                        group = groups[slot - 1] if slot <= len(groups) else "None"
                        row[f"r{slot}_leaving"] = "" if group == "None" else group
                    chem_types = parse_chem_types(processed.GetProp("m_chem_types"))
            else:
                normalized, is_chuckles, norm_err = normalize_input(
                    row.get("input", "").strip()
                )
                if norm_err:
                    raise ActivationError(norm_err)
                if is_chuckles:
                    row["chuckles"] = normalized
                    chem_types = _row_chemistry(row, parse_template_smiles(normalized))
                else:
                    overrides = _row_leaving(row) or None
                    result = pre_activate(normalized, leaving_overrides=overrides)
                    row["chuckles"] = result.chuckles
                    row["activation_policy"] = result.policy
                    for slot in (
                        set(_LG_COLS) | set(_row_leaving(row)) | set(result.leaving)
                    ):
                        row[f"r{slot}_leaving"] = result.leaving.get(slot, "") or ""
                    chem_types = result.chem_types
            # Keep the existing programmatic return field; persist a text field
            # as well so the next unchanged CSV build has the same declarations.
            row["_chem_types"] = chem_types
            row["chem_types"] = format_chem_types(chem_types)
        except (ActivationError, ValueError) as error:
            errors.append(f"{token}: {error}")
            warnings.warn(f"Skipping monomer '{token}': {error}")

    if errors:
        warnings.warn(f"{len(errors)} monomer(s) failed:\n" + "\n".join(errors))
    return rows, errors


def write_sdf(rows, output_sdf):
    """Write validated monomer definitions and return (written_count, errors)."""
    errors, written = [], 0
    with SDWriter(str(output_sdf)) as writer:
        for row in rows:
            token = row.get("token", "").strip()
            chuckles = row.get("chuckles", "").strip()
            if not chuckles:
                continue
            try:
                chem_types = row.get("_chem_types")
                if chem_types is None:
                    chem_types = parse_chem_types(row.get("chem_types", "")) or None
                mol = monomer_record(
                    chuckles,
                    token,
                    _row_leaving(row),
                    chem_types,
                    name=row.get("name", token),
                    m_type=row.get("type", "aa"),
                    m_subtype=row.get("subtype") or ("natural" if row.get("type", "") == "aa" else "cap"),
                    strict_metadata=False,
                    activation_policy=row.get("activation_policy"),
                )
            except ValueError as error:
                errors.append(f"{token}: {error}")
                continue
            writer.write(mol)
            written += 1
    return written, errors


def build_library_from_csv(csv_path, output_sdf, rebuild=False):
    """Build SDF and preserve authoring fields for subsequent CSV rebuilds."""
    csv_path = Path(csv_path)
    rows, errors = derive_monomers(csv_path, rebuild=rebuild)
    written, write_errors = write_sdf(rows, output_sdf)
    # Preserve legacy columns, author metadata, and newly discovered slots.
    fields = list(_CSV_COLUMNS)
    for row in rows:
        row.pop("_chem_types", None)
        fields.extend(key for key in row if key not in fields)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer_csv = csv.DictWriter(f, fieldnames=fields)
        writer_csv.writeheader()
        writer_csv.writerows(rows)
    return written, list(dict.fromkeys(errors + write_errors))


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
            require_supported_stereo(mol)
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


# Private helpers

def _find_leaving_atoms(mol, attach_idx, leaving_smiles, exclude=None, *, allow_isotopes=False):
    """Return atom indices to remove from mol for this leaving group.

    :param exclude: set of atom indices already claimed by an earlier slot on the
        same or a different attachment atom.  Used when the same heavy atom appears
        in two slots (e.g. backbone N in R1 and backbone-N mod slot) so each slot
        claims a DIFFERENT H neighbour.
    """
    exclude = exclude or set()
    lg_mol = Chem.MolFromSmiles(leaving_smiles)
    lg_root_num = lg_mol.GetAtomWithIdx(0).GetAtomicNum()
    isotope = lg_mol.GetAtomWithIdx(0).GetIsotope()
    attach_atom = mol.GetAtomWithIdx(attach_idx)

    for nb in sorted(attach_atom.GetNeighbors(), key=lambda atom: atom.GetIsotope()):
        if nb.GetAtomicNum() != lg_root_num:
            continue
        if nb.GetIdx() in exclude:
            continue  # already claimed by an earlier slot
        if nb.GetIsotope() != isotope and not (allow_isotopes and isotope == 0):
            continue
        if lg_root_num == 8:
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
