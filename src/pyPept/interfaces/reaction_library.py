"""
SMIRKS reaction library for pyPept assembly.

Loaded once at import from data/reactions.yaml.  Provides:
  - REACTIONS: dict of all entries keyed by id
  - REACTION_INDEX: (chem_type_a, chem_type_b) -> reaction entry
  - infer_chem_type(mol, attach_idx) -> str
  - run_bond_smirks(frag1, frag2, iso1, iso2, entry, intramolecular) -> ROMol

Intramolecular ring closure uses RDKit's grouped-reactant SMIRKS syntax:
wrap the bimolecular reactant side in parentheses and call RunReactants with
a single molecule.  E.g. for disulfide ring closure:

    ([997*][S:2].[998*][S:4]) >> [S:2][S:4]   →  RunReactants((assembled_mol,))

The globally-unique isotope labels ensure only the intended atom pair reacts.
"""

import yaml
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

_YAML_PATH = Path(__file__).parent.parent / 'data' / 'reactions.yaml'



def _load_reactions():
    with open(_YAML_PATH, encoding='utf-8') as f:
        entries = yaml.safe_load(f)
    return {e['id']: e for e in entries if e.get('id')}


REACTIONS = _load_reactions()

# ── Reaction routing — built from reactant_pairs in reactions.yaml ────────────
# Each YAML entry with a reactant_pairs field defines which chem_type pairs
# trigger it.  No separate _BOND_TABLE needed — add a reaction to YAML and the
# routing is automatic.
REACTION_INDEX: dict = {}
for _entry in REACTIONS.values():
    for _pair in _entry.get('reactant_pairs', []):
        _ct_a, _ct_b = _pair
        REACTION_INDEX[(_ct_a, _ct_b)] = _entry
        if _ct_a != _ct_b:
            REACTION_INDEX[(_ct_b, _ct_a)] = _entry

# ── Functional-group registry — single source of truth ───────────────────────
# Derives both _EXOTIC_SMARTS (assembly-time infer_chem_type) and
# _SIDECHAIN_RULES (pre_activate monomer labelling).
#
# Two SMARTS columns are needed for types where the dummy atom REPLACES an H on
# the attachment atom — the H-count differs between raw SMILES (pre) and CHUCKLES
# (post).  For all other types pre_smarts == infer_smarts.
#
# pre_smarts  — used in _SIDECHAIN_RULES; run against raw SMILES (no dummies).
#               Must require H ≥ 1 on the attachment atom so there IS an H to
#               replace when placing the dummy.
# infer_smarts — used in _EXOTIC_SMARTS; run against CHUCKLES (dummy present).
#               May accept H = 0 because the dummy already replaced the H.
#
# label_only=True: slot reserved in pre_activate but no bond-table reaction.
# Ordering determines priority — first match per atom wins in both derived lists.
#
# (chem_type, pre_smarts, pre_lg, infer_smarts, label_only)
_CHEM_TYPE_REGISTRY = [
    # ── Sulfur / selenium ────────────────────────────────────────────────────
    ('thiol',             '[SX2H1:1]',                              '[H]',  '[SX2;H0,H1:1]',                          False),
    ('selenol',           '[SeX2H1:1]',                             '[H]',  '[SeX2;H0,H1:1]',                         False),
    # ── Alkyl halide — !H0 covers CH2 (pre) and CH1 after dummy took one H ──
    ('alkyl_halide_c',    '[CX4;!H0:1][Cl,Br,I]',                  '[H]',  '[CX4;!H0:1][Cl,Br,I]',                  False),
    # ── Nitrogen nucleophiles — most-specific first ───────────────────────────
    # aminooxy:  NH2 (pre) → NH1 after dummy
    ('aminooxy',          '[NH2:1][OX2H0]',                         '[H]',  '[NH1:1][OX2H0]',                         False),
    # hydrazide: NH1 (pre) → NH0 after dummy; C(=O) guard vs plain hydrazine
    ('hydrazide',         '[NX3H1:1][NX3H2]',                      '[H]',  '[NX3H0:1]([NX3H2])C(=O)',                False),
    ('amine_primary',     '[NX3;H2:1]',                             '[H]',  '[NX3;H2:1]',                             False),
    ('guanidinium',       '[NX3;H1:1][CX3](=N)',                   '[H]',  '[NX3;H1:1][CX3](=N)',                    True),
    ('guanidinium_imine', '[NX2H1:1]=[CX3]([NX3])[NX3]',           '[H]',  '[NX2;H0,H1:1]=[CX3]([NX3])[NX3]',       False),
    # ── Carboxyl / oxygen ────────────────────────────────────────────────────
    # Sidechain COOH: dummy on carbonyl C (same convention as backbone R2).
    # LG=[OH] removes the hydroxyl.  infer_smarts is intentionally None —
    # carboxyl vs aldehyde C atoms are structurally identical in CHUCKLES,
    # so we disambiguate via leaving-group metadata in infer_chem_type.
    ('carboxyl',          '[CX3:1](=O)[OX2H1]',                   '[OH]',  None,                                    False),
    ('hydroxyl_phenolic', '[OX2H1:1][c]',                           '[H]',  '[OX2H1:1][c]',                           True),
    ('hydroxyl',          '[OX2H1:1][CX4]',                        '[H]',  '[OX2H1:1][CX4]',                         False),
    # ── Aromatic / amide N-H (label-only) ────────────────────────────────────
    ('aromatic_nh',       '[nH:1]',                                 '[H]',  '[nH:1]',                                 True),
    ('amide_nh',          '[NX3;H1:1][CX3]=O',                     '[H]',  '[NX3;H1:1][CX3]=O',                     True),
    # ── Phosphate (P(V) electrophile) ────────────────────────────────────────
    ('phosphate_p',       '[P:1](=O)([OH])[OH]',                  '[OH]',  '[P:1](=O)([OH])[OH]',                   False),
    # ── Bioorthogonal click ───────────────────────────────────────────────────
    ('cyclooctyne_c',     '[CX4;!H0:1][C;r]#[C;r]',               '[H]',  '[CX4;!H0:1][C;r]#[C;r]',                False),
    ('alkyne_c',          '[CX4;!H0:1]C#[CH]',                    '[H]',  '[CX4;!H0:1]C#[CH]',                     False),
    ('azide_alpha_c',     '[CX4;!H0:1][N]=[N+]=[N-]',             '[H]',  '[CX4;!H0:1][N]=[N+]=[N-]',              False),
    # terminal_alkene: internal vinyl C — CH1 (pre) → CH0 after dummy
    ('terminal_alkene',   '[CH1:1]=[CH2]',                         '[H]',  '[CH0,CH1:1]=[CH2]',                     False),
    ('tetrazine_c',       '[CX4;!H0:1][c]1[n][n][c][n][n]1',     '[H]',  '[CX4;!H0:1][c]1[n][n][c][n][n]1',       False),
    # tco_c: !r in both — pre_activate never labels in-ring C; no hand-crafted entries
    ('tco_c',             '[CX4;!H0;!r:1][C;r]=[C;r]',           '[H]',  '[CX4;!H0;!r:1][C;r]=[C;r]',             False),
    # ── Condensation bioorthogonal ────────────────────────────────────────────
    # aldehyde: CH1 (pre, requires H to distinguish from ketone) → CH0 after dummy + chain + O
    ('aldehyde',          '[CX3H1:1](=O)[!#7]',                   '[H]',  '[CX3;H0,H1:1](=O)[!#7]',                False),
    # formamide_c: formyl (N-CHO) — distinct from amide (N-CO-C, no H)
    ('formamide_c',       '[CX3H1:1](=O)[NX3]',                   '[H]',  '[CX3;H0,H1:1](=O)[NX3]',                False),
    ('nhs_ester',         '[CX4;!H0:1]C(=O)ON1C(=O)CCC1=O',      '[H]',  '[CX4;!H0:1]C(=O)ON1C(=O)CCC1=O',       False),
    # maleimide_c: ring alkene CH1 (pre) → CH0 after dummy
    # infer_smarts adds ([*]) so the dummy-bearing vinyl C is match[0] not the other vinyl C
    ('maleimide_c',       '[CH1:1]1=[C][C](=[O])[N][C]1=[O]',    '[H]',  '[CH0,CH1:1]([*])1=[C][C](=[O])[N][C]1=[O]', False),
]

# Derived detection lists — do not edit directly.
# _EXOTIC_SMARTS: infer_chem_type at assembly time (label_only excluded, infer_smarts column).
_EXOTIC_SMARTS = [
    (Chem.MolFromSmarts(infer_smarts) if infer_smarts else None, ct)
    for ct, _pre, _lg, infer_smarts, label_only in _CHEM_TYPE_REGISTRY
    if not label_only
]

# Enforce every _BOND_TABLE chem_type has SMARTS detection or a known element heuristic.
_HEURISTIC_TYPES = frozenset({
    'backbone_n', 'backbone_c', 'backbone_o', 'backbone_n_mod',
    'amine_secondary',  # falls through to amine_primary element heuristic
    'carbon',           # plain sp3 C without carbonyl; heuristic fallback at end of infer_chem_type
    'carboxyl',         # C-attachment COOH; infer_smarts=None, disambiguated via LG in heuristic
})
_registry_types = {ct for ct, *_, lo in _CHEM_TYPE_REGISTRY if not lo}
_bond_types = {ct for e in REACTIONS.values() for pair in e.get('reactant_pairs', []) for ct in pair}
_missing = _bond_types - _registry_types - _HEURISTIC_TYPES
assert not _missing, f"chem_types in _BOND_TABLE without SMARTS detection: {_missing}"


def infer_chem_type(mol, attach_idx: int, slot: int = None,
                    leaving: str = None) -> str:
    """
    Infer the chemistry type of an attachment atom from the original monomer mol.

    First tries SMARTS-based exotic pattern matching, then falls back to
    element + slot heuristics.

    :param mol: original monomer ROMol (with isotope-labelled [n*] dummies).
    :param attach_idx: local atom index of the attachment atom.
    :param slot: 0-based R-group slot override.  When provided, avoids calling
                 _slot_for_attachment (which is ambiguous for atoms bonded to
                 multiple dummies, e.g. backbone N with [1*] and [3*]).
    :param leaving: leaving-group SMILES (e.g. '[OH]', '[H]') from m_Rgroups.
                    Used to disambiguate carboxyl C from aldehyde C — they are
                    structurally identical in CHUCKLES but differ in LG.
    :returns: chem_type string matching reaction library keys.
    """
    from pyPept.sequence import _slot_for_attachment

    atom = mol.GetAtomWithIdx(attach_idx)
    sym = atom.GetAtomicNum()
    if slot is None:
        slot = _slot_for_attachment(mol, attach_idx)  # 0-based

    # ── Backbone slots (1-indexed: slot 1 = R1, slot 2 = R2, slot 3 = R3) ────
    if slot == 1:
        if sym == 7: return 'backbone_n'
        if sym == 8: return 'backbone_o'
    if slot == 2:
        if sym == 6:
            # Only classify as backbone_c (carbonyl) if the carbon has a =O neighbor.
            # Non-carbonyl carbons at slot 2 (e.g. benzyl caps) fall through to
            # the 'carbon' heuristic and use generic bond SMIRKS.
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(attach_idx, nb.GetIdx())
                    if bond and bond.GetBondTypeAsDouble() == 2.0:
                        return 'backbone_c'

    # ── SMARTS-based detection (covers thiol, selenol, and all exotic types) ───
    for patt, ct in _EXOTIC_SMARTS:
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            if match[0] == attach_idx:
                return ct

    # ── Heuristic fallbacks for C, N, O ──────────────────────────────────────
    if sym == 7:
        has_r1_dummy = any(
            nb.GetAtomicNum() == 0 and nb.GetIsotope() == 1
            for nb in atom.GetNeighbors()
        )
        if slot == 3 and has_r1_dummy:
            return 'backbone_n_mod'
        return 'amine_primary'

    if sym == 6:
        has_carbonyl = False
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(attach_idx, nb.GetIdx())
                if bond.GetBondTypeAsDouble() == 2.0:
                    has_carbonyl = True
                    break
        if has_carbonyl:
            if leaving == '[OH]':
                return 'carboxyl'
            if leaving == '[H]':
                return 'aldehyde'
            return 'carboxyl'
        return 'carbon'

    if sym == 8:
        return 'hydroxyl'

    return f'element_{sym}'


def _generic_bond_smirks(iso1: int, sym1: int, iso2: int, sym2: int) -> str:
    """Fallback single-bond SMIRKS using globally unique dummy isotopes."""
    pt = Chem.GetPeriodicTable()
    e1 = pt.GetElementSymbol(sym1)
    e2 = pt.GetElementSymbol(sym2)
    return f'[{iso1}*][{e1}:1].[{iso2}*][{e2}:2] >> [{e1}:1][{e2}:2]'


def _group_smirks_for_intramol(smirks: str) -> str:
    """
    Convert a bimolecular SMIRKS to the grouped (intramolecular) form by
    wrapping the entire reactant side in parentheses.

    '[A].[B] >> [P]'  →  '([A].[B]) >> [P]'

    RDKit treats the grouped form as a single-molecule reaction and produces
    a ring rather than a dimer.
    """
    reactants, _, products = smirks.partition('>>')
    return f'({reactants.strip()}) >> {products.strip()}'


def run_bond_smirks(frag1, frag2,
                    iso1: int, iso2: int,
                    entry: dict, intramolecular: bool):
    """
    Form a bond using the SMIRKS steps in *entry*, then return the product.

    For intermolecular bonds the SMIRKS dummy isotopes are replaced with the
    globally-unique iso1/iso2 so only the correct atom pair is targeted.

    For intramolecular ring closure the bimolecular step-0 SMIRKS is converted
    to the grouped form '([A].[B]) >> [P]' and called as RunReactants((mol,))
    — RDKit's native syntax for ring-forming reactions on a single molecule.
    Subsequent unimolecular steps (e.g. IEDDA retro-[4+2]) run unchanged.

    :param frag1: first fragment (carries iso1 dummy); equals frag2 when intramolecular.
    :param frag2: second fragment (carries iso2 dummy).
    :param iso1: globally-unique isotope of the dummy in frag1.
    :param iso2: globally-unique isotope of the dummy in frag2.
    :param entry: reaction YAML entry dict.
    :param intramolecular: True when both dummies are in the same molecule.
    :returns: product ROMol with the bond formed.
    """
    steps = entry['steps']
    take_largest = entry.get('take_largest', False)
    slot_a_iso = entry.get('slot_a')
    slot_b_iso = entry.get('slot_b')

    current: Chem.ROMol = None

    for step_i, smirks in enumerate(steps):
        if step_i == 0:
            targeted = smirks
            if slot_a_iso is not None and slot_b_iso is not None:
                targeted = targeted.replace(f'[{slot_a_iso}*]', f'[{iso1}*]', 1)
                targeted = targeted.replace(f'[{slot_b_iso}*]', f'[{iso2}*]', 1)

            if intramolecular:
                grouped = _group_smirks_for_intramol(targeted)
                rxn = AllChem.ReactionFromSmarts(grouped)
                if rxn is None:
                    raise ValueError(f"Bad grouped SMIRKS in '{entry['id']}' step 1: {grouped!r}")
                products = rxn.RunReactants((frag1,))
            else:
                rxn = AllChem.ReactionFromSmarts(targeted)
                if rxn is None:
                    raise ValueError(f"Bad SMIRKS in '{entry['id']}' step 1: {targeted!r}")

                products = rxn.RunReactants((frag1, frag2))
                if not products:
                    products = rxn.RunReactants((frag2, frag1))

                # Try with iso1/iso2 swapped — needed when the caller placed
                # isotopes in the reverse of what the SMIRKS template expects.
                if not products and slot_a_iso is not None and slot_b_iso is not None:
                    targeted_swap = smirks.replace(f'[{slot_a_iso}*]', f'[{iso2}*]', 1)
                    targeted_swap = targeted_swap.replace(f'[{slot_b_iso}*]', f'[{iso1}*]', 1)
                    rxn_swap = AllChem.ReactionFromSmarts(targeted_swap)
                    if rxn_swap is not None:
                        products = rxn_swap.RunReactants((frag1, frag2))
                        if not products:
                            products = rxn_swap.RunReactants((frag2, frag1))

        else:
            rxn = AllChem.ReactionFromSmarts(smirks)
            if rxn is None:
                raise ValueError(f"Bad SMIRKS in '{entry['id']}' step {step_i+1}: {smirks!r}")
            products = rxn.RunReactants((current,))

        if not products:
            raise ValueError(
                f"SMIRKS step {step_i + 1} of '{entry['id']}' produced no products. "
                f"SMIRKS: {smirks!r}"
            )

        product_mols = list(products[0])
        sanitized = []
        for m in product_mols:
            try:
                Chem.SanitizeMol(m)
                sanitized.append(m)
            except Exception as exc:
                raise ValueError(
                    f"Sanitization failed on step {step_i+1} product of "
                    f"'{entry['id']}': {exc}"
                ) from exc

        if take_largest and len(sanitized) > 1:
            current = max(sanitized, key=lambda m: m.GetNumHeavyAtoms())
        elif len(sanitized) == 1:
            current = sanitized[0]
        else:
            combined = sanitized[0]
            for frag in sanitized[1:]:
                combined = Chem.CombineMols(combined, frag)
            current = combined

    return current
