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

# Route declared chem_type pairs through their reactions in reactions.yaml.
REACTION_INDEX: dict = {}
for _entry in REACTIONS.values():
    for _pair in _entry.get('reactant_pairs', []):
        _ct_a, _ct_b = _pair
        REACTION_INDEX[(_ct_a, _ct_b)] = _entry
        if _ct_a != _ct_b:
            REACTION_INDEX[(_ct_b, _ct_a)] = _entry

# Functional-group registry
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
# label_only=True preserves the activation scanner's historical non-protecting
# match: other atoms in that match may have their own sites. Supported bonds are
# determined by REACTION_INDEX, independently of whether a type has reactions.
# Ordering determines priority — first match per atom wins in both derived lists.
#
# (chem_type, pre_smarts, pre_lg, infer_smarts, label_only)
_CHEM_TYPE_REGISTRY = [
    # Sulfur / selenium
    ('thiol',             '[SX2H1:1]',                              '[H]',  '[SX2;H0,H1:1]',                          False),
    ('selenol',           '[SeX2H1:1]',                             '[H]',  '[SeX2;H0,H1:1]',                         False),
    # The precursor contains halide; an activated substitution handle replaces
    # that halide. A retained C-halogen bond does not identify the selected site.
    ('alkyl_halide_c',    '[CX4;!H0:1][Cl,Br,I]',                  None,   None,                                    False),
    # Nitrogen nucleophiles — most-specific first
    # aminooxy:  NH2 (pre) → NH1 after dummy
    ('aminooxy',          '[NH2:1][OX2H0]',                         '[H]',  '[NX3;H0,H1,H2:1][OX2H0]',                False),
    # hydrazide: NH1 (pre) → NH0 after dummy; C(=O) guard vs plain hydrazine
    ('hydrazide',         '[NX3H1:1][NX3H2]',                      '[H]',  '[NX3;H0,H1:1]([NX3H2])C(=O)',            False),
    ('amine_primary',     '[NX3;H2:1]',                             '[H]',  '[NX3;H2:1]',                             False),
    ('guanidinium',       '[NX3;H1:1][CX3](=N)',                   '[H]',  '[NX3:1][CX3](=[#7])',                    True),
    ('guanidinium_imine', '[NX2H1:1]=[CX3]([NX3])[NX3]',           '[H]',  '[NX2;H0,H1:1]=[CX3]([NX3])[NX3]',       False),
    # Carboxyl / oxygen
    # Sidechain COOH: dummy on carbonyl C (same convention as backbone R2).
    # LG=[OH] removes the hydroxyl.  infer_smarts is intentionally None —
    # carboxyl vs aldehyde C atoms are structurally identical in CHUCKLES,
    # so we disambiguate via leaving-group metadata in infer_chem_type.
    # aryl_amide_c: C(=O) bonded to aryl (and OH leaving group in free form) —
    # used by aryl-scaffold arms (PhosOxScaffold) where the amide forms with
    # adjacent residue's αN. Must precede `carboxyl` so aryl-attached amide-C
    # routes through aryl_amide_to_backbone_n reaction, not isopeptide.
    ('aryl_amide_c',      '[CX3:1](=O)([OX2H1])[c]',              '[OH]',  '[CX3:1](=O)[c]',                         False),
    ('carboxyl',          '[CX3:1](=O)[OX2H1]',                   '[OH]',  None,                                    False),
    # Default phenols retain the documented label-only type.
    ('hydroxyl_phenolic', '[OX2H1:1][c]',                           '[H]',  '[OX2H1:1][c]',                           True),
    # aryl_phenol_o: aryl-OH (phenol) where the O forms a covalent (ether/amide)
    # bond to adjacent residue. Distinct from `hydroxyl_phenolic` (label_only;
    # used for residue identification but no bond reaction).
    ('aryl_phenol_o',     '[OX2H1:1][c]',                           '[H]',  '[OX2;H0,H1:1][c]',                       False),
    ('hydroxyl',          '[OX2H1:1][CX4]',                        '[H]',  '[OX2H1:1][CX4]',                         False),
    # quat_c_anchor: sp3 carbon with 3 C neighbours and 0 explicit H (quaternary
    # in the bonded form). Used by chondramide-family bridge-anchor monomers
    # (XlQuat). Precedes generic `carbon` so quaternary C with 3 C neighbours is
    # detected as the bridge anchor.
    ('quat_c_anchor',     '[CX4;H1:1]([CX4])([CX4])[CX4]',        '[H]',  '[CX4;H0,H1:1]([C])([C])([C])',           False),
    # aryl_c_anchor: substituted aromatic ring carbon (H=0) — used by hetero-
    # cycle scaffolds (ImzScaffold) where the aryl C is the anchor for chain
    # extension. Specific enough not to over-match (requires H=0, i.e. fully
    # substituted aryl C).
    ('aryl_c_anchor',     '[cX3;H0:1]([c])[c]',                    '[H]',  '[cX3;H0:1]',                             False),
    # sp3_c_anchor: sp3 carbon that forms a direct C-C bond to an aryl anchor.
    # Used by ValAryl-style residues whose alpha-C bonds directly to an
    # imidazole-C2 (CP01557) or similar aryl junction. The chain exit at R2
    # is this sp3 alpha-C — paired with aryl_c_anchor via aryl_c_c_bond.
    ('sp3_c_anchor',      '[CX4;H1:1]([NX3])[c]',                  '[H]',  '[CX4;H0,H1:1]([NX3])[c]',                False),
    # backbone_c_red: sp3 C without C=O, at slot 2 of a reduced-amide residue
    # (CH2-NH style backbone, no carbonyl). Detected by element + slot in
    # infer_chem_type heuristic (SMARTS can't disambiguate from generic methyls).
    ('backbone_c_red',    '[CX4;H2,H3:1][NX3;!$(N-C=O)]',         '[H]',  '[CX4;H1,H2,H3:1][NX3;!$(N-C=O)]',        False),
    # Aromatic / amide N-H (label-only)
    ('aromatic_nh',       '[nH:1]',                                 '[H]',  '[n:1]',                                  True),
    ('amide_nh',          '[NX3;H1:1][CX3]=O',                     '[H]',  '[NX3:1][CX3]=[O,S]',                    True),
    # Phosphate (P(V) electrophile)
    ('phosphate_p',       '[P:1](=O)([OH])[OH]',                  '[OH]',  '[PX4:1](=[OX1])([O])([O])[#0]',          False),
    # Bioorthogonal click
    ('cyclooctyne_c',     '[CX4;!H0:1][C;r]#[C;r]',               '[H]',  '[CX4;!H0:1][C;r]#[C;r]',                False),
    ('alkyne_c',          '[CX4;!H0:1]C#[CH]',                    '[H]',  '[CX4;!H0:1]C#[CH]',                     False),
    ('azide_alpha_c',     '[CX4;!H0:1][N]=[N+]=[N-]',             '[H]',  '[CX4;!H0:1][N]=[N+]=[N-]',              False),
    # terminal_alkene: internal vinyl C — CH1 (pre) → CH0 after dummy
    ('terminal_alkene',   '[CH1:1]=[CH2]',                         '[H]',  '[CH0,CH1:1]=[CH2]',                     False),
    ('tetrazine_c',       '[CX4;!H0:1][c]1[n][n][c][n][n]1',     '[H]',  '[CX4;!H0:1][c]1[n][n][c][n][n]1',       False),
    # tco_c: !r in both — pre_activate never labels in-ring C; no hand-crafted entries
    ('tco_c',             '[CX4;!H0;!r:1][C;r]=[C;r]',           '[H]',  '[CX4;!H0;!r:1][C;r]=[C;r]',             False),
    # Condensation bioorthogonal
    # aldehyde: CH1 (pre, requires H to distinguish from ketone) → CH0 after dummy + chain + O
    ('aldehyde',          '[CX3H1:1](=O)[!#7;!#1;!#0]',           '[H]',  '[CX3;H0,H1:1](=O)[!#7;!#1;!#0]',        False),
    # formamide_c: formyl (N-CHO) — distinct from amide (N-CO-C, no H)
    ('formamide_c',       '[CX3H1:1](=O)[#7X3]',                  '[H]',  '[CX3;H0,H1:1](=O)[#7X3]',               False),
    ('nhs_ester',         '[CX4;!H0:1]C(=O)ON1C(=O)CCC1=O',      '[H]',  '[CX4;!H0:1]C(=O)ON1C(=O)CCC1=O',       False),
    # maleimide_c: ring alkene CH1 (pre) → CH0 after dummy
    # infer_smarts adds ([*]) so the dummy-bearing vinyl C is match[0] not the other vinyl C
    ('maleimide_c',       '[CH1:1]1=[C][C](=[O])[N][C]1=[O]',    '[H]',  '[CH0,CH1:1]([*])1=[C][C](=[O])[N][C]1=[O]', False),
    # thia_michael_c: beta-C of acrylamide; dummy on terminal CH2 (pre) → sp3 CH1 after activation
    ('thia_michael_c',    '[CH2:1]=[CH]C(=O)[NX3]',              '[H]',  '[CX4:1]([*])CC(=O)[NX3]',                   False),
]

# Derived detection lists — do not edit directly.
# Nitrogen context is shared by raw activation and activated-site inference.
# Keeping the candidate-H scan separate preserves established slot ordering.
_INFERENCE_SMARTS = {
    ct: Chem.MolFromSmarts(infer_smarts) if infer_smarts else None
    for ct, _pre, _lg, infer_smarts, _label_only in _CHEM_TYPE_REGISTRY
}
_EXOTIC_SMARTS = [
    (_INFERENCE_SMARTS[ct], ct)
    for ct, _pre, _lg, infer_smarts, label_only in _CHEM_TYPE_REGISTRY
    if not label_only
]

# Every reaction type needs a SMARTS detector or a supported element heuristic.
_HEURISTIC_TYPES = frozenset({
    'backbone_n', 'backbone_c', 'backbone_o', 'backbone_n_mod',
    'amine_secondary',  # nitrogen context counts non-dummy substituents
    'carbon',           # plain sp3 C without carbonyl; heuristic fallback at end of infer_chem_type
    'carboxyl',         # C-attachment COOH; infer_smarts=None, disambiguated via LG in heuristic
    'element_16',       # sulfonyl/sulfonate S; detected by element number heuristic (element_{sym})
    # Chondramide aryl-ether bridge endpoints — trust monomer declaration; reverse-parser
    # needs no new SMARTS pattern because each side is identifiable from element + slot
    # (R3=O / R4=C) plus the residue context (β-OH-Tyr / α-quat-C anchor).
    'aryl_phenol_o',
    'quat_c_anchor',
    # Reduced-amide (CH2-NH) backbone — C-side endpoint that has no carbonyl.
    # Used in polyamine peptidomimetics (CP02627 bis-cyclam-p-xylene).
    'backbone_c_red',
    # Phosphine-oxide scaffold arm: aryl-C(=O) attaching to backbone N via
    # amide. Slot-4 carbonyl C with [OH] leaving group (carboxylic acid in
    # free form); element-heuristic detection in infer_chem_type.
    'aryl_amide_c',
    # Aryl-C anchor: an aromatic C atom on a ring carbon (e.g. imidazole C2/
    # C4/C5 of the CP01557 ImzScaffold) that bonds to a non-amide partner
    # via a direct C-C bond. Three reactions consume this type, paired with
    # backbone_c (aryl→sp3 αC), backbone_c_red (aryl→sp3 CH2), and
    # aryl_amide_c (aryl→aryl-amide-C).
    'aryl_c_anchor',
    # Reduced-N scaffold anchor: N atom on a scaffold (XylBridge) connecting
    # to a polyamine chain via benzylamine C-N bond.
    'reduced_n',
})
_registry_types = {ct for ct, *_ in _CHEM_TYPE_REGISTRY}
_bond_types = {ct for e in REACTIONS.values() for pair in e.get('reactant_pairs', []) for ct in pair}
_missing = _bond_types - _registry_types - _HEURISTIC_TYPES
assert not _missing, f"chem_types in _BOND_TABLE without SMARTS detection: {_missing}"


def nitrogen_chem_type(mol, attach_idx):
    """Recognize N functionality before assigning an amine or backbone role.

    Patterns accept both the free N-H and its dummy-substituted form. An
    unsupported aromatic N remains aromatic; the reaction index decides which
    connections exist. The dummy count is not an amine substitution count.
    """
    atom = mol.GetAtomWithIdx(attach_idx)
    if atom.GetAtomicNum() != 7:
        raise ValueError("Nitrogen chemistry requires a nitrogen atom")
    for kind in ('aromatic_nh', 'aminooxy', 'hydrazide', 'amide_nh',
                 'guanidinium_imine', 'guanidinium'):
        if any(match[0] == attach_idx
               for match in mol.GetSubstructMatches(_INFERENCE_SMARTS[kind])):
            return kind
    if any(bond.GetBondTypeAsDouble() != 1.0 for bond in atom.GetBonds()):
        return 'element_7'
    real_neighbors = sum(nb.GetAtomicNum() not in (0, 1)
                         for nb in atom.GetNeighbors())
    if real_neighbors <= 1:
        return 'amine_primary'
    return 'amine_secondary' if real_neighbors == 2 else 'element_7'


def _declared_site_type(mol, slot):
    if not mol.HasProp('m_chem_types'):
        return None
    for item in mol.GetProp('m_chem_types').split(','):
        key, separator, value = item.partition(':')
        if separator and key.strip() == str(slot):
            return value.strip()
    return None


def infer_chem_type(mol, attach_idx: int, slot: int = None,
                    leaving: str = None) -> str:
    """
    Infer the chemistry type of an attachment atom from the original monomer mol.

    Recognizes nitrogen functionality and validated attachment roles before
    functional-group SMARTS and element fallbacks. Slot roles never turn an
    aromatic or conjugated nitrogen into a primary amine.

    :param mol: original monomer ROMol (with isotope-labelled [n*] dummies).
    :param attach_idx: local atom index of the attachment atom.
    :param slot: 1-based R-group slot override.  When provided, avoids calling
                 _slot_for_attachment (which is ambiguous for atoms bonded to
                 multiple dummies, e.g. backbone N with [1*] and [3*]).
    :param leaving: leaving-group SMILES (e.g. '[OH]', '[H]') from m_Rgroups.
                    Used to disambiguate carboxyl C from aldehyde C — they are
                    structurally identical in CHUCKLES but differ in LG.
    :returns: chem_type string matching reaction library keys.
    """
    from pyPept.attachments import _slot_for_attachment

    atom = mol.GetAtomWithIdx(attach_idx)
    atomic_number = atom.GetAtomicNum()
    if slot is None:
        slot = _slot_for_attachment(mol, attach_idx)
    declared = _declared_site_type(mol, slot)

    if atomic_number == 7:
        kind = nitrogen_chem_type(mol, attach_idx)
        if kind not in ('amine_primary', 'amine_secondary'):
            return kind
        if slot == 1 or declared == 'backbone_n':
            return 'backbone_n'
        if slot == 3 and any(nb.GetAtomicNum() == 0 and nb.GetIsotope() == 1
                             for nb in atom.GetNeighbors()):
            return 'backbone_n_mod'
        return kind

    # These roles describe the intended partner, which is absent from the free
    # monomer. Accept only declarations on an actual saturated carbon site with
    # an H leaving group, never a declaration contradicting its atom/valence.
    if (atomic_number == 6 and not atom.GetIsAromatic()
            and atom.GetHybridization() == Chem.HybridizationType.SP3
            and leaving in (None, '', '[H]')
            and any(nb.GetAtomicNum() == 0 and nb.GetIsotope() == slot
                    for nb in atom.GetNeighbors())):
        if declared == 'backbone_c_red':
            return declared
        if declared == 'sp3_c_anchor' and any(
                nb.GetAtomicNum() == 7 for nb in atom.GetNeighbors()):
            return declared

    has_carbonyl = atomic_number == 6 and any(
        neighbor.GetAtomicNum() == 8
        and mol.GetBondBetweenAtoms(attach_idx, neighbor.GetIdx()).GetBondTypeAsDouble()
        == 2.0
        for neighbor in atom.GetNeighbors()
    )
    if slot == 1 and atomic_number == 8:
        return 'backbone_o'
    if slot == 2 and has_carbonyl:
        return 'backbone_c'

    # Carboxyl and aldehyde sites have the same activated C(=O) graph. Resolve
    # OH leaving groups before SMARTS can misidentify a carboxyl as an aldehyde.
    if has_carbonyl and leaving == '[OH]':
        if any(neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic()
               for neighbor in atom.GetNeighbors()):
            return 'aryl_amide_c'
        return 'carboxyl'

    if (atomic_number == 6 and atom.GetHybridization() == Chem.HybridizationType.SP3
            and leaving in ('[Cl]', '[Br]', '[I]', 'Cl', 'Br', 'I')):
        return 'alkyl_halide_c'

    # SMARTS-based detection (covers thiol, selenol, and all exotic types)
    for pattern, chem_type in _EXOTIC_SMARTS:
        if pattern is None:
            continue
        # An aryl aldehyde has the same activated C(=O)-aryl graph as an
        # aryl carboxyl slot; its H leaving group identifies the aldehyde.
        if chem_type == 'aryl_amide_c' and leaving == '[H]':
            continue
        for match in mol.GetSubstructMatches(pattern):
            if match[0] == attach_idx:
                return chem_type

    # Heuristic fallbacks for C and O
    if atomic_number == 6:
        if has_carbonyl:
            if leaving == '[H]' or (leaving is None and atom.GetTotalNumHs()):
                return 'aldehyde'
            return 'carboxyl'
        return 'carbon'

    if atomic_number == 8:
        return 'hydroxyl'

    return f'element_{atomic_number}'


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


def _target_smirks(smirks, slot_a, slot_b, iso1, iso2):
    """Substitute each placeholder once, even when slot labels coincide."""
    import re

    replacements = {}
    for slot, isotope in ((slot_a, iso1), (slot_b, iso2)):
        replacements.setdefault(str(slot), []).append(isotope)

    def replace(match):
        values = replacements.get(match.group(1))
        return f'[{values.pop(0)}#0]' if values else match.group()

    return re.sub(r'\[(\d+)\*\]', replace, smirks)


def _inherit_residue_ownership(product, reactants):
    """Carry lineage from the actual reactant order after every reaction step.

    Mapped product atoms lose custom properties in RDKit. Its reactant indices
    refer to this step's inputs, not necessarily the original monomer pair.
    Restore ownership before those inputs are replaced by the next product.
    """
    for atom in product.GetAtoms():
        # Unmapped atoms already retain their original owner. Only mapped
        # junction atoms lose it; avoid reassigning the entire growing chain.
        if atom.HasProp('_residue_idx'):
            continue
        if atom.HasProp('react_idx') and atom.HasProp('react_atom_idx'):
            source = reactants[atom.GetIntProp('react_idx')].GetAtomWithIdx(
                atom.GetIntProp('react_atom_idx'))
            if source.HasProp('_residue_idx'):
                atom.SetIntProp('_residue_idx', source.GetIntProp('_residue_idx'))


def run_bond_smirks(frag1, frag2,
                    iso1: int, iso2: int,
                    entry: dict, intramolecular: bool):
    """Execute a targeted reaction and preserve ownership through every step.

    Numbered dummies identify the two sites, independently of the slot labels
    in a reaction definition. The first successful template/fragment orientation
    determines lineage. Alternatives remain internal for symmetric reactions.
    Grouped reactants close a ring within one fragment; later steps consume the
    preceding product. The public result remains a single RDKit molecule.
    """
    current = None
    slot_a, slot_b = entry.get('slot_a'), entry.get('slot_b')
    for step_i, smirks in enumerate(entry['steps']):
        if step_i == 0:
            targets = [smirks]
            if slot_a is not None and slot_b is not None:
                targets = [
                    _target_smirks(smirks, slot_a, slot_b, first, second)
                    for first, second in ((iso1, iso2), (iso2, iso1))
                ]
            if intramolecular:
                targets = [_group_smirks_for_intramol(target) for target in targets]
                orders = [(frag1,)]
            else:
                orders = [(frag1, frag2), (frag2, frag1)]
        else:
            targets, orders = [smirks], [(current,)]

        products = ()
        for target in targets:
            reaction = AllChem.ReactionFromSmarts(target)
            if reaction is None:
                raise ValueError(
                    f"Bad SMIRKS in '{entry['id']}' step {step_i + 1}: {target!r}")
            for reactants in orders:
                products = reaction.RunReactants(reactants)
                if products:
                    break
            if products:
                break
        if not products:
            raise ValueError(
                f"SMIRKS step {step_i + 1} of '{entry['id']}' produced no products. "
                f"SMIRKS: {smirks!r}"
            )

        sanitized = []
        for product in products[0]:
            _inherit_residue_ownership(product, reactants)
            try:
                Chem.SanitizeMol(product)
            except Exception as exc:
                raise ValueError(
                    f"Sanitization failed on step {step_i + 1} product of "
                    f"'{entry['id']}': {exc}"
                ) from exc
            sanitized.append(product)
        if entry.get('take_largest', False) and len(sanitized) > 1:
            current = max(sanitized, key=lambda molecule: molecule.GetNumHeavyAtoms())
        else:
            current = sanitized[0]
            for fragment in sanitized[1:]:
                current = Chem.CombineMols(current, fragment)
    return current
