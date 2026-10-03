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
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from pyPept.site_chemistry import CHEMISTRY_TYPES, infer_chem_type, nitrogen_chem_type

_YAML_PATH = Path(__file__).parent.parent / 'data' / 'reactions.yaml'


@lru_cache(maxsize=512)
def _compiled_reaction(smirks):
    reaction = AllChem.ReactionFromSmarts(smirks)
    if reaction is None:
        raise ValueError(f'Invalid reaction SMARTS: {smirks}')
    _, errors = reaction.Validate(silent=True)
    if errors:
        raise ValueError(f'Invalid reaction mapping: {smirks}')
    reaction.Initialize()
    return reaction


def compile_reactions(entries):
    """Validate routes and resolve explicit aliases, independent of file order."""
    definitions, reactions, index = {}, {}, {}
    for entry in entries:
        name = entry.get('id')
        if not name or name in definitions:
            raise ValueError(f'Missing or duplicate reaction id: {name!r}')
        definitions[name] = entry

    def resolve(name, visiting=()):
        if name in reactions:
            return reactions[name]
        if name not in definitions or name in visiting:
            raise ValueError(f'Missing or cyclic reaction alias: {name}')
        entry = definitions[name]
        if 'alias_of' in entry:
            if set(entry) - {'id', 'alias_of', 'description', 'notes'}:
                raise ValueError(f'Reaction alias {name} cannot redefine chemistry')
            target = resolve(entry['alias_of'], (*visiting, name))
            entry = {**target, 'id': name, 'alias_of': target.get('alias_of', target['id']),
                     'reactant_pairs': [], 'description': entry.get('description', target.get('description', ''))}
        else:
            if not entry.get('steps'):
                raise ValueError(f'Reaction {name} has no steps')
            for step in entry['steps']:
                _compiled_reaction(step)
        reactions[name] = entry
        return entry

    for name in definitions:
        resolve(name)
    for name, entry in reactions.items():
        for pair in entry.get('reactant_pairs', []):
            if len(pair) != 2 or set(pair) - CHEMISTRY_TYPES:
                raise ValueError(f'Unknown attachment chemistry in {name}: {pair}')
            left, right = pair
            existing = index.get((left, right))
            if existing is not None:
                raise ValueError(f'Duplicate reaction route {left}/{right}: {existing["id"]}, {name}')
            index[left, right] = entry
            index[right, left] = entry
    return reactions, index

with _YAML_PATH.open(encoding='utf-8') as _source:
    REACTIONS, REACTION_INDEX = compile_reactions(yaml.safe_load(_source))
REACTION_FINGERPRINT = sha256(json.dumps(
    REACTIONS, sort_keys=True, separators=(',', ':')
).encode()).digest()


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


def _inherit_residue_ownership(product, reactants, connection_id=None):
    """Carry lineage from the actual reactant order after every reaction step.

    Mapped product atoms lose custom properties in RDKit. Its reactant indices
    refer to this step's inputs, not necessarily the original monomer pair.
    Restore ownership before those inputs are replaced by the next product.
    """
    previous_bonds = [
        (side, bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetIntProp('_connection_idx'))
        for side, molecule in enumerate(reactants) for bond in molecule.GetBonds()
        if bond.HasProp('_connection_idx')
    ]
    endpoints = {(side, index) for side, left, right, _ in previous_bonds for index in (left, right)}
    locations, mapped = {}, set()
    for atom in product.GetAtoms():
        if endpoints and atom.HasProp('react_idx') and atom.HasProp('react_atom_idx'):
            origin = atom.GetIntProp('react_idx'), atom.GetIntProp('react_atom_idx')
            if origin in endpoints:
                locations[origin] = atom.GetIdx()
        # Unmapped atoms already retain their original owner. Only mapped
        # junction atoms lose it; avoid reassigning the entire growing chain.
        if atom.HasProp('_residue_idx') and atom.HasProp('_template_atom'):
            continue
        if connection_id is not None and atom.HasProp("old_mapno") and not atom.HasProp('_residue_idx'):
            mapped.add(atom.GetIdx())
        if atom.HasProp('react_idx') and atom.HasProp('react_atom_idx'):
            source = reactants[atom.GetIntProp('react_idx')].GetAtomWithIdx(
                atom.GetIntProp('react_atom_idx'))
            for key in ('_residue_idx', '_template_atom'):
                if not atom.HasProp(key) and source.HasProp(key):
                    atom.SetIntProp(key, source.GetIntProp(key))

    # Ring closures can rebuild spectator bonds too. Restore every surviving
    # tagged bond through this step's atom correspondence, then tag new bonds.
    for side, left, right, label in previous_bonds:
        if (side, left) in locations and (side, right) in locations:
            bond = product.GetBondBetweenAtoms(locations[side, left], locations[side, right])
            if bond is not None:
                bond.SetIntProp('_connection_idx', label)
    for index in mapped:
        for bond in product.GetAtomWithIdx(index).GetBonds():
            other = bond.GetOtherAtomIdx(index)
            if other not in mapped or other < index:
                continue
            left, right = bond.GetBeginAtom(), bond.GetEndAtom()
            previous = None
            if (
                left.HasProp("react_idx")
                and right.HasProp("react_idx")
                and left.GetIntProp("react_idx") == right.GetIntProp("react_idx")
            ):
                previous = reactants[left.GetIntProp("react_idx")].GetBondBetweenAtoms(
                    left.GetIntProp("react_atom_idx"),
                    right.GetIntProp("react_atom_idx"),
                )
            if previous is None:
                bond.SetIntProp("_connection_idx", connection_id)
            elif previous.HasProp("_connection_idx"):
                bond.SetIntProp(
                    "_connection_idx", previous.GetIntProp("_connection_idx")
                )


def run_bond_smirks(
    frag1,
    frag2,
    iso1: int,
    iso2: int,
    entry: dict,
    intramolecular: bool,
    *,
    connection_id=None,
):
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
            reaction = _compiled_reaction(target)
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
            _inherit_residue_ownership(product, reactants, connection_id)
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
