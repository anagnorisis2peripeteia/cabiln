"""Chemical perception shared by monomer ingestion and numbered attachments.

Rule precedence resolves overlapping chemical descriptions. It never assigns
R numbers. Context claims identify atoms belonging to one reactive handle;
merely inspecting an atom does not consume an attachment on that atom.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
from itertools import groupby
import json
from types import MappingProxyType

from rdkit import Chem
from rdkit.Chem import rdqueries


@dataclass(frozen=True)
class SiteRule:
    name: str
    precedence: int
    raw_smarts: str | None
    leaving: str | None
    site_smarts: str | None
    protect_context: bool = True
    infer: bool = True


# Explicit precedence is independent of declaration and execution order.
SITE_RULES = (
    SiteRule('protected_amine', 0, None, None, '[#7:1]-[CX4]([c])([c])[c]', infer=False),
    SiteRule('thiol', 10, '[SX2H1:1]', '[H]', '[SX2;H1,$(S-[#0]):1]'),
    SiteRule('disulfide', 11, None, None, '[SX2:1]-[SX2]'),
    SiteRule('thioester', 12, None, None, '[SX2:1]-[CX3](=[OX1])'),
    SiteRule('thioether', 13, None, None, '[SX2:1]([#6])[#6]'),
    SiteRule('sulfonyl_s', 14, None, None, '[SX4:1](=O)(=O)[#0]'),
    SiteRule('selenol', 20, '[SeX2H1:1]', '[H]', '[SeX2;H1,$([Se]-[#0]):1]'),
    SiteRule('diselenide', 21, None, None, '[SeX2:1]-[SeX2]'),
    SiteRule('selenoester', 22, None, None, '[SeX2:1]-[CX3](=[OX1])'),
    SiteRule('selenoether', 23, None, None, '[SeX2:1]([#6])[#6]'),
    SiteRule('alkyl_halide_c', 30, '[CX4;!H0:1][Cl,Br,I]', None, None),
    SiteRule('aminooxy', 40, '[NH2:1][OX2H0]', '[H]', '[NX3;+0;!$(N-[#6,#7,#15,#16]):1][OX2H0]'),
    SiteRule('substituted_aminooxy', 41, None, None, '[NX3;+0:1]([#6])[OX2H0]', infer=False),
    SiteRule('hydrazide', 50, '[NX3H1:1][NX3H2]', '[H]', '[NX3;H0,H1:1]([NX3H2])C(=O)'),
    SiteRule('amine_primary', 60, '[NX3;+0;H2:1]', '[H]', '[NX3;+0;H2:1]'),
    SiteRule('guanidinium', 70, '[NX3;H1:1][CX3](=N)', '[H]', '[NX3:1][CX3](=[#7])[#7]', protect_context=False, infer=False),
    SiteRule('guanidinium_imine', 80, '[NX2H1:1]=[CX3]([NX3])[NX3]', '[H]', '[#7:1]=[CX3]([NX3])[NX3]'),
    SiteRule('amidine_nh', 81, None, None, '[NX3:1][CX3]=[#7]', infer=False),
    SiteRule('amidine_imine', 82, None, None, '[#7:1]=[CX3][NX3]', infer=False),
    SiteRule('aryl_amide_c', 90, '[CX3:1](=O)([OX2H1])[c]', '[OH]', '[CX3:1](=O)[c]'),
    SiteRule('carboxyl', 100, '[CX3:1](=O)[OX2H1]', '[OH]', None),
    SiteRule('aryl_phenol_o', 120, '[OX2H1:1][c]', '[H]', '[OX2;H1,$(O-[#0]):1][c]', protect_context=False),
    SiteRule('hydroxyl', 130, '[OX2H1:1][CX4]', '[H]', '[OX2;H1,$(O-[#0]):1][CX4]'),
    SiteRule('phosphate_ester_o', 131, None, None, '[OX2:1]([#6])[PX4](=O)'),
    SiteRule('ester_o', 132, None, None, '[OX2:1]([#6])[CX3](=[OX1])'),
    SiteRule('ether', 133, None, None, '[OX2:1]([#6])[#6]'),
    SiteRule('quat_c_anchor', 140, None, '[H]', '[CX4;H0,H1:1]([C])([C])([C])'),
    SiteRule('aryl_c_anchor', 150, None, '[H]', '[cX3;H0:1]'),
    SiteRule('sp3_c_anchor', 160, None, '[H]', '[CX4;H0,H1:1]([NX3])[c]'),
    SiteRule('backbone_c_red', 170, None, '[H]', '[CX4;H1,H2,H3:1][NX3;!$(N-C=O)]'),
    SiteRule('aromatic_nh', 180, '[nH:1]', '[H]', '[n:1]', protect_context=False, infer=False),
    SiteRule('sulfonamide_nh', 181, '[NX3;H1:1]S(=O)(=O)', '[H]', '[NX3:1]S(=O)(=O)', protect_context=False, infer=False),
    SiteRule('phosphoramide_nh', 182, None, None, '[NX3:1]P(=O)', infer=False),
    SiteRule('urea_nh', 183, None, None, '[NX3:1][CX3](=O)[NX3]', infer=False),
    SiteRule('carbamate_nh', 184, None, None, '[NX3:1][CX3](=O)[OX2]', infer=False),
    SiteRule('amide_nh', 190, '[NX3;H1:1][CX3]=O', '[H]', '[NX3:1][CX3]=[O,S]', protect_context=False, infer=False),
    SiteRule('phosphate_p', 200, '[P:1](=O)([OH])[OH]', '[OH]', '[PX4;+0:1](=[OX1])([O,#0])([O,#0])[#0]'),
    SiteRule('phosphate', 201, None, None, '[PX4:1](=[OX1])([O])([O])[O]'),
    SiteRule('cyclooctyne_c', 210, '[CX4;!H0:1][C;r]#[C;r]', '[H]', '[CX4;!H0:1][C;r]#[C;r]'),
    SiteRule('alkyne_c', 220, '[CX4;!H0:1]C#[CH]', '[H]', '[CX4;!H0:1]C#[CH]'),
    SiteRule('azide_alpha_c', 230, '[CX4;!H0:1][N]=[N+]=[N-]', '[H]', '[CX4;!H0:1][N]=[N+]=[N-]'),
    SiteRule('terminal_alkene', 240, '[CH1:1]=[CH2]', '[H]', '[CH0,CH1:1]=[CH2]'),
    SiteRule('tetrazine_c', 250, '[CX4;!H0:1][c]1[n][n][c][n][n]1', '[H]', '[CX4;!H0:1][c]1[n][n][c][n][n]1'),
    SiteRule('tco_c', 260, '[CX4;!H0;!r:1][C;r]=[C;r]', '[H]', '[CX4;!H0;!r:1][C;r]=[C;r]'),
    SiteRule('aldehyde', 270, '[CX3H1:1](=O)[!#7;!#1;!#0]', '[H]', '[CX3;H0,H1:1](=O)[!#7;!#1;!#0]'),
    SiteRule('formamide_c', 280, '[CX3H1:1](=O)[#7X3]', '[H]', '[CX3;H0,H1:1](=O)[#7X3]'),
    SiteRule('nhs_ester', 290, '[CX4;!H0:1]C(=O)ON1C(=O)CCC1=O', '[H]', '[CX4;!H0:1]C(=O)ON1C(=O)CCC1=O'),
    SiteRule('maleimide_c', 300, '[CH1:1]1=[C][C](=[O])[N][C]1=[O]', '[H]', '[CH0,CH1:1]([*])1=[C][C](=[O])[N][C]1=[O]'),
    SiteRule('thia_michael_c', 310, '[CH2:1]=[CH]C(=O)[NX3]', '[H]', '[CX4:1]([*])CC(=O)[NX3]'),
    SiteRule('amine_secondary', 1000, '[NX3;+0;H1;!$(N-C=[O,S,N]);!$(N-S(=O)):1]', '[H]', None, protect_context=False),
)


@dataclass(frozen=True)
class CompiledRule:
    rule: SiteRule
    raw: object
    site: object
    raw_element: int
    site_element: int


def compile_rules(rules):
    """Reject malformed or duplicate definitions before accepting any monomer."""
    compiled, names = [], set()
    for rule in rules:
        if rule.name in names:
            raise ValueError(f'Duplicate chemistry rule: {rule.name}')
        names.add(rule.name)
        patterns = []
        for text in (rule.raw_smarts, rule.site_smarts):
            pattern = Chem.MolFromSmarts(text) if text else None
            if text and pattern is None:
                raise ValueError(f'Invalid SMARTS for {rule.name}: {text}')
            if pattern is not None and pattern.GetAtomWithIdx(0).GetAtomMapNum() != 1:
                raise ValueError(f'{rule.name}: the attachment atom must be first and mapped :1')
            patterns.append(pattern)
        elements = [pattern.GetAtomWithIdx(0).GetAtomicNum() if pattern is not None else 0
                    for pattern in patterns]
        compiled.append(CompiledRule(rule, *patterns, *elements))
    compiled.sort(key=lambda item: item.rule.precedence)
    return MappingProxyType({item.rule.name: item for item in compiled})


COMPILED_RULES = compile_rules(SITE_RULES)
# Bump the convention when classification semantics outside these rules change.
CHEMISTRY_CONVENTION = 'cabiln-site-chemistry-v4'
CHEMISTRY_FINGERPRINT = sha256(json.dumps(
    [CHEMISTRY_CONVENTION, [asdict(rule) for rule in sorted(SITE_RULES, key=lambda rule: rule.name)]],
    sort_keys=True, separators=(',', ':'),
).encode()).hexdigest()
CHEMISTRY_TYPES = frozenset(rule.name for rule in SITE_RULES) | {
    'backbone_n', 'backbone_c', 'backbone_o', 'backbone_n_mod',
    'carbon', 'element_7', 'element_16', 'lactone_c', 'reduced_n',
    'hydroxyl_phenolic',  # legacy declaration for aryl_phenol_o
}


def canonical_atom_order(mol, roles=None):
    """Canonical traversal positions, retaining the caller's atom indices.

    Atom maps are authoring metadata, not a chemical numbering convention.
    Selected backbone roles distinguish otherwise symmetric molecular arms.
    Removing ordinary explicit H keeps numbering identical to implicit-H input.
    """
    if not roles and not any(atom.GetAtomicNum() == 1 or atom.GetAtomMapNum() for atom in mol.GetAtoms()):
        Chem.MolToSmiles(mol)
        return {index: position for position, index in enumerate(json.loads(mol.GetProp('_smilesAtomOutputOrder')))}
    view = Chem.Mol(mol)
    for atom in view.GetAtoms():
        atom.SetIntProp('_site_source', atom.GetIdx())
        atom.SetAtomMapNum((roles or {}).get(atom.GetIdx(), 0))
    view = Chem.RemoveHs(view)
    Chem.MolToSmiles(view)
    return {
        view.GetAtomWithIdx(index).GetIntProp('_site_source'): position
        for position, index in enumerate(json.loads(view.GetProp('_smilesAtomOutputOrder')))
    }


@dataclass(frozen=True)
class DetectedSite:
    atom: int
    leaving: str | None
    functionality: str
    capacity: int = 1


@dataclass(frozen=True)
class SiteChemistry:
    functionality: str
    role: str
    reaction_type: str


class Perception:
    """Cache perception for one graph, either globally or at selected anchors."""

    def __init__(self, mol, *, targeted=False):
        self.mol = mol
        self.rules = COMPILED_RULES
        self._matches = {}
        self._anchors = {}
        self._at_matches = {}
        self._nitrogens = {}
        self.targeted = targeted
        self.elements = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
        if targeted:
            for atom in mol.GetAtoms():
                atom.SetIntProp('_site_query_index', atom.GetIdx())

    def matches(self, name, *, raw=False):
        key = name, raw
        if key not in self._matches:
            compiled = self.rules[name]
            pattern = compiled.raw if raw else compiled.site
            element = compiled.raw_element if raw else compiled.site_element
            possible = pattern is not None and (element == 0 or element in self.elements)
            # A symmetric group can have several attachment atoms in the same
            # atom set. Preserve each mapping so every anchor is classified.
            self._matches[key] = self.mol.GetSubstructMatches(
                pattern, uniquify=False, maxMatches=0
            ) if possible else ()
        return self._matches[key]

    def matches_at(self, name, index):
        """Test a known anchor in its full context, without enumerating other sites."""
        if not self.targeted:
            if name not in self._anchors:
                self._anchors[name] = {match[0] for match in self.matches(name)}
            return index in self._anchors[name]
        key = name, index
        if key not in self._at_matches:
            compiled = self.rules[name]
            pattern = compiled.site
            possible = pattern is not None and compiled.site_element in (
                0, self.mol.GetAtomWithIdx(index).GetAtomicNum()
            )
            if possible:
                query = Chem.Mol(pattern)
                query.GetAtomWithIdx(0).ExpandQuery(
                    rdqueries.HasIntPropWithValueQueryAtom('_site_query_index', index)
                )
                possible = self.mol.HasSubstructMatch(query)
            self._at_matches[key] = possible
        return self._at_matches[key]

    def nitrogen(self, index):
        if index not in self._nitrogens:
            self._nitrogens[index] = self._nitrogen(index)
        return self._nitrogens[index]

    def _nitrogen(self, index):
        atom = self.mol.GetAtomWithIdx(index)
        if atom.GetAtomicNum() != 7:
            raise ValueError('Nitrogen chemistry requires a nitrogen atom')
        # Conjugation and aromaticity take precedence over substitution count.
        for kind in ('protected_amine', 'aromatic_nh', 'sulfonamide_nh', 'phosphoramide_nh',
                     'urea_nh', 'carbamate_nh', 'hydrazide', 'amide_nh', 'aminooxy', 'substituted_aminooxy',
                     'guanidinium_imine', 'guanidinium', 'amidine_imine', 'amidine_nh'):
            if self.matches_at(kind, index):
                return kind
        if atom.GetFormalCharge() or any(bond.GetBondTypeAsDouble() != 1.0 for bond in atom.GetBonds()):
            return 'element_7'
        neighbors = sum(nb.GetAtomicNum() not in (0, 1) for nb in atom.GetNeighbors())
        return 'amine_primary' if neighbors <= 1 else 'amine_secondary' if neighbors == 2 else 'element_7'

    def raw_sites(self, excluded, order):
        """Resolve reactive handles before numbering, independently of scan order.

        A handle can claim its context (e.g. an NHS ester). Competing rules on
        the same atom need distinct explicit precedence. Equivalent matches of
        one handle use canonical traversal, never input atom order.
        """
        candidates = [
            (compiled.rule, match)
            for compiled in self.rules.values()
            for match in self.matches(compiled.rule.name, raw=True)
            if match[0] not in excluded
            and (not self.mol.GetAtomWithIdx(match[0]).GetIsAromatic()
                 or compiled.rule.name == 'aromatic_nh')
        ]
        candidates.sort(key=lambda item: (item[0].precedence, tuple(order[index] for index in item[1])))
        consumed, sites = set(excluded), []
        for _, group in groupby(candidates, key=lambda item: item[0].precedence):
            group = [(rule, match) for rule, match in group if match[0] not in consumed]
            claims = {}
            for rule, match in group:
                claims.setdefault(match[0], set()).add(rule.name)
            if any(len(names) > 1 for names in claims.values()):
                raise ValueError('Overlapping attachment rules need explicit precedence')
            for rule, match in group:
                index = match[0]
                if index in consumed:
                    continue
                atom = self.mol.GetAtomWithIdx(index)
                kind = self.nitrogen(index) if atom.GetAtomicNum() == 7 else rule.name
                if kind == 'protected_amine' or kind.startswith('element_'):
                    continue
                # Multiple replaceable H atoms share chemistry, not a fake
                # primary/secondary distinction between otherwise equal ports.
                capacity = 2 if rule.name == 'amine_primary' else 1
                sites.append(DetectedSite(index, rule.leaving, kind, capacity))
                consumed.add(index)
                if rule.protect_context:
                    consumed.update(match[1:])
        return sorted(sites, key=lambda site: order[site.atom])

    def classify(self, index, slot=None, leaving=None, declared=None):
        from pyPept.leaving_groups import leaving_family

        leaving = leaving_family(leaving)
        atom = self.mol.GetAtomWithIdx(index)
        number = atom.GetAtomicNum()
        has_carbonyl = number == 6 and any(
            neighbor.GetAtomicNum() == 8
            and self.mol.GetBondBetweenAtoms(index, neighbor.GetIdx()).GetBondTypeAsDouble() == 2
            for neighbor in atom.GetNeighbors()
        )
        if number == 7:
            kind = self.nitrogen(index)
        elif (number == 6 and not atom.GetIsAromatic()
              and atom.GetHybridization() == Chem.HybridizationType.SP3
              and leaving in (None, '', '[H]')
              and any(nb.GetAtomicNum() == 0 and nb.GetIsotope() == slot for nb in atom.GetNeighbors())
              and (declared == 'backbone_c_red' or
                   (declared == 'sp3_c_anchor' and any(nb.GetAtomicNum() == 7 for nb in atom.GetNeighbors())))):
            kind = declared
        elif has_carbonyl and leaving == '[OH]':
            kind = 'aryl_amide_c' if any(nb.GetAtomicNum() == 6 and nb.GetIsAromatic()
                                       for nb in atom.GetNeighbors()) else 'carboxyl'
        elif (number == 6 and atom.GetHybridization() == Chem.HybridizationType.SP3
              and leaving in ('[Cl]', '[Br]', '[I]', 'Cl', 'Br', 'I')):
            kind = 'alkyl_halide_c'
        else:
            selected = None
            for item in self.rules.values():
                rule = item.rule
                if selected is not None and rule.precedence > selected.precedence:
                    break
                if (rule.infer and item.site is not None
                        and item.site_element in (0, number)
                        and not (rule.name == 'aryl_amide_c' and leaving == '[H]')
                        and self.matches_at(rule.name, index)):
                    if selected is not None:
                        raise ValueError('Overlapping attachment types need explicit precedence')
                    selected = rule
            if selected is not None:
                kind = selected.name
            elif number == 6:
                kind = ('aldehyde' if leaving == '[H]' or
                        (leaving is None and atom.GetTotalNumHs()) else 'carboxyl') if has_carbonyl else 'carbon'
            else:
                kind = 'hydroxyl' if number == 8 and (
                    atom.GetTotalNumHs() or any(nb.GetAtomicNum() == 0 for nb in atom.GetNeighbors())
                ) else f'element_{number}'

        role = 'sidechain'
        if number == 7:
            if slot == 1 or declared == 'backbone_n':
                role = 'backbone_n'
            elif slot == 3 and any(nb.GetAtomicNum() == 0 and nb.GetIsotope() == 1
                                   for nb in atom.GetNeighbors()):
                role = 'backbone_n_mod'
        elif number == 8 and slot == 1:
            role = 'backbone_o'
        elif slot == 2 and has_carbonyl:
            role = 'backbone_c'
        reaction_type = kind
        if role != 'sidechain' and (number != 7 or kind in ('amine_primary', 'amine_secondary')):
            reaction_type = role
        return SiteChemistry(kind, role, reaction_type)


def nitrogen_chem_type(mol, attach_idx):
    return Perception(mol).nitrogen(attach_idx)


def infer_chem_type(mol, attach_idx, slot=None, leaving=None):
    """Compatibility interface for callers inspecting one numbered attachment."""
    if slot is None:
        slot = next((nb.GetIsotope() for nb in mol.GetAtomWithIdx(attach_idx).GetNeighbors()
                     if nb.GetAtomicNum() == 0), None)
    declared = {}
    if mol.HasProp('m_chem_types'):
        from pyPept.monomer_store import parse_chem_types
        declared = parse_chem_types(mol.GetProp('m_chem_types'))
    return Perception(mol).classify(attach_idx, slot, leaving, declared.get(slot)).reaction_type
