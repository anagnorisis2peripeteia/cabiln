"""Explicit reaction hypotheses for whole-molecule monomer recognition.

These are matching proposals, not normalized answers. Every proposed assembly
must still be checked against the untouched source molecule. Only the declared
triazole and terminal-alkene metathesis transformations are inverted here.
Metathesis precursor carbons are virtual: they never acquire source ownership.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator

from rdkit import Chem
from rdkit.Chem import AllChem

from pyPept.interfaces.reaction_library import REACTIONS


class ReactionProposalLimitError(RuntimeError):
    """Reaction alternatives exceeded a budget; absence of a cover is unproven."""


@dataclass(frozen=True)
class ReactionJunction:
    """A product witness connecting two declared attachment chemical types.

    Endpoint indices refer to the proposal molecule. They retain their source
    indices for the supported transformations. Resolve actual R slots from the
    selected library candidates; YAML dummy isotope numbers are not those slots.
    """

    reaction_id: str
    endpoint_a: int
    endpoint_b: int
    chem_type_a: str
    chem_type_b: str
    source_atoms: frozenset[int]
    virtual_atoms: tuple[int, ...] = ()
    byproduct_smiles: str | None = None


@dataclass(frozen=True)
class BondChange:
    atoms: tuple[int, int]
    before: str | None
    after: str | None


@dataclass(frozen=True)
class AtomChange:
    atom: int
    # (element, isotope, charge, explicit H, aromatic, no-implicit-H)
    before: tuple | None
    after: tuple


@dataclass(frozen=True)
class ReactionChange:
    reaction_id: str
    bonds: tuple[BondChange, ...]
    atoms: tuple[AtomChange, ...]


@dataclass(frozen=True)
class RecognitionVariant:
    molecule: Chem.Mol
    atom_origins: tuple[int | None, ...]
    junctions: tuple[ReactionJunction, ...] = ()
    changes: tuple[ReactionChange, ...] = ()


@dataclass(frozen=True)
class ReactionProposals:
    variants: tuple[RecognitionVariant, ...]
    truncated: bool = False
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class _Pattern:
    reaction_id: str
    query: Chem.Mol
    endpoints: tuple[int, int]
    chem_types: tuple[str, str]
    kind: str
    # Triazole query indices: C4, C5, N1, N2, N3.
    ring: tuple[int, ...] = ()
    terminal_h: int = 0
    terminal_no_implicit: bool = False
    byproduct_smiles: str | None = None


@dataclass(frozen=True)
class _Site:
    pattern: _Pattern
    match: tuple[int, ...]
    edited_atoms: frozenset[int]


def _mapped_index(molecule, number):
    return next(
        (a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomMapNum() == number),
        None,
    )


def _attachment_map(reactant):
    anchors = [
        a.GetNeighbors()[0].GetAtomMapNum()
        for a in reactant.GetAtoms()
        if a.GetAtomicNum() == 0 and a.GetIsotope() and a.GetDegree() == 1
    ]
    return anchors[0] if len(anchors) == 1 and anchors[0] else None


def _triazole_pattern(reaction_id, reaction, endpoints, chem_types):
    product = Chem.Mol(reaction.GetProductTemplate(0))
    Chem.FastFindRings(product)
    for ring in product.GetRingInfo().AtomRings():
        carbons = [i for i in ring if product.GetAtomWithIdx(i).GetAtomicNum() == 6]
        nitrogens = [i for i in ring if product.GetAtomWithIdx(i).GetAtomicNum() == 7]
        if len(ring) != 5 or len(carbons) != 2 or len(nitrogens) != 3:
            continue
        c4 = next(
            (i for i in carbons if product.GetBondBetweenAtoms(i, endpoints[0])), None
        )
        n1 = next(
            (i for i in nitrogens if product.GetBondBetweenAtoms(i, endpoints[1])), None
        )
        if c4 is None or n1 is None:
            continue
        c5 = next(i for i in carbons if i != c4)
        n3 = next((i for i in nitrogens if product.GetBondBetweenAtoms(i, c4)), None)
        if n3 is None or n3 == n1 or not product.GetBondBetweenAtoms(c5, n1):
            continue
        n2 = next(i for i in nitrogens if i not in (n1, n3))
        if not all(
            product.GetBondBetweenAtoms(a, b) for a, b in ((n1, n2), (n2, n3), (c4, c5))
        ):
            continue
        reactant = reaction.GetReactantTemplate(0)
        carbon_map = product.GetAtomWithIdx(c5).GetAtomMapNum()
        terminal = _mapped_index(reactant, carbon_map)
        if terminal is None:
            continue
        atom = reactant.GetAtomWithIdx(terminal)
        return _Pattern(
            reaction_id,
            product,
            endpoints,
            chem_types,
            "triazole",
            (c4, c5, n1, n2, n3),
            atom.GetNumExplicitHs(),
            atom.GetNoImplicit(),
        )
    return None


def _metathesis_pattern(reaction_id, reaction, endpoints, chem_types, entry):
    """Accept only the declared two-terminal-carbons -> ethylene rule shape."""
    if not entry.get("take_largest") or reaction.GetNumProductTemplates() != 2:
        return None
    product, byproduct = (Chem.Mol(reaction.GetProductTemplate(i)) for i in (0, 1))
    if any(m.GetNumAtoms() != 2 or m.GetNumBonds() != 1 for m in (product, byproduct)):
        return None
    if any(a.GetAtomicNum() != 6 for m in (product, byproduct) for a in m.GetAtoms()):
        return None
    if any(
        m.GetBondWithIdx(0).GetBondType() != Chem.BondType.DOUBLE
        for m in (product, byproduct)
    ):
        return None
    lost_maps = set()
    for reactant in (reaction.GetReactantTemplate(i) for i in (0, 1)):
        if reactant.GetNumAtoms() != 3:
            return None
        anchor = _mapped_index(reactant, _attachment_map(reactant))
        if anchor is None:
            return None
        terminals = [
            a
            for a in reactant.GetAtoms()
            if a.GetAtomicNum() == 6 and a.GetIdx() != anchor
        ]
        if len(terminals) != 1 or terminals[0].GetNumExplicitHs() != 2:
            return None
        bond = reactant.GetBondBetweenAtoms(anchor, terminals[0].GetIdx())
        if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
            return None
        lost_maps.add(terminals[0].GetAtomMapNum())
    if lost_maps != {a.GetAtomMapNum() for a in byproduct.GetAtoms()} or 0 in lost_maps:
        return None
    return _Pattern(
        reaction_id,
        product,
        endpoints,
        chem_types,
        "metathesis",
        byproduct_smiles="C=C",
    )


def _patterns():
    patterns = []
    for reaction_id, entry in REACTIONS.items():
        pairs = [tuple(p) for p in entry.get("reactant_pairs", ())]
        supported = [
            p
            for p in pairs
            if p
            in {
                ("alkyne_c", "azide_alpha_c"),
                ("cyclooctyne_c", "azide_alpha_c"),
                ("terminal_alkene", "terminal_alkene"),
            }
        ]
        steps = entry.get("steps", ())
        if not supported or len(steps) != 1:
            continue
        reaction = AllChem.ReactionFromSmarts(steps[0])
        if (
            reaction is None
            or reaction.GetNumReactantTemplates() != 2
            or not reaction.GetNumProductTemplates()
        ):
            continue
        product = reaction.GetProductTemplate(0)
        endpoints = tuple(
            _mapped_index(product, _attachment_map(reaction.GetReactantTemplate(i)))
            for i in (0, 1)
        )
        if None in endpoints:
            continue
        for chem_types in supported:
            if chem_types == ("terminal_alkene", "terminal_alkene"):
                pattern = _metathesis_pattern(
                    reaction_id, reaction, endpoints, chem_types, entry
                )
            else:
                pattern = _triazole_pattern(
                    reaction_id, reaction, endpoints, chem_types
                )
            if pattern is not None:
                patterns.append(pattern)
    return patterns


def _sites(source, max_matches):
    sites = []
    for pattern in _patterns():
        matches = source.GetSubstructMatches(
            pattern.query, uniquify=True, maxMatches=max_matches + 1
        )
        if len(matches) > max_matches:
            raise ReactionProposalLimitError(
                "Reaction product matching exceeded max_matches."
            )
        for match in matches:
            if pattern.kind == "triazole" and pattern.terminal_no_implicit:
                terminal = source.GetAtomWithIdx(match[pattern.ring[1]])
                if terminal.GetTotalNumHs(includeNeighbors=True) != pattern.terminal_h:
                    continue
            edited = pattern.ring if pattern.kind == "triazole" else pattern.endpoints
            sites.append(_Site(pattern, match, frozenset(match[i] for i in edited)))
            if len(sites) > max_matches:
                raise ReactionProposalLimitError(
                    "Reaction product matching exceeded max_matches."
                )
    return sites


def _atom_state(atom):
    return (
        atom.GetAtomicNum(),
        atom.GetIsotope(),
        atom.GetFormalCharge(),
        atom.GetNumExplicitHs(),
        atom.GetIsAromatic(),
        atom.GetNoImplicit(),
    )


def _bond_states(molecule):
    return {
        tuple(
            sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        ): f"{bond.GetBondType()}:{bond.GetStereo()}:{tuple(bond.GetStereoAtoms())}"
        for bond in molecule.GetBonds()
    }


def _apply(source, selected):
    molecule = Chem.RWMol(source)
    origins = list(range(source.GetNumAtoms()))
    junctions, changes = [], []
    for site in selected:
        pattern, match = site.pattern, site.match
        before_atoms = [_atom_state(a) for a in molecule.GetAtoms()]
        before_bonds = _bond_states(molecule)
        a, b = (match[i] for i in pattern.endpoints)
        virtual = ()
        if pattern.kind == "triazole":
            c4, c5, n1, n2, n3 = (match[i] for i in pattern.ring)
            molecule.RemoveBond(c5, n1)
            molecule.RemoveBond(c4, n3)
            for x, y, kind in (
                (c4, c5, Chem.BondType.TRIPLE),
                (n1, n2, Chem.BondType.DOUBLE),
                (n2, n3, Chem.BondType.DOUBLE),
            ):
                bond = molecule.GetBondBetweenAtoms(x, y)
                bond.SetBondType(kind)
                bond.SetIsAromatic(False)
            for index, charge in ((c4, 0), (c5, 0), (n1, 0), (n2, 1), (n3, -1)):
                atom = molecule.GetAtomWithIdx(index)
                atom.SetIsAromatic(False)
                atom.SetFormalCharge(charge)
                atom.SetNumExplicitHs(pattern.terminal_h if index == c5 else 0)
                atom.SetNoImplicit(
                    pattern.terminal_no_implicit
                    if index == c5
                    else index in (n1, n2, n3)
                )
        else:
            molecule.RemoveBond(a, b)
            added = []
            for anchor in (a, b):
                carbon = Chem.Atom(6)
                carbon.SetNumExplicitHs(2)
                carbon.SetNoImplicit(True)
                index = molecule.AddAtom(carbon)
                origins.append(None)
                added.append(index)
                molecule.AddBond(anchor, index, Chem.BondType.DOUBLE)
            virtual = tuple(added)
        Chem.SanitizeMol(molecule)
        if pattern.chem_types[0] == "cyclooctyne_c" and not all(
            molecule.GetAtomWithIdx(match[i]).IsInRing() for i in pattern.ring[:2]
        ):
            raise ValueError("The restored alkyne is not cyclic.")
        after_bonds = _bond_states(molecule)
        bond_changes = tuple(
            BondChange(pair, before_bonds.get(pair), after_bonds.get(pair))
            for pair in sorted(before_bonds.keys() | after_bonds.keys())
            if before_bonds.get(pair) != after_bonds.get(pair)
        )
        atom_changes = tuple(
            AtomChange(
                i, before_atoms[i] if i < len(before_atoms) else None, _atom_state(atom)
            )
            for i, atom in enumerate(molecule.GetAtoms())
            if i >= len(before_atoms) or before_atoms[i] != _atom_state(atom)
        )
        junctions.append(
            ReactionJunction(
                pattern.reaction_id,
                a,
                b,
                *pattern.chem_types,
                frozenset(match),
                virtual,
                pattern.byproduct_smiles,
            )
        )
        changes.append(ReactionChange(pattern.reaction_id, bond_changes, atom_changes))
    return RecognitionVariant(
        molecule.GetMol(), tuple(origins), tuple(junctions), tuple(changes)
    )


def _iter_variants(
    source: Chem.Mol,
    *,
    max_variants: int = 128,
    max_matches: int = 512,
) -> Iterator[RecognitionVariant]:
    """Yield bounded combinations, with the untouched proposal first."""
    if source is None:
        raise ValueError("A source molecule is required.")
    if max_variants < 1 or max_matches < 1:
        raise ValueError("Reaction proposal limits must be positive.")
    yield RecognitionVariant(Chem.Mol(source), tuple(range(source.GetNumAtoms())))
    sites = _sites(source, max_matches)
    stack = [(0, (), frozenset())]
    count = 1
    states = 0
    while stack:
        states += 1
        if states > max_variants * (len(sites) + 1) * 4:
            raise ReactionProposalLimitError(
                "Reaction alternative search exceeded its state budget."
            )
        index, selected, edited = stack.pop()
        if index == len(sites):
            if not selected:
                continue
            try:
                variant = _apply(source, selected)
            except (ValueError, RuntimeError):
                # A motif match can be incompatible with surrounding valence.
                # Reject that proposal; the original source remains untouched.
                continue
            if count >= max_variants:
                raise ReactionProposalLimitError(
                    "Reaction proposals exceeded max_variants."
                )
            count += 1
            yield variant
            continue
        site = sites[index]
        stack.append((index + 1, selected, edited))
        if not edited & site.edited_atoms:
            stack.append((index + 1, selected + (site,), edited | site.edited_atoms))


def recognition_variants(
    source: Chem.Mol,
    *,
    max_proposals: int = 128,
    max_matches: int = 512,
) -> ReactionProposals:
    """Return original first, then compatible combinations of reaction proposals.

    All original atoms retain their indices, element and isotope. Added atoms
    have ``atom_origins=None`` and a reaction/byproduct witness. Multiple sites
    of the same reaction are considered; overlapping edits are alternatives.
    Exhaustion is request-local metadata, retaining proposals already generated.

    Unsupported transformations (for example IEDDA and NHS ester loss) receive
    no inverse proposal. A product motif alone never proves reaction history.
    """
    variants = []
    try:
        variants.extend(
            _iter_variants(source, max_variants=max_proposals, max_matches=max_matches)
        )
    except ReactionProposalLimitError as error:
        return ReactionProposals(tuple(variants), True, (str(error),))
    return ReactionProposals(tuple(variants))
