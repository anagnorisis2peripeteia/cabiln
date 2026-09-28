"""Library-derived ownership candidates for direct-bond decomposition.

This module proposes partitions of the molecule passed to it. Reaction adapters
may pass a precursor proposal with additional atoms; its indices are retained
verbatim here. Callers must project those indices to their original source and
validate emitted notation by independent forward assembly.
"""

from __future__ import annotations

import heapq
import threading
from collections import OrderedDict
from collections.abc import Callable
from dataclasses import dataclass

from rdkit import Chem, rdBase
from rdkit.Chem import rdqueries

from pyPept import monomer_store
from pyPept.attachments import attachment_sites
from pyPept.leaving_groups import restore_leaving_groups
from pyPept.structure import require_supported_stereo

_bit_count = getattr(int, "bit_count", lambda value: bin(value).count("1"))


@dataclass(frozen=True)
class Candidate:
    """One template occurrence, using atom indices in the proposal molecule.

    ``ports`` are (owned anchor, external neighbor, consumed slot). ``anchors``
    includes every template slot, even unused slots, as (owned anchor, slot).
    Different slots on the same atom remain distinct. ``attachment_types`` uses
    the same effective chemistry as assembly and the palette.
    """

    symbol: str
    atoms: frozenset[int]
    ports: tuple[tuple[int, int, int], ...]
    anchors: tuple[tuple[int, int], ...]
    has_backbone: bool
    kind: str
    template_size: int
    attachment_types: tuple[tuple[int, str], ...] = ()


@dataclass(frozen=True)
class RecognitionBudgets:
    """Work limits; frontier/assignment limits apply to admitted partition search."""

    max_patterns: int = 50000
    max_matches_per_pattern: int = 1000
    max_candidates: int = 20000
    max_states: int = 50000
    max_solutions: int = 128
    max_assignments_per_partition: int = 8
    max_frontier_states: int = 2000
    max_covers: int = 32
    allow_unspecified_stereo: bool = True

    def __post_init__(self):
        for name in (
            "max_patterns",
            "max_matches_per_pattern",
            "max_candidates",
            "max_states",
            "max_solutions",
            "max_assignments_per_partition",
            "max_frontier_states",
            "max_covers",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")


@dataclass(frozen=True)
class RecognitionSearchResult:
    """Best observed known-atom partitions, with complete covers preferred.

    Uncovered atoms remain unknown; the caller decides whether their boundaries
    can be represented. ``exhausted`` means a state, examined-solution, assignment,
    frontier, or returned-alternative limit prevented complete enumeration.
    ``match_truncated`` means candidate generation was incomplete. No cutoff or
    first cover proves uniqueness.
    """

    covers: tuple[tuple[Candidate, ...], ...]
    exhausted: bool
    match_truncated: bool
    states: int
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class RejectedRegion:
    """A connected unknown region impossible for every boundary slot assignment.

    The caller must establish that the local chemistry cannot be emitted
    regardless of slot choices. The search independently proves the region is
    sealed by selected ownership groups and cannot gain any recognized atoms
    while those owners remain. Only then may it learn a reusable conflict.
    Atom IDs refer to the proposal molecule, including adapter-added atoms.

    ``requires_complete_recognition`` additionally promises that every nonempty
    unknown subset is impossible. For example, no subset of a carbonyl-free
    region can satisfy an emitter that requires a backbone carbonyl. The search
    can then use failure of the candidate union to cover the entire region.
    """

    atoms: frozenset[int]
    reason: str
    requires_complete_recognition: bool = False


@dataclass(frozen=True)
class _Pattern:
    symbol: str
    query: object
    owners: tuple[int, ...]
    atom_states: tuple[tuple, ...]
    ports: tuple[tuple[int, int, int], ...]
    anchors: tuple[tuple[int, int], ...]
    has_backbone: bool
    kind: str
    template_size: int
    attachment_types: tuple[tuple[int, str], ...]


_pattern_cache: OrderedDict = OrderedDict()
_core_cache: OrderedDict = OrderedDict()
_pattern_lock = threading.RLock()
_MAX_CACHED_QUERY_ATOMS = 50000
_ANCHORS_PROP = "_cabiln_recognition_slots"
_PORT_PROP = "_cabiln_recognition_port"


def _atom_state(atom):
    return (
        atom.GetAtomicNum(),
        atom.GetFormalCharge(),
        atom.GetIsotope(),
        atom.GetIsAromatic(),
        atom.GetNumRadicalElectrons(),
        atom.GetTotalNumHs(),
    )


def _candidate_key(candidate):
    return (
        -len(candidate.atoms),
        not candidate.has_backbone,
        candidate.kind != "aa",
        len(candidate.symbol),
        candidate.symbol,
        candidate.ports,
        candidate.anchors,
        tuple(sorted(candidate.atoms)),
    )


def _snapshot():
    # An external atomic replacement may occur between stat and the loader.
    # Retry rather than cache molecules under a stamp from another file version.
    for _ in range(3):
        path = monomer_store.library_path()
        before = path.stat()
        molecules, _ = monomer_store._load_sdf()
        after = path.stat()
        stamp = (str(path), after.st_mtime_ns, after.st_size)
        if (
            before.st_mtime_ns == after.st_mtime_ns
            and before.st_size == after.st_size
            and path == monomer_store.library_path()
        ):
            return stamp, tuple(molecules)
    raise ValueError("Monomer library changed repeatedly during recognition")


def _compile_patterns(molecules, limit):
    patterns, diagnostics = [], []
    truncated = False
    attempted = 0
    params = Chem.AdjustQueryParameters.NoAdjustments()
    params.adjustDegree = True
    params.adjustDegreeFlags = Chem.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES
    params.makeDummiesQueries = True
    for original in sorted(molecules, key=lambda mol: mol.GetProp("symbol")):
        symbol = original.GetProp("symbol")
        try:
            require_supported_stereo(original)
            sites = attachment_sites(original)
            slots = tuple(site["slot"] for site in sites)
            if len(set(slots)) != len(slots):
                raise ValueError("attachment slot numbers must be unique")
            leaving = {site["slot"]: site["leaving"] for site in sites}
            types = tuple((site["slot"], site["chem_type"]) for site in sites)
            by_slot = dict(types)
            has_backbone = by_slot.get(1) in ("backbone_n", "backbone_o") and (
                by_slot.get(2) in ("backbone_c", "backbone_c_red", "sp3_c_anchor")
            )
            tagged = Chem.Mol(original)
            anchors = {}
            for atom in tagged.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    anchors.setdefault(atom.GetNeighbors()[0].GetIdx(), []).append(
                        atom.GetIsotope()
                    )
            for index, anchor_slots in anchors.items():
                tagged.GetAtomWithIdx(index).SetProp(
                    _ANCHORS_PROP, ",".join(map(str, sorted(anchor_slots)))
                )
            kind = original.GetProp("m_type") if original.HasProp("m_type") else "aa"
            size = original.GetNumHeavyAtoms()
        except ValueError as error:
            diagnostics.append(f"Template {symbol}: {error}")
            continue
        for mask in range(1 << len(slots)):
            if attempted >= limit:
                truncated = True
                break
            attempted += 1
            used = {slot for index, slot in enumerate(slots) if mask & (1 << index)}
            editable = Chem.RWMol(tagged)
            for atom in editable.GetAtoms():
                if atom.GetAtomicNum() == 0 and atom.GetIsotope() in used:
                    atom.SetIntProp(_PORT_PROP, atom.GetIsotope())
                    # Protect used dummies from terminal restoration. A bonded
                    # labeled hydrogen has the same valence contribution here.
                    atom.SetAtomicNum(1)
                    atom.SetIsotope(999)
            try:
                with rdBase.BlockLogs():
                    restored = restore_leaving_groups(editable, leaving)
                ports, restored_anchors = [], []
                for atom in restored.GetAtoms():
                    if atom.HasProp(_PORT_PROP):
                        ports.append(
                            (
                                atom.GetNeighbors()[0].GetIdx(),
                                atom.GetIdx(),
                                atom.GetIntProp(_PORT_PROP),
                            )
                        )
                        atom.SetAtomicNum(0)
                        atom.SetIsotope(0)
                    if atom.HasProp(_ANCHORS_PROP):
                        restored_anchors.extend(
                            (atom.GetIdx(), int(slot))
                            for slot in atom.GetProp(_ANCHORS_PROP).split(",")
                        )
                restored.UpdatePropertyCache(strict=False)
                owners = tuple(
                    atom.GetIdx()
                    for atom in restored.GetAtoms()
                    if atom.GetAtomicNum() != 0
                )
                if not owners:
                    raise ValueError("template state has no owned atoms")
                query = Chem.AdjustQueryProperties(restored, params)
                for index in owners:
                    source_atom = restored.GetAtomWithIdx(index)
                    query_atom = query.GetAtomWithIdx(index)
                    query_atom.ExpandQuery(
                        rdqueries.FormalChargeEqualsQueryAtom(
                            source_atom.GetFormalCharge()
                        )
                    )
                    query_atom.ExpandQuery(
                        rdqueries.IsotopeEqualsQueryAtom(source_atom.GetIsotope())
                    )
                    query_atom.ExpandQuery(
                        rdqueries.NumRadicalElectronsEqualsQueryAtom(
                            source_atom.GetNumRadicalElectrons()
                        )
                    )
                patterns.append(
                    _Pattern(
                        symbol,
                        query,
                        owners,
                        tuple(_atom_state(restored.GetAtomWithIdx(i)) for i in owners),
                        tuple(ports),
                        tuple(restored_anchors),
                        has_backbone,
                        kind,
                        size,
                        types,
                    )
                )
            except ValueError as error:
                diagnostics.append(f"Template {symbol}, slot mask {mask}: {error}")
        if truncated:
            break
    if truncated:
        diagnostics.append(
            f"Pattern limit ({limit}) reached; library states are incomplete"
        )
    if len(diagnostics) > 6:
        diagnostics = diagnostics[:5] + [
            f"{len(diagnostics) - 5} further template compilation diagnostics"
        ]
    return tuple(patterns), truncated, tuple(diagnostics)


def _template_core(molecule):
    """Return a necessary heavy-atom subgraph, without state-specific chemistry.

    Slot restoration changes dummies and hydrogens, but preserves every other
    atom's element and every bond between those atoms. ANY-bond queries also
    tolerate Kekulization/aromaticity changes. Charge, isotope, stereo, degree,
    valence and hydrogen counts deliberately impose no prefilter constraints.
    Exact compiled-state matching still checks those properties afterward.
    """
    core = Chem.RWMol()
    indices = {}
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() > 1:
            indices[atom.GetIdx()] = core.AddAtom(
                rdqueries.AtomNumEqualsQueryAtom(atom.GetAtomicNum())
            )
    if not indices:
        return None  # No necessary heavy atoms: every proposal remains viable.
    any_bond = Chem.BondFromSmarts("~")
    for bond in molecule.GetBonds():
        first, second = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if first in indices and second in indices:
            core.AddBond(indices[first], indices[second], Chem.BondType.UNSPECIFIED)
            core.ReplaceBond(core.GetNumBonds() - 1, any_bond)
    return core.GetMol()


def _patterns(budgets, molecule):
    stamp, molecules = _snapshot()
    with _pattern_lock:
        if stamp not in _core_cache:
            _core_cache[stamp] = tuple(_template_core(item) for item in molecules)
            while len(_core_cache) > 2:
                _core_cache.popitem(last=False)
        _core_cache.move_to_end(stamp)
        viable = tuple(
            index
            for index, query in enumerate(_core_cache[stamp])
            if query is None or molecule.HasSubstructMatch(query, useChirality=False)
        )
        key = stamp, budgets.max_patterns, viable
        if key in _pattern_cache:
            _pattern_cache.move_to_end(key)
            return _pattern_cache[key][0]
        compiled = _compile_patterns(
            tuple(molecules[index] for index in viable), budgets.max_patterns
        )
        query_atoms = sum(pattern.query.GetNumAtoms() for pattern in compiled[0])
        # A small input must not retain every library state. Subset caches also
        # have an aggregate size bound; a larger necessary working set is used
        # for this call but not retained for unrelated later inputs.
        if query_atoms <= _MAX_CACHED_QUERY_ATOMS:
            while _pattern_cache and (
                len(_pattern_cache) >= 2
                or sum(entry[1] for entry in _pattern_cache.values()) + query_atoms
                > _MAX_CACHED_QUERY_ATOMS
            ):
                _pattern_cache.popitem(last=False)
            _pattern_cache[key] = compiled, query_atoms
        return compiled


def _enumerate_candidates(molecule, patterns, budgets):
    candidates, diagnostics = set(), []
    truncated = False
    params = Chem.SubstructMatchParameters()
    params.useChirality = True
    params.specifiedStereoQueryMatchesUnspecified = budgets.allow_unspecified_stereo
    params.uniquify = False  # Symmetric atoms can carry distinct attachment slots.
    params.maxMatches = budgets.max_matches_per_pattern + 1
    source_states = tuple(_atom_state(atom) for atom in molecule.GetAtoms())
    for pattern in patterns:
        matches = molecule.GetSubstructMatches(pattern.query, params)
        if len(matches) > budgets.max_matches_per_pattern:
            truncated = True
            matches = matches[: budgets.max_matches_per_pattern]
        for match in matches:
            # Molecular substructure queries do not constrain neutral charge or
            # isotope zero. Explicit checks prevent false known-piece ownership.
            if any(
                state != source_states[match[index]]
                for index, state in zip(pattern.owners, pattern.atom_states)
            ):
                continue
            if any(
                molecule.GetAtomWithIdx(match[index]).GetChiralTag()
                != Chem.ChiralType.CHI_UNSPECIFIED
                and pattern.query.GetAtomWithIdx(index).GetChiralTag()
                == Chem.ChiralType.CHI_UNSPECIFIED
                for index in pattern.owners
            ):
                continue
            if any(
                molecule.GetBondBetweenAtoms(
                    match[bond.GetBeginAtomIdx()], match[bond.GetEndAtomIdx()]
                ).GetStereo()
                not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
                and bond.GetStereo()
                in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
                for bond in pattern.query.GetBonds()
            ):
                continue
            atoms = frozenset(match[index] for index in pattern.owners)
            ports = tuple(
                sorted(
                    (match[inside], match[outside], slot)
                    for inside, outside, slot in pattern.ports
                )
            )
            if any(outside in atoms for _, outside, _ in ports):
                continue
            candidate = Candidate(
                pattern.symbol,
                atoms,
                ports,
                tuple(sorted((match[index], slot) for index, slot in pattern.anchors)),
                pattern.has_backbone,
                pattern.kind,
                pattern.template_size,
                pattern.attachment_types,
            )
            candidates.add(candidate)
            if len(candidates) > budgets.max_candidates:
                diagnostics.append(
                    f"Candidate limit ({budgets.max_candidates}) reached; matches are incomplete"
                )
                return (
                    tuple(
                        sorted(candidates, key=_candidate_key)[: budgets.max_candidates]
                    ),
                    True,
                    tuple(diagnostics),
                )
    if truncated:
        diagnostics.append(
            f"Match limit ({budgets.max_matches_per_pattern} per pattern) reached; matches are incomplete"
        )
    return tuple(sorted(candidates, key=_candidate_key)), truncated, tuple(diagnostics)


def _cover_key(cover):
    owners = {atom: candidate for candidate in cover for atom in candidate.atoms}
    backbone_edges = 0
    for candidate in cover:
        for inside, outside, slot in candidate.ports:
            other = owners.get(outside)
            if slot == 2 and other is not None and (outside, inside, 1) in other.ports:
                backbone_edges += 1
    return -len(owners), -backbone_edges, tuple(_candidate_key(item) for item in cover)


def _search(molecule, candidates, budgets, accept_cover=None):
    # Fast proposals avoid chemistry callbacks while retaining source ownership.
    # If none assembles acceptably, admitted search must explore different atom
    # partitions before exhausting their interchangeable slot/name assignments.
    if accept_cover is not None:
        return _search_admissible(molecule, candidates, budgets, accept_cover)
    atom_count = molecule.GetNumAtoms()
    all_atoms = (1 << atom_count) - 1
    masks = tuple(
        sum(1 << atom for atom in candidate.atoms) for candidate in candidates
    )
    candidate_atoms = 0
    for mask in masks:
        candidate_atoms |= mask
    coverage_upper_bound = _bit_count(candidate_atoms)
    by_atom = [[] for _ in range(atom_count)]
    port_pairs = tuple(
        frozenset((inside, outside) for inside, outside, _ in candidate.ports)
        for candidate in candidates
    )
    backbone_out = tuple(
        frozenset(
            (inside, outside) for inside, outside, slot in candidate.ports if slot == 2
        )
        for candidate in candidates
    )
    backbone_in = tuple(
        frozenset(
            (outside, inside) for inside, outside, slot in candidate.ports if slot == 1
        )
        for candidate in candidates
    )
    for index, candidate in enumerate(candidates):
        for atom in candidate.atoms:
            by_atom[atom].append(index)

    def compatible(index, chosen, owner):
        for inside, outside in port_pairs[index]:
            other = owner.get(outside)
            if other is not None and (outside, inside) not in port_pairs[other]:
                return False
        for other in chosen:
            for inside, outside in port_pairs[other]:
                if (
                    outside in candidates[index].atoms
                    and (outside, inside) not in port_pairs[index]
                ):
                    return False
        return True

    # Explicit stack avoids recursion limits for large or mostly unknown input.
    stack = [((), 0, 0)]
    results = {}
    best_coverage = -1
    states = 0
    alternatives_truncated = False
    solutions = set()
    solutions_truncated = False

    def record(chosen, used):
        nonlocal best_coverage, alternatives_truncated
        coverage = _bit_count(used)
        if coverage < best_coverage:
            return False
        key = tuple(sorted(chosen))
        cover = tuple(candidates[index] for index in key)
        if coverage > best_coverage:
            best_coverage = coverage
            results.clear()
            alternatives_truncated = False
            solutions.clear()
        results[key] = cover
        if len(results) > budgets.max_covers:
            worst = max(results, key=lambda value: _cover_key(results[value]))
            del results[worst]
            alternatives_truncated = True
        return True

    while stack and states < budgets.max_states:
        chosen, used, excluded = stack.pop()
        states += 1
        if atom_count - _bit_count(excluded) < best_coverage:
            continue
        # Consistent partial progress remains useful when the proposal state
        # budget stops before a terminal partition.
        record(chosen, used)
        owner = {atom: index for index in chosen for atom in candidates[index].atoms}
        available = {
            index
            for index, mask in enumerate(masks)
            if not mask & (used | excluded) and compatible(index, chosen, owner)
        }
        remaining = all_atoms & ~(used | excluded)
        options_by_atom = []
        forced_unknown = 0
        while remaining:
            bit = remaining & -remaining
            atom = bit.bit_length() - 1
            remaining ^= bit
            options = [index for index in by_atom[atom] if index in available]
            if options:
                options_by_atom.append((len(options), atom, options))
            else:
                forced_unknown |= bit
        excluded |= forced_unknown
        if atom_count - _bit_count(excluded) < best_coverage:
            continue
        if not options_by_atom:
            solutions.add(tuple(sorted(chosen)))
            if (
                len(solutions) >= budgets.max_solutions
                and best_coverage == coverage_upper_bound
                and stack
            ):
                solutions_truncated = True
                break
            continue
        possible_out = set().union(
            *(backbone_out[index] for index in available),
            *(backbone_out[index] for index in chosen),
        )
        possible_in = set().union(
            *(backbone_in[index] for index in available),
            *(backbone_in[index] for index in chosen),
        )
        possible_edges = possible_out & possible_in
        if (
            atom_count - _bit_count(excluded) == best_coverage
            and len(results) >= budgets.max_covers
        ):
            worst_backbone_count = min(
                -_cover_key(cover)[1] for cover in results.values()
            )
            # Remaining candidates may overlap each other, so this deliberately
            # overestimates attainable edges. A lower bound would prune unsafely.
            if len(possible_edges) < worst_backbone_count:
                continue
        _, atom, options = min(options_by_atom, key=lambda item: (item[0], item[1]))
        options.sort(
            key=lambda index: (
                -len((backbone_in[index] | backbone_out[index]) & possible_edges),
                index,
            )
        )
        # Unknown is explored after known candidates, and is bounded out once a
        # complete known cover has been found. It permits useful partial covers.
        stack.append((chosen, used, excluded | (1 << atom)))
        for index in reversed(options):
            stack.append((chosen + (index,), used | masks[index], excluded))
    exhausted = bool(stack) or alternatives_truncated or solutions_truncated
    diagnostics = []
    if stack and not solutions_truncated:
        diagnostics.append(f"Search state limit ({budgets.max_states}) reached")
    if solutions_truncated:
        diagnostics.append(f"Examined solution limit ({budgets.max_solutions}) reached")
    if alternatives_truncated:
        diagnostics.append(f"Alternative limit ({budgets.max_covers}) reached")
    return (
        tuple(sorted(results.values(), key=_cover_key)),
        exhausted,
        states,
        tuple(diagnostics),
    )


@dataclass(frozen=True)
class _OwnershipGroup:
    atoms: frozenset[int]
    boundary: frozenset[tuple[int, int]]
    choices: tuple[Candidate, ...]
    mask: int


def _ownership_groups(candidates):
    grouped = {}
    for candidate in candidates:
        boundary = frozenset((a, b) for a, b, _ in candidate.ports)
        grouped.setdefault((candidate.atoms, boundary), []).append(candidate)
    return tuple(
        _OwnershipGroup(
            atoms,
            boundary,
            tuple(choices),
            sum(1 << atom for atom in atoms),
        )
        for (atoms, boundary), choices in grouped.items()
    )


def _assignments(groups, limit):
    """Keep slot/name alternatives separate from the owned-atom partition.

    Best-first enumeration of an assignment grid is bounded, not a proof of a
    global optimum: coupled slot choices may become useful together. Its cutoff
    is reported separately so an invalid assignment cannot monopolize ownership
    search through a combinatorial product of unrelated aliases.
    """
    possible_out = {
        (a, b)
        for group in groups
        for choice in group.choices
        for a, b, slot in choice.ports
        if slot == 2
    }
    possible_in = {
        (b, a)
        for group in groups
        for choice in group.choices
        for a, b, slot in choice.ports
        if slot == 1
    }
    possible_edges = possible_out & possible_in

    def choice_key(choice):
        edges = {
            (a, b) if slot == 2 else (b, a)
            for a, b, slot in choice.ports
            if slot in (1, 2)
        }
        return -len(edges & possible_edges), _candidate_key(choice)

    domains = tuple(tuple(sorted(group.choices, key=choice_key)) for group in groups)

    def cover_at(indices):
        return tuple(
            sorted(
                (domain[index] for domain, index in zip(domains, indices)),
                key=_candidate_key,
            )
        )

    first = (0,) * len(domains)
    first_cover = cover_at(first)
    pending = [(_cover_key(first_cover), first, first_cover)]
    seen = {first}
    results = []
    while pending and len(results) < limit:
        _, indices, cover = heapq.heappop(pending)
        results.append(cover)
        for axis, domain in enumerate(domains):
            if indices[axis] + 1 >= len(domain):
                continue
            neighbor = indices[:axis] + (indices[axis] + 1,) + indices[axis + 1 :]
            if neighbor in seen:
                continue
            seen.add(neighbor)
            other = cover_at(neighbor)
            heapq.heappush(pending, (_cover_key(other), neighbor, other))
    return tuple(results), bool(pending)


def _search_admissible(molecule, candidates, budgets, accept_cover):
    """Search ownership before interpreting each partition's slots and names."""
    groups = _ownership_groups(candidates)
    all_atoms = (1 << molecule.GetNumAtoms()) - 1
    by_atom = [[] for _ in range(molecule.GetNumAtoms())]
    for index, group in enumerate(groups):
        for atom in group.atoms:
            by_atom[atom].append(index)
    pending, queued, visited = [], set(), set()
    serial = 0
    nogoods = []
    frontier_truncated = False

    def violates_conflict(chosen):
        chosen_set = frozenset(chosen)
        return any(conflict <= chosen_set for conflict in nogoods)

    def learn_region(chosen, rejection):
        region = rejection.atoms
        if not region or any(
            atom < 0 or atom >= molecule.GetNumAtoms() for atom in region
        ):
            return False
        owner = {atom: index for index in chosen for atom in groups[index].atoms}
        if any(atom in owner for atom in region):
            return False
        # Check connectedness and that every boundary neighbor is already owned.
        reached, frontier, boundary_owners = set(), [min(region)], set()
        while frontier:
            atom = frontier.pop()
            if atom in reached:
                continue
            reached.add(atom)
            for neighbor in molecule.GetAtomWithIdx(atom).GetNeighbors():
                other = neighbor.GetIdx()
                if other in region:
                    frontier.append(other)
                elif other in owner:
                    boundary_owners.add(owner[other])
                else:
                    return False
        if reached != set(region):
            return False
        boundary_mask = 0
        for index in boundary_owners:
            boundary_mask |= groups[index].mask
        region_mask = sum(1 << atom for atom in region)
        # Use ALL compiled candidate groups, not just those left in this search
        # state. If one could claim any region atom without these boundary owners,
        # rejecting a distant choice might still repair the partition.
        potential = 0
        for group in groups:
            if not group.mask & boundary_mask:
                potential |= group.mask & region_mask
        if rejection.requires_complete_recognition:
            # The caller additionally proves every nonempty unknown subset is
            # impossible (e.g. none can contain a required carbonyl). Leaving
            # even one atom outside the candidate union then proves failure.
            if potential == region_mask:
                return False
        elif potential:
            return False
        conflict = frozenset(boundary_owners)
        if not any(existing <= conflict for existing in nogoods):
            nogoods[:] = [existing for existing in nogoods if not conflict <= existing]
            nogoods.append(conflict)
        return True

    def push(chosen, used, excluded):
        nonlocal serial, frontier_truncated
        chosen = tuple(sorted(chosen))
        if violates_conflict(chosen):
            return
        owner = {atom: index for index in chosen for atom in groups[index].atoms}

        def compatible(index):
            for a, b in groups[index].boundary:
                other = owner.get(b)
                if other is not None and (b, a) not in groups[other].boundary:
                    return False
            return all(
                (b, a) in groups[index].boundary
                for other in chosen
                for a, b in groups[other].boundary
                if b in groups[index].atoms
            )

        available = tuple(
            index
            for index, group in enumerate(groups)
            if not group.mask & (used | excluded)
            and not violates_conflict(chosen + (index,))
            and compatible(index)
        )
        possible = used
        available_mask = 0
        for index in available:
            possible |= groups[index].mask
            available_mask |= 1 << index
        excluded |= all_atoms & ~possible
        key = chosen, excluded
        if key in queued or key in visited:
            return
        serial += 1
        # The union over available, possibly overlapping groups is a sound
        # coverage upper bound. A one-atom unknown alternative precedes a larger
        # residue drop; name and slot aliases cannot create extra queue states.
        entry = (
            -_bit_count(possible),
            -_bit_count(used),
            serial,
            chosen,
            used,
            excluded,
            available_mask,
        )
        if len(pending) >= budgets.max_frontier_states:
            frontier_truncated = True
            worst = max(range(len(pending)), key=lambda index: pending[index][:3])
            if entry[:3] >= pending[worst][:3]:
                return
            removed = pending[worst]
            queued.remove((removed[3], removed[5]))
            pending[worst] = entry
            heapq.heapify(pending)
        else:
            heapq.heappush(pending, entry)
        queued.add(key)

    push((), 0, 0)
    results, decisions = {}, {}
    best_coverage = -1
    states = 0
    admitted_solutions = set()
    assignments_truncated = False
    alternatives_truncated = False
    solutions_truncated = False
    while pending and states < budgets.max_states:
        negative_upper, _, _, chosen, used, excluded, available_mask = heapq.heappop(
            pending
        )
        queued.remove((chosen, excluded))
        visited.add((chosen, excluded))
        states += 1
        if -negative_upper < best_coverage or violates_conflict(chosen):
            continue
        remaining = all_atoms & ~(used | excluded)
        if remaining:
            options = []
            while remaining:
                bit = remaining & -remaining
                remaining ^= bit
                atom = bit.bit_length() - 1
                values = tuple(
                    index for index in by_atom[atom] if available_mask & (1 << index)
                )
                options.append((len(values), atom, values))
            _, atom, values = min(options)
            for index in values:
                push(chosen + (index,), used | groups[index].mask, excluded)
            push(chosen, used, excluded | (1 << atom))
            continue

        coverage = _bit_count(used)
        if coverage < best_coverage:
            continue
        allowance = min(
            budgets.max_assignments_per_partition, budgets.max_states - states
        )
        if not allowance:
            # Retain the unexamined partition so the final status reflects the
            # state budget, rather than silently reporting a completed search.
            pending.append(
                (negative_upper, 0, serial, chosen, used, excluded, available_mask)
            )
            break
        assignments, truncated = _assignments(
            tuple(groups[index] for index in chosen), allowance
        )
        partition_rejected = False
        for cover in assignments:
            states += 1
            if cover not in decisions:
                decisions[cover] = accept_cover(cover)
            decision = decisions[cover]
            if isinstance(decision, RejectedRegion):
                if learn_region(chosen, decision):
                    partition_rejected = True
                    break
                continue
            if not decision:
                continue
            if coverage > best_coverage:
                best_coverage = coverage
                results.clear()
                admitted_solutions.clear()
                alternatives_truncated = False
            results[cover] = cover
            admitted_solutions.add(cover)
            if len(results) > budgets.max_covers:
                worst = max(results, key=_cover_key)
                del results[worst]
                alternatives_truncated = True
        assignments_truncated |= truncated and not partition_rejected
        if (
            len(admitted_solutions) >= budgets.max_solutions
            and pending
            and -pending[0][0] <= best_coverage
        ):
            solutions_truncated = True
            break

    diagnostics = []
    if pending and not solutions_truncated:
        diagnostics.append(f"Search state limit ({budgets.max_states}) reached")
    if solutions_truncated:
        diagnostics.append(
            f"Examined admitted-solution limit ({budgets.max_solutions}) reached"
        )
    if assignments_truncated:
        diagnostics.append(
            f"Slot/name assignment limit ({budgets.max_assignments_per_partition} per ownership partition) reached"
        )
    if alternatives_truncated:
        diagnostics.append(f"Alternative limit ({budgets.max_covers}) reached")
    if frontier_truncated:
        diagnostics.append(
            f"Ownership frontier limit ({budgets.max_frontier_states}) reached"
        )
    exhausted = (
        bool(pending)
        or assignments_truncated
        or alternatives_truncated
        or frontier_truncated
    )
    return (
        tuple(sorted(results.values(), key=_cover_key)),
        exhausted,
        states,
        tuple(diagnostics),
    )


def recognize(
    molecule: Chem.Mol,
    *,
    budgets: RecognitionBudgets | None = None,
    accept_cover: (
        Callable[[tuple[Candidate, ...]], bool | RejectedRegion] | None
    ) = None,
) -> RecognitionSearchResult:
    """Propose complete known covers, or the best observed partial known covers.

    No source atoms are modified and no production notation is emitted here.
    Direct-bond template states are recognized; chemistry that transforms the
    owned core requires a separate reaction adapter and final assembly check.
    ``allow_unspecified_stereo`` permits library-supplied stereochemistry only
    where the proposal leaves it unspecified; callers must disclose inference.
    An optional ``accept_cover`` validates terminal partitions before they may
    establish a coverage/ranking bound. Rejections do not consume the admitted-
    solution limit; predicate calls are memoized per distinct cover. Its errors
    propagate, so the caller should handle expected chemical failures itself.
    ``RejectedRegion`` can supply a slot-independent local failure; the search
    learns a conflict only after proving that the selected boundary owners seal
    the region and exclude every candidate that might recognize any of its atoms,
    or leave an atom unrecognizable under a hereditary local-failure certificate.
    """
    if molecule is None or molecule.GetNumAtoms() == 0:
        raise ValueError("Recognition requires a non-empty molecule")
    require_supported_stereo(molecule)
    budgets = budgets or RecognitionBudgets()
    patterns, compile_truncated, compile_warnings = _patterns(budgets, molecule)
    candidates, match_truncated, match_warnings = _enumerate_candidates(
        molecule, patterns, budgets
    )
    covers, exhausted, states, search_warnings = _search(
        molecule, candidates, budgets, accept_cover
    )
    return RecognitionSearchResult(
        covers,
        exhausted,
        compile_truncated or match_truncated,
        states,
        compile_warnings + match_warnings + search_warnings,
    )
