"""Recognize a molecular structure using the current monomer library.

Recognition assigns source atoms and slots across the complete molecule. Chains,
branches and crosslink notation are presentations of that assignment. Every
returned interpretation is independently reassembled and compared to the input.
"""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass, replace

from rdkit import Chem

from pyPept.monomer_store import library_version
from pyPept.recognition import RecognitionBudgets, RejectedRegion, _prepare_recognition
from pyPept.recognition_notation import (
    Interpretation,
    MonomerAssignment,
    UnknownRegionError,
    emit_interpretation,
    opaque_encodings,
)
from pyPept.recognition_reactions import recognition_variants
from pyPept.structure import compare_structures, require_supported_stereo


class StereoInferenceWarning(UserWarning):
    """The input omitted stereo that the selected library monomers supply."""


@dataclass(frozen=True)
class ConversionResult:
    """One verified library interpretation, with explicit recognition quality.

    ``assignments`` maps emitted occurrences to original input atom indices.
    ``details`` retains the legacy primary-backbone view. ``recognition_status``
    is complete (known templates cover every atom), partial (local unknown pieces)
    or unresolved (opaque component encoding). None asserts a unique synthesis
    history. ``search_complete`` reports whether candidate/proposal search hit a
    budget; it is not a chemical uniqueness claim.
    """

    cabiln: str
    details: list[tuple[str, float, int]]
    inferred_stereo: bool
    synthetic_components: tuple[int, ...]
    warnings: tuple[str, ...]
    assignments: tuple[MonomerAssignment, ...] = ()
    recognition_status: str = "unresolved"
    search_complete: bool = True


_STEREO_INFERENCE_MESSAGE = (
    "Input leaves stereochemistry unspecified; the monomer library supplies "
    "additional stereochemistry."
)


def _library_stamp():
    return library_version()


def _assemble(cabiln, diagnostics=None):
    from pyPept.molecule import Molecule
    from pyPept.sequence import Sequence

    sink = diagnostics.append if diagnostics is not None else lambda _message: None
    try:
        return Molecule(
            Sequence(cabiln, warning_sink=sink), depiction=None
        ).get_molecule(fmt="ROMol")
    except (ValueError, RuntimeError):
        return None


def _reaction_edges(variant, cover):
    """Resolve witnessed reaction endpoints against actual selected library slots."""
    owners = {atom: node for node in cover for atom in node.atoms}
    occupied = {(id(node), slot) for node in cover for _, _, slot in node.ports}
    options = []
    for junction in variant.junctions:
        left = owners.get(junction.endpoint_a)
        right = owners.get(junction.endpoint_b)
        if left is None or right is None:
            return
        left_types, right_types = dict(left.attachment_types), dict(
            right.attachment_types
        )
        left_slots = [
            slot
            for atom, slot in left.anchors
            if atom == junction.endpoint_a
            and left_types.get(slot) == junction.chem_type_a
            and (id(left), slot) not in occupied
        ]
        right_slots = [
            slot
            for atom, slot in right.anchors
            if atom == junction.endpoint_b
            and right_types.get(slot) == junction.chem_type_b
            and (id(right), slot) not in occupied
        ]
        alternatives = [
            (junction.endpoint_a, a, junction.endpoint_b, b)
            for a in left_slots
            for b in right_slots
            if left is not right or a != b
        ]
        if not alternatives:
            return
        options.append(alternatives)
    yield from itertools.product(*options)


def _recognize_component(source, origins, budgets, output_notation):
    proposals = recognition_variants(source)
    limitations = list(proposals.warnings)
    search_complete = not proposals.truncated
    best = None
    best_score = (-1, -1, -1)
    active_variant, problem, checked = None, None, {}
    # Try all reaction proposals before expensive partial-cover exploration.
    # A missing core transformation must not spend the entire budget on direct
    # matches before its supported precursor proposal can be considered.
    for admissible_search in (False, True):
        for variant_index, variant in enumerate(proposals.variants):
            variant_origins = tuple(
                origins[index] if index is not None else None
                for index in variant.atom_origins
            )
            # Retain one proposal's candidates and verification decisions. The
            # usual single-proposal import reuses both phases; switching proposals
            # drops the old state so memory cannot multiply by proposal count.
            if active_variant != variant_index:
                problem = _prepare_recognition(variant.molecule, budgets)
                checked = {}
                active_variant = variant_index

            def accept_cover(cover):
                # Invalid interpretations must not establish the search's coverage
                # threshold: a tiny submatch can otherwise hide a useful unknown
                # residue and prune every correct, slightly smaller known cover.
                nonlocal best, best_score
                accepted = False
                for extra_edges in _reaction_edges(variant, cover):
                    try:
                        interpretation = emit_interpretation(
                            variant.molecule,
                            cover,
                            extra_edges=extra_edges,
                            atom_origins=variant_origins,
                            notation=output_notation,
                        )
                    except UnknownRegionError as error:
                        return RejectedRegion(
                            error.atoms,
                            str(error),
                            requires_complete_recognition=True,
                        )
                    except (ValueError, RuntimeError):
                        continue
                    notation = interpretation.cabiln
                    if notation in checked:
                        accepted |= checked[notation]
                        continue
                    diagnostics = []
                    rebuilt = _assemble(notation, diagnostics)
                    comparison = compare_structures(source, rebuilt)
                    checked[notation] = comparison.compatible
                    if not comparison.compatible:
                        continue
                    accepted = True
                    score = (
                        interpretation.recognized_atoms,
                        -interpretation.unknown_pieces,
                        interpretation.backbone_connections,
                    )
                    if score > best_score:
                        best_score = score
                        best = interpretation, rebuilt, diagnostics
                return accepted

            if admissible_search:
                search = problem.search(accept_cover)
            else:
                # Fast proposals handle fully known structures and reaction products
                # without rebuilding every equivalent template interpretation. This
                # phase cannot rule out an admissible lower-coverage partition.
                search = problem.search()
                for cover in search.covers:
                    accept_cover(cover)
                    if best is not None and not best[0].unknown_pieces:
                        break
            limitations.extend(search.warnings)
            search_complete &= not (search.exhausted or search.match_truncated)
            # This selects one verified interpretation; it does not prove a unique
            # synthesis history or identify which equivalent library alias was used.
            if best is not None and not best[0].unknown_pieces:
                return best, search_complete, tuple(dict.fromkeys(limitations))
            # An unknown-only component can still have uniquely detected free termini
            # even when small cap submatches do not form usable partial covers.
            if not variant.junctions and best is None:
                try:
                    interpretation = emit_interpretation(
                        variant.molecule,
                        (),
                        atom_origins=variant_origins,
                        notation=output_notation,
                    )
                except (ValueError, RuntimeError):
                    continue
                diagnostics = []
                rebuilt = _assemble(interpretation.cabiln, diagnostics)
                if compare_structures(source, rebuilt).compatible:
                    best = interpretation, rebuilt, diagnostics
                    best_score = (0, -interpretation.unknown_pieces, 0)
    if best is not None:
        return best, search_complete, tuple(dict.fromkeys(limitations))

    for notation, _ in opaque_encodings(source):
        diagnostics = []
        rebuilt = _assemble(notation, diagnostics)
        if compare_structures(source, rebuilt).compatible:
            interpretation = Interpretation(notation, (), (), 1, 0, 0)
            diagnostics.append(
                "Could not identify every library residue; preserved the component "
                "as coarse synthetic CABILN. Residue match details are unavailable; "
                "the synthetic slots do not establish peptide residue boundaries."
            )
            return (
                (interpretation, rebuilt, diagnostics),
                search_complete,
                tuple(dict.fromkeys(limitations)),
            )
    raise ValueError(
        "Cannot represent this SMILES component in CABILN without changing its "
        "structure; register the missing monomer or use another input structure"
    )


def convert_smiles(
    smiles: str,
    *,
    notation: str = "percent",
    budgets: RecognitionBudgets | None = None,
) -> ConversionResult:
    """Return a verified, library-driven interpretation of every input component.

    The selection policy prefers recognized atom coverage, then useful local
    unknowns, then R2-to-R1 connections. Canonical source order and a stable
    library order prevent input serialization from selecting different residues.
    Source charge, isotopes, connectivity and specified stereo are mandatory;
    omitted stereo can be supplied by the selected library with a warning.
    """
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("SMILES must be a non-empty string")
    if notation not in ("percent", "bracket"):
        raise ValueError("Notation must be 'percent' or 'bracket'")
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or not molecule.GetNumAtoms():
        raise ValueError(f"Invalid SMILES: {smiles[:80]}")
    require_supported_stereo(molecule)
    stamp = _library_stamp()
    mappings = []
    fragments = Chem.GetMolFrags(molecule, asMols=True, fragsMolAtomMapping=mappings)
    parts, details, diagnostics, synthetic_components, assignments = [], [], [], [], []
    statuses = []
    search_complete = True
    next_bond_id = 1
    next_residue_index = 0
    for component_index, (fragment, input_indices) in enumerate(
        zip(fragments, mappings)
    ):
        ranks = Chem.CanonicalRankAtoms(fragment, breakTies=True)
        order = sorted(range(fragment.GetNumAtoms()), key=lambda index: ranks[index])
        canonical = Chem.RenumberAtoms(fragment, order)
        origins = tuple(input_indices[index] for index in order)
        (interpretation, rebuilt, messages), complete, limits = _recognize_component(
            canonical, origins, budgets, notation
        )
        search_complete &= complete
        messages = list(messages) + list(limits)
        if interpretation.unknown_pieces:
            synthetic_components.append(component_index)
            if interpretation.assignments:
                statuses.append("partial")
                messages.append(
                    f"Preserved {interpretation.unknown_pieces} unmatched region(s) "
                    "as local synthetic CABILN with source attachment sites; "
                    "no exact library match is claimed for those regions."
                )
            else:
                statuses.append("unresolved")
        else:
            statuses.append("complete")
        diagnostics.extend(
            (
                f"Component {component_index + 1}: {message}"
                if len(fragments) > 1
                else message
            )
            for message in messages
        )
        bond_ids = {}

        def renumber(match):
            nonlocal next_bond_id
            key = match.group(0)
            if key not in bond_ids:
                bond_ids[key] = f"!{next_bond_id}"
                next_bond_id += 1
            return bond_ids[key]

        parts.append(re.sub(r"![A-Za-z0-9_]+", renumber, interpretation.cabiln))
        details.extend(interpretation.details)
        assignments.extend(
            replace(item, residue_index=item.residue_index + next_residue_index)
            for item in interpretation.assignments
        )
        from pyPept.sequence import Sequence

        next_residue_index += len(
            Sequence(
                interpretation.cabiln, warning_sink=lambda _message: None
            ).s_monomers
        )
    cabiln = "%".join(parts)
    if len(parts) > 1:
        rebuilt = _assemble(cabiln)
    comparison = compare_structures(molecule, rebuilt)
    if not comparison.compatible:
        raise ValueError("Cannot preserve all SMILES components in one CABILN string")
    if _library_stamp() != stamp:
        raise ValueError(
            "The monomer library changed during conversion; retry the input"
        )
    if not comparison.exact:
        diagnostics.append(_STEREO_INFERENCE_MESSAGE)
    status = (
        "unresolved"
        if "unresolved" in statuses
        else "partial" if "partial" in statuses else "complete"
    )
    return ConversionResult(
        cabiln=cabiln,
        details=details,
        inferred_stereo=not comparison.exact,
        synthetic_components=tuple(synthetic_components),
        warnings=tuple(dict.fromkeys(diagnostics)),
        assignments=tuple(assignments),
        recognition_status=status,
        search_complete=search_complete,
    )


def smiles_to_cabiln_core(smiles: str):
    """Compatibility interface returning notation/details and Python warnings."""
    import warnings

    result = convert_smiles(smiles)
    for message in result.warnings:
        warnings.warn(
            message,
            (
                StereoInferenceWarning
                if message == _STEREO_INFERENCE_MESSAGE
                else UserWarning
            ),
            stacklevel=2,
        )
    return result.cabiln, result.details
