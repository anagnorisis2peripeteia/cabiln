"""Reprocess authored monomers without losing their intended attachment sites."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace

from rdkit import Chem

from pyPept.attachments import attachment_sites
from pyPept.interfaces.monomer_pipeline import (
    activate_sites,
    backbone_exclusions,
    find_sidechain_slots,
)
from pyPept.leaving_groups import restore_leaving_groups
from pyPept.monomer_store import monomer_record
from pyPept.site_chemistry import canonical_atom_order


@dataclass(frozen=True)
class ReprocessedMonomer:
    molecule: object
    slot_mapping: dict[int, int]
    authored_sites: tuple[int, ...]


def restored_sites(molecule):
    """Restore leaving groups while retaining each site's original anchor."""
    mol = Chem.Mol(molecule)
    sites = attachment_sites(mol)
    anchors = {}
    for atom in mol.GetAtoms():
        atom.SetIntProp("_migration_origin", atom.GetIdx())
        if atom.GetAtomicNum() == 0:
            slot = atom.GetIsotope()
            if slot in anchors:
                raise ValueError(f"Duplicate attachment R{slot}")
            anchors[slot] = atom.GetNeighbors()[0].GetIdx()
    leaving = {site["slot"]: site["leaving"] or "[H]" for site in sites}
    free = restore_leaving_groups(mol, leaving)
    positions = {
        atom.GetIntProp("_migration_origin"): atom.GetIdx() for atom in free.GetAtoms()
    }
    return free, {slot: positions[index] for slot, index in anchors.items()}, sites


def reprocess_monomer(molecule):
    """Use the current detector, with the stored definition supplying intent.

    Existing R1/R2 anchors choose the backbone or cap. Existing sidechain ports
    remain mandatory, including authored handles outside automatic perception.
    Detection adds newly supported sites. R3 is reserved for a second backbone
    N hydrogen; all other ports follow the same canonical order as raw input.
    Every old port must retain its anchor, leaving group and reaction chemistry.
    """
    free, anchors, original = restored_sites(molecule)
    mol = Chem.AddHs(free)
    backbone = {slot: index for slot, index in anchors.items() if slot in (1, 2)}
    order = canonical_atom_order(
        free, {index: slot for slot, index in backbone.items()}
    )
    known = {site["slot"]: site for site in original}
    # (anchor, leaving group, reaction chemistry); duplicate H ports are real.
    ports = {
        slot: (index, known[slot]["leaving"] or "[H]", known[slot]["chem_type"])
        for slot, index in anchors.items()
    }
    reserved = dict(backbone)
    if set(backbone) == {1, 2} or (1 in backbone and anchors.get(3) == backbone[1]):
        atom = mol.GetAtomWithIdx(backbone[1])
        hydrogens = sum(nb.GetAtomicNum() == 1 for nb in atom.GetNeighbors())
        if atom.GetAtomicNum() == 7 and hydrogens >= 2:
            reserved[3] = backbone[1]
    assigned, mapping = {}, {}
    for slot, anchor in reserved.items():
        kind = known[slot]["chem_type"] if slot in backbone else "backbone_n_mod"
        old_slot = next(
            (
                old
                for old, port in ports.items()
                if old not in mapping
                and port[0] == anchor
                and (old == slot or slot == 3 and port[1] == "[H]")
            ),
            None,
        )
        if old_slot is not None:
            assigned[slot] = ports[old_slot]
            mapping[old_slot] = slot
        else:
            assigned[slot] = (anchor, "[H]", kind)

    detected = list(
        find_sidechain_slots(
            mol,
            backbone_exclusions(mol, backbone),
            atom_order=order,
        ).values()
    )
    remaining = [(old, port) for old, port in ports.items() if old not in mapping]
    remaining, alignment = _align_symmetric_sites(free, remaining, detected, backbone)
    # Combine capacities rather than creating an extra dummy for the same H.
    capacities = Counter((index, group) for index, group, _ in detected)
    authored = []
    for old, port in remaining:
        key = port[:2]
        if capacities[key]:
            capacities[key] -= 1
            index = next(
                i for i, candidate in enumerate(detected) if candidate[:2] == key
            )
            detected.pop(index)
        else:
            authored.append(old)
    combined = remaining + [(None, port) for port in detected]
    # Additional capacity on one anchor follows its existing equivalent ports.
    combined.sort(
        key=lambda item: (order[item[1][0]], item[1][1], item[1][2],
                          item[0] is None, item[0] or 0)
    )
    for slot, (old, port) in enumerate(combined, start=4):
        assigned[slot] = port
        if old is not None:
            mapping[old] = slot
    # Keep the authored template itself: restoring a prochiral C-H port first
    # would erase the configuration specified for its later bonded product.
    authored_template = Chem.Mol(molecule)
    for atom in authored_template.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetIsotope(mapping[atom.GetIsotope()])
    inverse = {new: old for old, new in enumerate(alignment)}
    result = activate_sites(
        Chem.AddHs(authored_template),
        {
            slot: free.GetAtomWithIdx(inverse[port[0]]).GetIntProp("_migration_origin")
            for slot, port in assigned.items()
        },
        {slot: port[1] for slot, port in assigned.items()},
        {slot: port[2] for slot, port in assigned.items()},
    )
    for old, new in mapping.items():
        if result.chem_types[new] != known[old]["chem_type"]:
            raise ValueError(f"R{old} → R{new} changes attachment chemistry")
    if set(mapping) != set(ports):
        raise ValueError("Reprocessing lost an attachment site")
    processed = monomer_record(
        result.chuckles,
        molecule.GetProp("symbol"),
        result.leaving,
        result.chem_types,
        name=molecule.GetProp("m_name"),
        m_type=molecule.GetProp("m_type"),
        m_subtype=molecule.GetProp("m_subtype"),
        activation_policy=result.policy,
    )
    additional = set(assigned) - set(mapping.values())
    previous = restore_leaving_groups(processed, result.leaving, slots=additional)
    if Chem.MolToSmiles(previous) != Chem.MolToSmiles(authored_template):
        raise ValueError("Reprocessing changed the authored attachment template")
    return ReprocessedMonomer(
        processed, mapping, tuple(sorted(mapping[slot] for slot in authored))
    )


def _align_symmetric_sites(free, original, detected, backbone):
    """Match equivalent authored handles before merging detector capacities.

    A maleimide's two equivalent alkene carbons must not become two handles
    merely because its author and detector chose opposite sides of the ring.
    One whole-graph symmetry maps every site together, with the backbone fixed.
    """
    capacities = Counter(port[:2] for port in detected)

    def overlap(ports):
        return sum((Counter(port[:2] for _, port in ports) & capacities).values())

    best, score = original, overlap(original)
    alignment = tuple(range(free.GetNumAtoms()))
    if score == len(original):
        return original, alignment
    matches = free.GetSubstructMatches(
        free, useChirality=True, uniquify=False, maxMatches=10001
    )
    for match in matches:
        if any(match[index] != index for index in backbone.values()):
            continue
        candidate = [
            (old, (match[index], group, kind)) for old, (index, group, kind) in original
        ]
        candidate_score = overlap(candidate)
        if candidate_score > score:
            best, score, alignment = candidate, candidate_score, match
        if score == len(original):
            return best, alignment
    if len(matches) > 10000:
        raise ValueError(
            "Molecular symmetry exceeds the reprocessing limit; keep the explicit template"
        )
    return best, alignment


def migrate_sequence(source, old_library, new_library, slot_mappings):
    """Rewrite CABILN connections against two explicit libraries; prove identity.

    Unchanged connections retain their spelling. Renumbered connections are
    serialized from the resolved graph, including nested and cross-chain bonds.
    The caller must supply the library under which the old source was written.
    """
    from pyPept.molecule import Molecule
    from pyPept.peptide import Connection, Endpoint, Peptide, serialize
    from pyPept.sequence import Sequence
    from pathlib import Path

    old_library, new_library = Path(old_library).resolve(), Path(new_library).resolve()
    if not old_library.is_file() or not new_library.is_file():
        raise ValueError("Migration requires both original and resulting SDF files")

    def parse(text, library):
        return Sequence(
            text,
            monomer_lib=str(library),
            track_source=True,
            warning_sink=lambda _: None,
        )

    old = parse(source, old_library)
    peptide = Peptide.from_sequence(old)
    nodes, mappings, definitions = [], {}, {}
    for node in peptide.occurrences:
        symbol = node.definition.symbol
        mappings[node.id] = {
            int(k): int(v) for k, v in slot_mappings.get(symbol, {}).items()
        }
        key = symbol, node.token
        if key not in definitions:
            definitions[key] = Peptide.from_sequence(
                parse(node.token, new_library)
            ).occurrences[0]
        nodes.append(replace(definitions[key], id=node.id, token=node.token))

    def endpoint(old):
        return Endpoint(
            old.occurrence_id, mappings[old.occurrence_id].get(old.slot, old.slot)
        )

    edges = tuple(
        Connection(endpoint(edge.left), endpoint(edge.right), edge.label)
        for edge in peptide.connections
    )
    expected = Peptide(tuple(nodes), edges)
    text = source
    if edges != peptide.connections:
        text = serialize(
            expected, "bracket" if ".[" in source or ".{" in source else "percent"
        ).text
    actual = parse(text, new_library)
    before = Chem.MolToSmiles(Molecule(old, depiction=None).mol)
    intended = Chem.MolToSmiles(Molecule(expected, depiction=None).mol)
    after = Chem.MolToSmiles(Molecule(actual, depiction=None).mol)
    if before != intended or before != after:
        raise ValueError("Site migration changed the assembled structure")
    return text


def reprocess_library(source, destination):
    """Write a new SDF, CSV and migration manifest; never replace the input."""
    import csv
    import json
    from hashlib import sha256
    from pathlib import Path

    from pyPept.library_quality import definition_hash

    source, destination = Path(source).resolve(), Path(destination).resolve()
    if destination.exists():
        raise ValueError("Choose a new output directory for the migrated library")
    originals = list(Chem.SDMolSupplier(str(source), removeHs=False))
    if not originals or any(mol is None for mol in originals):
        raise ValueError(
            "The source library is empty or contains unreadable definitions"
        )
    symbols = [mol.GetProp("symbol") for mol in originals]
    if len(set(symbols)) != len(symbols):
        raise ValueError("The source library has duplicate symbols")
    results = [reprocess_monomer(mol) for mol in originals]
    records = {
        name: {
            "before": definition_hash(old),
            "after": definition_hash(result.molecule),
            "slots": result.slot_mapping,
            "authored_sites": result.authored_sites,
        }
        for name, old, result in zip(symbols, originals, results)
    }
    rows = {}
    companion = source.with_name("monomers.csv")
    if companion.exists():
        with companion.open(newline="", encoding="utf-8") as stream:
            rows = {row["token"]: row for row in csv.DictReader(stream)}
    output_rows = []
    for name, old, result in zip(symbols, originals, results):
        mol = result.molecule
        row = dict(rows.get(name, {}))
        row.update(
            token=name,
            name=mol.GetProp("m_name"),
            type=mol.GetProp("m_type"),
            subtype=mol.GetProp("m_subtype"),
            chuckles=Chem.MolToSmiles(mol),
            chem_types=mol.GetProp("m_chem_types"),
            activation_policy=mol.GetProp("m_activation_policy"),
        )
        if not row.get("input"):
            row["input"] = row["chuckles"]
        # Explicit inputs carry site numbers too; retain their authored intent.
        if "*" in row["input"]:
            row["input"] = row["chuckles"]
        for key in row:
            if key.startswith("r") and key.endswith("_leaving"):
                row[key] = ""
        for slot, group in enumerate(mol.GetProp("m_Rgroups").split(","), start=1):
            row[f"r{slot}_leaving"] = "" if group == "None" else group
        output_rows.append(row)
    destination.mkdir(parents=True)
    sdf = destination / "monomers.sdf"
    with Chem.SDWriter(str(sdf)) as writer:
        for result in results:
            writer.write(result.molecule)
    fields = list(dict.fromkeys(key for row in output_rows for key in row))
    with (destination / "monomers.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output_rows)
    manifest = {
        "version": 1,
        "policy": results[0].molecule.GetProp("m_activation_policy"),
        "source_sha256": sha256(source.read_bytes()).hexdigest(),
        "target_sha256": sha256(sdf.read_bytes()).hexdigest(),
        "source_aliases_sha256": (
            sha256(companion.read_bytes()).hexdigest() if companion.exists() else None
        ),
        "target_aliases_sha256": sha256(
            (destination / "monomers.csv").read_bytes()
        ).hexdigest(),
        "records": records,
    }
    (destination / "site-migration.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    return manifest


def main():
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", required=True, type=Path, help="Original monomers.sdf"
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="New library directory"
    )
    parser.add_argument(
        "--sequence", help="Migrate this CABILN using an existing output directory"
    )
    args = parser.parse_args()
    if args.sequence is None:
        manifest = reprocess_library(args.source, args.output)
        print(f"Reprocessed {len(manifest['records'])} monomers into {args.output}")
    else:
        from hashlib import sha256

        manifest = json.loads((args.output / "site-migration.json").read_text())
        target = args.output / "monomers.sdf"
        if (
            manifest["source_sha256"] != sha256(args.source.read_bytes()).hexdigest()
            or manifest["target_sha256"] != sha256(target.read_bytes()).hexdigest()
        ):
            raise ValueError("The libraries do not match the migration manifest")
        for path, key in (
            (args.source.with_name("monomers.csv"), "source_aliases_sha256"),
            (args.output / "monomers.csv", "target_aliases_sha256"),
        ):
            actual = sha256(path.read_bytes()).hexdigest() if path.exists() else None
            if manifest.get(key) != actual:
                raise ValueError("The aliases do not match the migration manifest")
        mappings = {
            name: record["slots"] for name, record in manifest["records"].items()
        }
        print(
            migrate_sequence(
                args.sequence, args.source.resolve(), target.resolve(), mappings
            )
        )


if __name__ == "__main__":
    main()
