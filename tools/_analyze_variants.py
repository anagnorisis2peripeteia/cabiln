#!/usr/bin/env python3
"""Analyze monomers.sdf for directional variants (N-term / C-term pairs).

Reads SDF properties: m_abbr, m_type, m_subtype, m_Rgroups, m_chem_types
Identifies groups that appear to be directional variants of the same base monomer.
"""

import re
from pathlib import Path
from collections import defaultdict

SDF_PATH = Path(__file__).resolve().parent.parent / "src" / "pyPept" / "data" / "monomers.sdf"


def parse_sdf(path: Path) -> list[dict]:
    """Parse an SDF file into a list of monomer property dicts."""
    monomers = []
    current_props = {}
    current_prop_name = None
    reading_value = False

    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")

            # Detect property header like ">  <m_abbr>  (1) "
            prop_match = re.match(r">\s+<(\w+)>\s+\(\d+\)", line)
            if prop_match:
                current_prop_name = prop_match.group(1)
                reading_value = True
                continue

            # Read property value (next non-empty line after header)
            if reading_value:
                if line.strip() == "":
                    if current_prop_name and current_prop_name in current_props:
                        reading_value = False
                        current_prop_name = None
                    continue
                current_props[current_prop_name] = line.strip()
                reading_value = False
                current_prop_name = None
                continue

            # End of molecule
            if line.startswith("$$$$"):
                if current_props:
                    monomers.append(current_props)
                current_props = {}
                continue

    # Don't forget last entry
    if current_props:
        monomers.append(current_props)

    return monomers


def strip_directional(abbr: str) -> tuple[str, str]:
    """Strip directional prefix/suffix and return (base_name, direction).

    Convention in pyPept:
      - Leading underscore  (_Foo)  = C-terminal cap
      - Trailing underscore (Foo_)  = N-terminal cap
    Also handles explicit Nterm/Cterm if present.
    """
    # Explicit Nterm/Cterm suffixes
    for tag in ("_Nterm", "_Cterm", "-Nterm", "-Cterm"):
        if abbr.endswith(tag):
            direction = "N-term" if "Nterm" in tag else "C-term"
            return abbr[: -len(tag)], direction

    # Leading underscore = C-term cap
    if abbr.startswith("_") and not abbr.endswith("_"):
        return abbr[1:], "C-term"

    # Trailing underscore = N-term cap
    if abbr.endswith("_") and not abbr.startswith("_"):
        return abbr[:-1], "N-term"

    return abbr, "none"


def main():
    monomers = parse_sdf(SDF_PATH)
    print(f"Total monomers parsed: {len(monomers)}\n")

    # ── 1. Extract all abbreviations and identify directional ones ──
    directional = []
    non_directional = []
    all_abbrs = []

    for m in monomers:
        abbr = m.get("m_abbr", "???")
        all_abbrs.append(abbr)
        base, direction = strip_directional(abbr)
        if direction != "none":
            directional.append((abbr, base, direction, m))
        else:
            non_directional.append((abbr, m))

    # Group directional monomers by base name
    groups = defaultdict(list)
    for abbr, base, direction, m in directional:
        groups[base].append((abbr, direction, m))

    print("=" * 80)
    print("DIRECTIONAL VARIANTS (monomers with leading/trailing underscore or Nterm/Cterm)")
    print("=" * 80)
    print(f"\nFound {len(directional)} directional monomers across {len(groups)} base names.\n")

    for base, members in sorted(groups.items()):
        has_pair = len(members) > 1
        pair_label = "PAIRED" if has_pair else "UNPAIRED (single direction only)"
        print(f"--- Base: {base!r}  [{pair_label}] ---")
        for abbr, direction, m in sorted(members, key=lambda x: x[1]):
            mtype = m.get("m_type", "?")
            msubtype = m.get("m_subtype", "?")
            rgroups = m.get("m_Rgroups", "?")
            chem_types = m.get("m_chem_types", "(not set)")
            mname = m.get("m_name", "?")
            print(f"  abbr:       {abbr}")
            print(f"  name:       {mname}")
            print(f"  direction:  {direction}")
            print(f"  m_type:     {mtype}")
            print(f"  m_subtype:  {msubtype}")
            print(f"  m_Rgroups:  {rgroups}")
            print(f"  m_chem_types: {chem_types}")

            # Parse R-groups to show which slots are active
            rg_parts = [x.strip() for x in rgroups.split(",")]
            active = [(f"R{i+1}", rg_parts[i]) for i in range(len(rg_parts)) if rg_parts[i] != "None"]
            print(f"  active R:   {active}")
            print()

    # ── 2. Find monomers that also exist as non-directional ──
    non_dir_set = {abbr for abbr, _ in non_directional}
    print("\n" + "=" * 80)
    print("BASE NAMES THAT ALSO EXIST AS STANDALONE (non-directional) MONOMERS")
    print("=" * 80 + "\n")
    for base in sorted(groups.keys()):
        if base in non_dir_set:
            # Get the standalone monomer details
            standalone = [m for a, m in non_directional if a == base][0]
            print(f"  {base}: standalone exists (type={standalone.get('m_type','?')}, "
                  f"subtype={standalone.get('m_subtype','?')}, "
                  f"Rgroups={standalone.get('m_Rgroups','?')})")
        else:
            print(f"  {base}: NO standalone entry")

    # ── 3. Caps analysis ──
    caps = [(m.get("m_abbr", "???"), m) for m in monomers if m.get("m_type") == "cap"]

    print("\n" + "=" * 80)
    print("ALL CAPS")
    print("=" * 80)
    print(f"\nTotal caps: {len(caps)}\n")

    nterm_caps = []
    cterm_caps = []
    ambiguous_caps = []

    for abbr, m in caps:
        _, direction = strip_directional(abbr)
        rgroups = m.get("m_Rgroups", "")
        rg_parts = [x.strip() for x in rgroups.split(",")]
        r1 = rg_parts[0] if len(rg_parts) > 0 else "None"
        r2 = rg_parts[1] if len(rg_parts) > 1 else "None"

        # Determine actual terminal attachment from R-groups
        has_r1 = r1 != "None"
        has_r2 = r2 != "None"

        if has_r1 and not has_r2:
            actual_dir = "C-term (R1 only)"
        elif has_r2 and not has_r1:
            actual_dir = "N-term (R2 only)"
        elif has_r1 and has_r2:
            actual_dir = "both R1+R2"
        else:
            actual_dir = "neither R1 nor R2"

        active = [(f"R{i+1}", rg_parts[i]) for i in range(len(rg_parts)) if rg_parts[i] != "None"]

        print(f"  {abbr:25s}  naming_dir={direction:8s}  rgroup_dir={actual_dir:20s}  active={active}")

        if "C-term" in actual_dir:
            cterm_caps.append(abbr)
        elif "N-term" in actual_dir:
            nterm_caps.append(abbr)
        else:
            ambiguous_caps.append(abbr)

    print(f"\n  C-term caps (R1 only, attach to chain C-terminus): {len(cterm_caps)}")
    print(f"  N-term caps (R2 only, attach to chain N-terminus): {len(nterm_caps)}")
    print(f"  Ambiguous (both or neither R1/R2): {len(ambiguous_caps)}")

    # ── 4. Summary: caps with directional variants ──
    dir_cap_bases = {base for base, members in groups.items()
                     if any(m.get("m_type") == "cap" for _, _, m in members)}

    print("\n" + "=" * 80)
    print("CAPS WITH DIRECTIONAL VARIANTS (paired N-term + C-term)")
    print("=" * 80 + "\n")

    for base in sorted(dir_cap_bases):
        members = groups[base]
        print(f"  {base}:")
        for abbr, direction, m in sorted(members, key=lambda x: x[1]):
            rgroups = m.get("m_Rgroups", "?")
            mname = m.get("m_name", "?")
            print(f"    {abbr:20s} {direction:8s}  Rgroups={rgroups}  name={mname}")
        print()

    # ── 5. R-group topology explanation ──
    print("=" * 80)
    print("WHY SEPARATE ENTRIES?")
    print("=" * 80 + "\n")
    print("Caps need separate N-term and C-term entries because they have different")
    print("R-group topologies:")
    print()
    print("  N-terminal caps (e.g., Ac_, Bn_):")
    print("    - Have R2 (C-terminal attachment) but NOT R1")
    print("    - R2 provides the bond to the first residue's nitrogen")
    print("    - Convention: trailing underscore (Foo_)")
    print()
    print("  C-terminal caps (e.g., _NH2, _Bn):")
    print("    - Have R1 (N-terminal attachment) but NOT R2")
    print("    - R1 provides the bond from the last residue's carbonyl")
    print("    - Convention: leading underscore (_Foo)")
    print()
    print("  A cap like benzylamine (NHBn / Bn) needs two SDF entries because")
    print("  the molecular structure differs depending on which terminus it caps:")
    print("    NHBn_  (N-term): PhCH2-NH- attaching via R2=[OH] to first residue")
    print("    _NHBn  (C-term): -NH-CH2Ph attaching via R1=[H] from last residue")
    print()
    print("  Similarly, toluene/benzyl (Bn):")
    print("    Bn_  (N-term): PhCH2- via R2=[OH]")
    print("    _Bn  (C-term): -CH2Ph via R1=[H]")
    print()

    # ── 6. All m_type values and counts ──
    type_counts = defaultdict(int)
    for m in monomers:
        type_counts[m.get("m_type", "???")] += 1

    print("=" * 80)
    print("MONOMER TYPE DISTRIBUTION")
    print("=" * 80 + "\n")
    for t, c in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"  {t:20s}: {c}")
    print()


if __name__ == "__main__":
    main()
