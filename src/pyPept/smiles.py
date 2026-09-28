"""Convert molecular SMILES to CABILN using the monomer library. No web dependencies."""

from __future__ import annotations

import re as _re
import threading
from collections import namedtuple as _nt2
from dataclasses import dataclass

from pyPept.monomer_store import _load_sdf, library_path
from pyPept.structure import compare_structures, require_supported_stereo


class StereoInferenceWarning(UserWarning):
    """The input omitted stereo that the selected library monomers supply."""


@dataclass(frozen=True)
class ConversionResult:
    """Request-local conversion output, including any limitations or inference.

    ``synthetic_components`` contains zero-based input component indices that
    needed synthetic notation. Warnings distinguish local residue templates
    from coarse component fallbacks; a coarse fallback has no residue details.
    """

    cabiln: str
    details: list[tuple[str, float, int]]
    inferred_stereo: bool
    synthetic_components: tuple[int, ...]
    warnings: tuple[str, ...]


_STEREO_INFERENCE_MESSAGE = (
    "Input leaves stereochemistry unspecified; the monomer library supplies "
    "additional stereochemistry."
)


def _s2c_library_stamp():
    """Keep every derived library cache in step with the selected SDF."""
    path = library_path()
    stat = path.stat()
    return str(path), stat.st_mtime_ns, stat.st_size


# Leaving-group SDF label (e.g. '[OH]') → SMILES token for regex substitution.
# '[H]' is the implicit fallback and is not listed here.
_LG_LABEL_TO_SMI = {"[OH]": "O", "[SH]": "S", "[Cl]": "Cl", "[Br]": "Br", "[I]": "I"}

# Leaving-group SDF label → atomic number for RWMol.ReplaceAtom paths.
_LG_LABEL_TO_ANUM = {"[OH]": 8, "[SH]": 16, "[Cl]": 17, "[Br]": 35, "[I]": 53}


def _resolve_lg_smiles(rg):
    """Return SMILES token for leaving-group label *rg*, or '[H]' as fallback."""
    return _LG_LABEL_TO_SMI.get(rg, "[H]")


def _resolve_dummy_rgroups(rw2, rg_list, skip_iso):
    """Replace or remove all dummy atoms in *rw2* except the one with isotope *skip_iso*.

    Iterates highest-index first so atom removal does not invalidate earlier indices.
    Mutates *rw2* in place.
    """
    from rdkit.Chem import Atom as _Atom

    dummies = [
        (a.GetIdx(), a.GetIsotope())
        for a in rw2.GetAtoms()
        if a.GetAtomicNum() == 0 and a.GetIsotope() != skip_iso
    ]
    for d_idx, iso in sorted(dummies, key=lambda x: -x[0]):
        rg = rg_list[iso - 1].strip() if 0 < iso <= len(rg_list) else None
        anum = _LG_LABEL_TO_ANUM.get(rg)
        if anum is not None:
            rw2.ReplaceAtom(d_idx, _Atom(anum))
        else:
            rw2.RemoveAtom(d_idx)


def _s2c_normalize(mol):
    """Normalize resonance forms (e.g. guanidinium) via InChI round-trip."""
    try:
        from rdkit import Chem as _Chem
        from rdkit.Chem.inchi import MolFromInchi, MolToInchi

        inchi = MolToInchi(mol)
        if inchi:
            m = MolFromInchi(inchi)
            if m:
                _Chem.SanitizeMol(m)
                return m
    except (ValueError, RuntimeError):
        pass
    return mol


def _s2c_cap_smiles(smi, rgroups_str):
    """Replace [n*] dummy atoms with leaving-group atoms; fallback to [H]."""
    rgroups = [r.strip() for r in rgroups_str.split(",")]

    def _rep(m):
        n = int(m.group(1))
        rg = rgroups[n - 1] if n - 1 < len(rgroups) else "None"
        return _resolve_lg_smiles(rg)

    return _re.sub(r"\[(\d+)\*\]", _rep, smi)


def _s2c_isolate_residues(m, aa_units):
    """Cut each residue from the peptide ring; cap C=O termini with -OH."""
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol

    aas, mappings = [], []
    for atom_idxs in aa_units:
        rw = RWMol()
        amap = {}
        for idx in atom_idxs:
            amap[idx] = rw.AddAtom(m.GetAtomWithIdx(idx))
        mappings.append(amap)
        # First pass: add all bonds (no stereo yet)
        for bond in m.GetBonds():
            bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bi in atom_idxs and ei in atom_idxs:
                rw.AddBond(amap[bi], amap[ei], bond.GetBondType())
            elif bi in atom_idxs or ei in atom_idxs:
                inside = bi if bi in atom_idxs else ei
                outside = ei if inside == bi else bi
                ia = m.GetAtomWithIdx(inside)
                oa = m.GetAtomWithIdx(outside)
                if ia.GetSymbol() == "C" and oa.GetSymbol() in ("N", "O"):
                    has_dbl_o = any(
                        nb.GetSymbol() == "O"
                        and m.GetBondBetweenAtoms(inside, nb.GetIdx()).GetBondType()
                        == _C.BondType.DOUBLE
                        for nb in ia.GetNeighbors()
                    )
                    if has_dbl_o:
                        o_idx = rw.AddAtom(_C.Atom("O"))
                        rw.AddBond(amap[inside], o_idx, _C.BondType.SINGLE)
        # Second pass: restore E/Z stereo on double bonds (stereoAtoms must be bonded first)
        from rdkit.Chem.rdchem import BondStereo as _BS

        for bond in m.GetBonds():
            bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if (
                bi in atom_idxs
                and ei in atom_idxs
                and bond.GetBondTypeAsDouble() == 2.0
            ):
                db_stereo = bond.GetStereo()
                if db_stereo not in (_BS.STEREONONE, _BS.STEREOANY):
                    sa = list(bond.GetStereoAtoms())
                    if len(sa) == 2 and sa[0] in atom_idxs and sa[1] in atom_idxs:
                        try:
                            nb = rw.GetBondBetweenAtoms(amap[bi], amap[ei])
                            nb.SetStereo(db_stereo)
                            nb.SetStereoAtoms(amap[sa[0]], amap[sa[1]])
                        except (ValueError, RuntimeError):
                            pass
        mol_out = rw.GetMol()
        try:
            _C.SanitizeMol(mol_out)
        except (ValueError, RuntimeError):
            # An aromatic ring atom (e.g. triazole N from a CuAAC crosslink) was
            # cut from its ring by the extraction, making it aromatic-but-not-in-ring.
            # Strip all aromaticity from the fragment and re-sanitize.
            _rw2 = _C.RWMol(mol_out)
            for _at2 in _rw2.GetAtoms():
                _at2.SetIsAromatic(False)
            for _bd2 in _rw2.GetBonds():
                if _bd2.GetBondTypeAsDouble() == 1.5:
                    _bd2.SetBondType(_C.BondType.SINGLE)
            mol_out = _rw2.GetMol()
            _C.SanitizeMol(
                mol_out,
                _C.SanitizeFlags.SANITIZE_ALL ^ _C.SanitizeFlags.SANITIZE_KEKULIZE,
            )
        # Fix chirality: AddAtom copies ChiralTags, but they are relative to
        # the original molecule's neighbor ordering.  In the fragment the
        # neighbor order may differ, flipping the effective CW/CCW meaning.
        # Correct each chiral centre so its CIP code matches the original.
        from rdkit.Chem import AllChem as _AC2
        from rdkit.Chem import rdchem as _RC2

        _AC2.AssignStereochemistry(mol_out, cleanIt=True, force=True)
        inv_map = {v: k for k, v in amap.items()}
        needs_fix = False
        for frag_idx, orig_idx in inv_map.items():
            orig_cip = m.GetAtomWithIdx(orig_idx).GetPropsAsDict().get("_CIPCode")
            if not orig_cip:
                continue
            frag_cip = mol_out.GetAtomWithIdx(frag_idx).GetPropsAsDict().get("_CIPCode")
            if frag_cip and frag_cip != orig_cip:
                needs_fix = True
                break
        if needs_fix:
            rw2 = _C.RWMol(mol_out)
            for frag_idx, orig_idx in inv_map.items():
                orig_cip = m.GetAtomWithIdx(orig_idx).GetPropsAsDict().get("_CIPCode")
                if not orig_cip:
                    continue
                frag_cip = (
                    mol_out.GetAtomWithIdx(frag_idx).GetPropsAsDict().get("_CIPCode")
                )
                if frag_cip and frag_cip != orig_cip:
                    ct = rw2.GetAtomWithIdx(frag_idx).GetChiralTag()
                    if ct == _RC2.ChiralType.CHI_TETRAHEDRAL_CW:
                        rw2.GetAtomWithIdx(frag_idx).SetChiralTag(
                            _RC2.ChiralType.CHI_TETRAHEDRAL_CCW
                        )
                    elif ct == _RC2.ChiralType.CHI_TETRAHEDRAL_CCW:
                        rw2.GetAtomWithIdx(frag_idx).SetChiralTag(
                            _RC2.ChiralType.CHI_TETRAHEDRAL_CW
                        )
            mol_out = rw2.GetMol()
        aas.append(mol_out)
    return aas, mappings


def _s2c_strip_n_cap(mol):
    """Remove N-terminal carbonyl-type cap (ac, fmoc, Boc) from backbone N.

    Only strips substituents whose first atom is a carbonyl C (has a =O).
    N-methyl groups (plain CH3) are preserved.
    Returns (stripped_mol, did_strip).
    """
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol

    patt = _C.MolFromSmarts("[N]-[C]-[C](=O)")
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return mol, False

    n_idx, ca_idx = matches[0][0], matches[0][1]
    n_atom = mol.GetAtomWithIdx(n_idx)

    cap_starts = []
    for nb in n_atom.GetNeighbors():
        if nb.GetIdx() == ca_idx:
            continue
        if nb.GetSymbol() == "C":
            has_carbonyl = any(
                x.GetSymbol() == "O"
                and mol.GetBondBetweenAtoms(
                    nb.GetIdx(), x.GetIdx()
                ).GetBondTypeAsDouble()
                == 2.0
                for x in nb.GetNeighbors()
            )
            if has_carbonyl:
                cap_starts.append(nb.GetIdx())

    if not cap_starts:
        return mol, False

    cap_atoms = set()
    queue = list(cap_starts)
    while queue:
        ai = queue.pop()
        if ai in cap_atoms or ai == n_idx:
            continue
        cap_atoms.add(ai)
        for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
            if nb.GetIdx() not in cap_atoms and nb.GetIdx() != n_idx:
                queue.append(nb.GetIdx())

    rw = RWMol(mol)
    for ai in sorted(cap_atoms, reverse=True):
        rw.RemoveAtom(ai)
    try:
        _C.SanitizeMol(rw)
        return rw.GetMol(), True
    except (ValueError, RuntimeError):
        return mol, False


def _s2c_strip_c_cap(mol):
    """Convert C-terminal cap to carboxylic acid (-COOH) for residue matching.

    Handles am (NH2), ester (OR), and other amide/ester C-caps by finding the
    terminal NCC=O, removing the cap substituent on the carbonyl C, and
    replacing it with OH.  Only touches the BACKBONE terminal C so
    Asn/Gln side-chain amides are safe.
    Returns (converted_mol, did_convert).
    """
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol

    patt = _C.MolFromSmarts("[N]-[C]-[C](=O)")
    matches = mol.GetSubstructMatches(patt)
    for match in matches:
        ca_idx, co_idx = match[1], match[2]
        co_atom = mol.GetAtomWithIdx(co_idx)
        already_acid = False
        cap_root = None
        for nb in co_atom.GetNeighbors():
            ni = nb.GetIdx()
            if ni == ca_idx:
                continue
            bd = mol.GetBondBetweenAtoms(co_idx, ni)
            if nb.GetSymbol() == "O" and bd.GetBondTypeAsDouble() == 2.0:
                continue  # carbonyl =O, skip
            if (
                nb.GetSymbol() == "O"
                and bd.GetBondTypeAsDouble() == 1.0
                and nb.GetTotalNumHs() > 0
            ):
                already_acid = True
                break
            cap_root = ni
        if already_acid or cap_root is None:
            continue
        # BFS from cap_root to collect all cap atoms
        cap_atoms: set = set()
        queue = [cap_root]
        while queue:
            ai = queue.pop()
            if ai in cap_atoms or ai == co_idx:
                continue
            cap_atoms.add(ai)
            for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni2 = nb.GetIdx()
                if ni2 not in cap_atoms and ni2 != co_idx:
                    queue.append(ni2)
        rw = RWMol(mol)
        rw.ReplaceAtom(cap_root, _C.Atom("O"))
        for ai in sorted(cap_atoms - {cap_root}, reverse=True):
            rw.RemoveAtom(ai)
        try:
            _C.SanitizeMol(rw)
            return rw.GetMol(), True
        except (ValueError, RuntimeError):
            continue

    return mol, False


# Capped library cache (invalidates when SDF changes)
_s2c_lib_cache: list = []
_s2c_lib_stamp: tuple | None = None
_s2c_lib_lock = threading.Lock()

# Full capped library (includes non-backbone monomers, for branch segment matching)
_s2c_full_lib_cache: list = []
_s2c_full_lib_stamp: tuple | None = None

# Scaffold crosslink patterns (multi-arm non-amide linkers, e.g. TBMB)
_s2c_scaffold_cache: list | None = None
_s2c_scaffold_stamp: tuple | None = None

# Ring-forming crosslink patterns derived from reactions.yaml product SMARTS
_s2c_ring_xlink_cache: list | None = None


def _s2c_get_lib():
    """Return (or rebuild) the capped monomer library."""
    from rdkit import Chem as _C

    global _s2c_lib_cache, _s2c_lib_stamp
    stamp = _s2c_library_stamp()
    with _s2c_lib_lock:
        if _s2c_lib_cache and _s2c_lib_stamp == stamp:
            return _s2c_lib_cache
        raw_mols, _ = _load_sdf()
        lib = []
        seen = set()
        for mol in raw_mols:
            props = mol.GetPropsAsDict()
            abbr = props.get("m_abbr", "").strip()
            if not abbr or abbr in seen:
                continue
            ct_str = props.get("m_chem_types", "")
            rg_str = props.get("m_Rgroups", "")
            cts = {}
            for item in ct_str.split(","):
                item = item.strip()
                if ":" in item:
                    k, v = item.split(":", 1)
                    try:
                        cts[int(k)] = v.strip()
                    except ValueError:
                        pass
            bn_r = next(
                (r for r, t in cts.items() if t in ("backbone_n", "backbone_o")), None
            )
            bc_r = next(
                (
                    r
                    for r, t in cts.items()
                    if t in ("backbone_c", "backbone_c_red", "sp3_c_anchor")
                ),
                None,
            )
            if bn_r is None or bc_r is None:
                continue
            bn_atom = next(
                (
                    a
                    for a in mol.GetAtoms()
                    if a.GetAtomicNum() == 0 and a.GetIsotope() == bn_r
                ),
                None,
            )
            bc_atom = next(
                (
                    a
                    for a in mol.GetAtoms()
                    if a.GetAtomicNum() == 0 and a.GetIsotope() == bc_r
                ),
                None,
            )
            if bn_atom is None or bc_atom is None:
                continue
            bn_nbrs = [nb.GetIdx() for nb in bn_atom.GetNeighbors()]
            bc_nbrs = [nb.GetIdx() for nb in bc_atom.GetNeighbors()]
            if not bn_nbrs or not bc_nbrs:
                continue
            _bb_path = _C.GetShortestPath(mol, bn_nbrs[0], bc_nbrs[0])
            bb_dist = len(_bb_path) - 1 if _bb_path else 0
            smi = _C.MolToSmiles(mol, allHsExplicit=False)
            capped_smi = _s2c_cap_smiles(smi, rg_str)
            try:
                cm = _C.MolFromSmiles(capped_smi)
            except (ValueError, RuntimeError):
                cm = None
            if cm is None:
                continue
            cm = _C.RemoveHs(cm)
            orig_cm = cm  # preserve stereochemistry
            norm_cm = _s2c_normalize(cm)  # resonance-normalised (strips stereo)
            if norm_cm is None:
                continue
            from rdkit.Chem import rdchem as _RC
            from rdkit.Chem import rdMolDescriptors as _RDM

            n_rings = _RDM.CalcNumRings(orig_cm)
            has_ez = any(
                b.GetStereo()
                not in (_RC.BondStereo.STEREONONE, _RC.BondStereo.STEREOANY)
                for b in orig_cm.GetBonds()
                if b.GetBondTypeAsDouble() == 2.0
            )
            # Count R-group slots beyond backbone (R3, R4, ...) — used by
            # _s2c_match as a tiebreaker so monomers with crosslink slots win
            # over plain aliases when both match the same atoms (e.g. XlQuat
            # over Ile for chondramide's quat-C residue).
            n_rgroups_extra = max(0, len(cts) - 2)
            lib.append(
                (
                    abbr,
                    orig_cm,
                    norm_cm,
                    norm_cm.GetNumAtoms(),
                    n_rings,
                    has_ez,
                    bb_dist,
                    n_rgroups_extra,
                )
            )
            seen.add(abbr)
        lib.sort(key=lambda x: x[3], reverse=True)
        _s2c_lib_cache = lib
        _s2c_lib_stamp = stamp
        return lib


def _s2c_get_raw_lib():
    """Return {abbr: raw_mol} with R-atom dummy atoms (isotope labels) preserved."""
    _, by_abbr = _load_sdf()
    return by_abbr


def _s2c_crosslink_r(abbr, raw_lib):
    """Return the non-backbone R-group number for a crosslink monomer.

    Parses m_chem_types and returns the first R-group that is not backbone_n,
    backbone_c, or backbone_n_mod.  Falls back to 4 if nothing found.
    """
    raw_mol = raw_lib.get(abbr)
    if raw_mol is None:
        return 4
    chem_types = raw_mol.GetPropsAsDict().get("m_chem_types", "")
    _backbone = {"backbone_n", "backbone_c", "backbone_n_mod"}
    for part in chem_types.split(","):
        part = part.strip()
        if ":" not in part:
            continue
        r_num, r_type = part.split(":", 1)
        if r_type.strip() not in _backbone:
            try:
                return int(r_num.strip())
            except ValueError:
                pass
    return 4


def _s2c_get_full_lib():
    """Capped library including non-backbone monomers (caps, linkers, etc.)."""
    from rdkit import Chem as _C

    global _s2c_full_lib_cache, _s2c_full_lib_stamp
    stamp = _s2c_library_stamp()
    with _s2c_lib_lock:
        if _s2c_full_lib_cache and _s2c_full_lib_stamp == stamp:
            return _s2c_full_lib_cache
        raw_mols, _ = _load_sdf()
        lib = []
        seen = set()
        for mol in raw_mols:
            props = mol.GetPropsAsDict()
            abbr = props.get("m_abbr", "").strip()
            if not abbr or abbr in seen:
                continue
            rg_str = props.get("m_Rgroups", "")
            smi = _C.MolToSmiles(mol, allHsExplicit=False)
            capped_smi = _s2c_cap_smiles(smi, rg_str)
            try:
                cm = _C.MolFromSmiles(capped_smi)
            except (ValueError, RuntimeError):
                cm = None
            if cm is None:
                continue
            cm = _C.RemoveHs(cm)
            norm_cm = _s2c_normalize(cm)
            if norm_cm is None:
                continue
            lib.append((abbr, cm, norm_cm, norm_cm.GetNumAtoms()))
            seen.add(abbr)
        lib.sort(key=lambda x: x[3], reverse=True)
        _s2c_full_lib_cache = lib
        _s2c_full_lib_stamp = stamp
        return lib


def _s2c_get_scaffold_patterns():
    """Return scaffold crosslink patterns for all multi-arm non-backbone linkers.

    A scaffold linker is any SDF monomer whose CHUCKLES has ≥ 2 dummy atoms with
    isotope ≥ 4 (non-backbone attachment points) and NO dummy atoms with isotopes
    1 or 2 (backbone N/C connections), making it a pure crosslink scaffold.

    Each entry: (abbr, smarts_mol, {smarts_atom_idx: r_group_num})
    where smarts_atom_idx is the index in the SMARTS mol for each attachment-point
    wildcard (the backbone atoms that bond into the scaffold in the assembled mol).

    Detection strategy: replace [n*] (n ≥ 4) in the CHUCKLES canonical SMILES with
    [*:n] (atom-mapped wildcard) so GetAtomMapNum() recovers the R-group number from
    the SMARTS match, then verify each matched atom sits inside a backbone node.
    """
    import re as _re2

    from rdkit import Chem as _C

    global _s2c_scaffold_cache, _s2c_scaffold_stamp
    stamp = _s2c_library_stamp()
    with _s2c_lib_lock:
        if _s2c_scaffold_cache is not None and _s2c_scaffold_stamp == stamp:
            return _s2c_scaffold_cache

    # Load raw lib OUTSIDE the lock — _s2c_get_raw_lib acquires the same
    # non-reentrant lock and would deadlock if called inside our with-block.
    raw_lib = _s2c_get_raw_lib()
    patterns = []

    for abbr, raw_mol in raw_lib.items():
        if raw_mol is None:
            continue
        atoms = list(raw_mol.GetAtoms())
        # Attachment points: dummy atoms (atomic num 0) with isotope ≥ 4
        attach = [
            (a.GetIdx(), a.GetIsotope())
            for a in atoms
            if a.GetAtomicNum() == 0 and a.GetIsotope() >= 4
        ]
        if len(attach) < 2:
            continue
        # Exclude backbone linkers (they have R1 or R2 = backbone N/C)
        if any(a.GetAtomicNum() == 0 and a.GetIsotope() in (1, 2) for a in atoms):
            continue

        try:
            chuckles = _C.MolToSmiles(raw_mol, canonical=True)
        except (ValueError, RuntimeError):
            continue

        # [n*] → [*:n] for n ≥ 4; leave backbone R-groups unchanged
        smarts_str = _re2.sub(
            r"\[(\d+)\*\]",
            lambda m: f"[*:{m.group(1)}]" if int(m.group(1)) >= 4 else m.group(0),
            chuckles,
        )
        smarts_mol = _C.MolFromSmarts(smarts_str)
        if smarts_mol is None:
            continue

        r_positions = {
            a.GetIdx(): a.GetAtomMapNum()
            for a in smarts_mol.GetAtoms()
            if a.GetAtomMapNum() >= 4
        }
        if len(r_positions) < 2:
            continue

        global_min_slot = min(iso for _, iso in attach)
        patterns.append((abbr, smarts_mol, r_positions, global_min_slot))

        # For thia_michael_c arms (no leaving group; unreacted form = vinyl CH2=CH-):
        # the full SMARTS CC[*:n] can't match the terminal vinyl because CH2= has no
        # non-H atom for [*:n] to bind.  Generate partial SMARTS for each such arm
        # by removing its [*:n] token, so partial products are detected correctly.
        mol_props = raw_mol.GetPropsAsDict()
        m_chem_str = mol_props.get("m_chem_types", "")
        arm_chem_types: dict = {}  # r_slot_int → chem_type_str
        for _part in m_chem_str.split(","):
            if ":" in _part:
                _slot_str, _ct = _part.strip().split(":", 1)
                try:
                    arm_chem_types[int(_slot_str)] = _ct.strip()
                except ValueError:
                    pass

        # Patch the full-match pattern tuple to include arm_chem_types and min_required.
        # Full pattern: allow N-1 arms to match (1 unreacted arm tolerated via relaxation).
        # is_partial=False: relaxation is allowed for the one unmatched arm.
        patterns[-1] = (
            abbr,
            smarts_mol,
            r_positions,
            global_min_slot,
            arm_chem_types,
            max(2, len(r_positions) - 1),
            False,
        )

        # For thia_michael_c arms the SMARTS token CC[*:n] requires a heavy atom after
        # the beta-carbon, but pyPept's mol builder H-caps unconnected dummy atoms, turning
        # the arm into a terminal propanoyl (-CH2-CH3) with no such heavy neighbour.
        # The full SMARTS therefore cannot match any unreacted thia_michael_c arm.
        # Fix: generate partial SMARTS for every non-trivial proper subset of thia_michael_c
        # arm slots removed (1 up to N-1 slots removed).  Each partial SMARTS replaces the
        # removed [*:n] tokens with nothing, so the H-capped arm tail matches exactly.
        # min_required for partial patterns = number of remaining [*:n] slots (all must hit).
        import itertools as _it

        thia_slots = sorted(
            s for s, ct in arm_chem_types.items() if ct == "thia_michael_c"
        )
        if thia_slots:
            for _k in range(1, len(thia_slots)):  # remove 1…N-1 slots
                for _remove in _it.combinations(thia_slots, _k):
                    _partial_str = smarts_str
                    for _slot in _remove:
                        _partial_str = _re2.sub(
                            r"\[\*:" + str(_slot) + r"\]", "", _partial_str
                        )
                    _partial_mol = _C.MolFromSmarts(_partial_str)
                    if _partial_mol is None:
                        continue
                    _partial_r_pos = {
                        a.GetIdx(): a.GetAtomMapNum()
                        for a in _partial_mol.GetAtoms()
                        if a.GetAtomMapNum() >= 4
                    }
                    if len(_partial_r_pos) < 1:
                        continue
                    # is_partial=True: all remaining arms must hit backbone; no relaxation.
                    patterns.append(
                        (
                            abbr,
                            _partial_mol,
                            _partial_r_pos,
                            global_min_slot,
                            arm_chem_types,
                            len(_partial_r_pos),
                            True,
                        )
                    )

    with _s2c_lib_lock:
        _s2c_scaffold_cache = patterns
        _s2c_scaffold_stamp = stamp
    return patterns


def _s2c_extract_triazole_derx_ops(prod_mol, idx_a, idx_b):
    """Extract bond-surgery ops to de-react a 1,2,3-triazole back to alkyne + azide.

    Used as a pre-processing step before backbone detection so the coverage library
    can recognise alkyne residues (e.g. Pra) and azide residues (e.g. AzK) whose
    reactive groups have been consumed into the triazole ring.

    prod_mol : RDKit Mol built from the product SMARTS (MolFromSmarts).
    idx_a    : atom index in prod_mol of the alkyne-side attachment (e.g. Pra Cβ, map 1).
    idx_b    : atom index in prod_mol of the azide-side attachment (e.g. AzK Cε, map 5).

    Returns a dict with keys:
      'new_bond_pairs'  – list of (prod_idx_a, prod_idx_b) ring closure bonds to REMOVE
      'bond_restores'   – list of (prod_idx_a, prod_idx_b, BondType) pre-reaction bonds
      'atom_restores'   – dict prod_idx → {formal_charge, explicit_H, is_aromatic}
    or None if prod_mol does not contain a 1,2,3-triazole ring.
    """
    from rdkit import Chem as _C

    _C.FastFindRings(prod_mol)
    ring_info = prod_mol.GetRingInfo()

    # Locate the unique 1,2,3-triazole ring: 5-membered, 3 N + 2 C
    triazole_ring = None
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        atoms = [prod_mol.GetAtomWithIdx(i) for i in ring]
        if (
            sum(1 for a in atoms if a.GetAtomicNum() == 7) == 3
            and sum(1 for a in atoms if a.GetAtomicNum() == 6) == 2
        ):
            triazole_ring = set(ring)
            break
    if triazole_ring is None:
        return None

    # C4: the ring carbon bonded to the exocyclic alkyne attachment (idx_a, e.g. Pra Cβ)
    # C5: the other ring carbon (bonded to N1)
    c4_idx = c5_idx = None
    for ai in triazole_ring:
        atom = prod_mol.GetAtomWithIdx(ai)
        if atom.GetAtomicNum() != 6:
            continue
        nb_idxs = {nb.GetIdx() for nb in atom.GetNeighbors()}
        if idx_a in nb_idxs:
            c4_idx = ai
        else:
            c5_idx = ai
    if c4_idx is None or c5_idx is None:
        return None

    # N1: ring nitrogen bonded to the azide attachment (idx_b, e.g. AzK Cε)
    n1_idx = next(
        (
            nb.GetIdx()
            for nb in prod_mol.GetAtomWithIdx(idx_b).GetNeighbors()
            if nb.GetIdx() in triazole_ring and nb.GetAtomicNum() == 7
        ),
        None,
    )
    # N3: ring nitrogen bonded to C4
    n3_idx = next(
        (
            nb.GetIdx()
            for nb in prod_mol.GetAtomWithIdx(c4_idx).GetNeighbors()
            if nb.GetIdx() in triazole_ring and nb.GetAtomicNum() == 7
        ),
        None,
    )
    if n1_idx is None or n3_idx is None:
        return None
    ring_n_idxs = [
        ai for ai in triazole_ring if prod_mol.GetAtomWithIdx(ai).GetAtomicNum() == 7
    ]
    n2_idx = next((ai for ai in ring_n_idxs if ai != n1_idx and ai != n3_idx), None)
    if n2_idx is None:
        return None

    BT = _C.BondType
    return {
        # These two bonds are the ring-closure bonds formed in the cycloaddition
        "new_bond_pairs": [(c5_idx, n1_idx), (c4_idx, n3_idx)],
        # Restore pre-reaction bond types (were aromatic in triazole, now restored)
        "bond_restores": [
            (c4_idx, c5_idx, BT.TRIPLE),  # alkyne triple bond
            (n1_idx, n2_idx, BT.DOUBLE),  # N1=N2 in azide
            (n2_idx, n3_idx, BT.DOUBLE),  # N2=N3 in azide
        ],
        # Restore formal charges, H counts, and aromaticity for ring atoms
        "atom_restores": {
            n1_idx: {"formal_charge": 0, "explicit_H": 0, "is_aromatic": False},
            n2_idx: {"formal_charge": 1, "explicit_H": 0, "is_aromatic": False},
            n3_idx: {"formal_charge": -1, "explicit_H": 0, "is_aromatic": False},
            c4_idx: {"formal_charge": 0, "explicit_H": 0, "is_aromatic": False},
            c5_idx: {"formal_charge": 0, "explicit_H": 1, "is_aromatic": False},
        },
    }


def _s2c_pre_detect_ring_xlinks(mol, patterns):
    """Pre-detect ring crosslinks and build a de-reacted mol for backbone detection.

    For each ring crosslink pattern with de-reaction ops (e.g. CuAAC triazole),
    applies bond surgery to restore the pre-reaction functional groups (alkyne,
    azide) so the coverage library can correctly identify all backbone residues.

    Since only bond types / formal charges are changed (no atoms added or removed),
    the de-reacted mol has identical atom indices to the original mol.

    Returns (xlink_pairs, derx_mol) where:
      xlink_pairs : list of (matom_a, matom_b, slot_a, slot_b) with original atom indices
      derx_mol    : mol with ring bonds de-reacted (same atom count/indices as mol)
                    Falls back to mol if de-reaction is not possible / fails.
    """
    from rdkit import Chem as _C

    xlink_pairs = []
    rw = None  # lazy-init RWMol

    for pat in patterns:
        rxn_id, prod_mol, idx_a, idx_b, slot_a, slot_b = pat[:6]
        derx_ops = pat[6] if len(pat) > 6 else None
        if derx_ops is None:
            continue

        for rmatch in mol.GetSubstructMatches(prod_mol, useChirality=False):
            matom_a = rmatch[idx_a]
            matom_b = rmatch[idx_b]
            xlink_pairs.append((matom_a, matom_b, slot_a, slot_b))

            if rw is None:
                rw = _C.RWMol(mol)

            # 1. Remove new ring-closure bonds (cross-reactant bonds formed in reaction)
            for pa, pb in derx_ops.get("new_bond_pairs", []):
                oa, ob = rmatch[pa], rmatch[pb]
                if rw.GetBondBetweenAtoms(oa, ob) is not None:
                    rw.RemoveBond(oa, ob)

            # 2. Restore pre-reaction bond types (e.g. aromatic → triple/double)
            for pa, pb, bt in derx_ops.get("bond_restores", []):
                oa, ob = rmatch[pa], rmatch[pb]
                bd = rw.GetBondBetweenAtoms(oa, ob)
                if bd is not None:
                    bd.SetBondType(bt)
                    bd.SetIsAromatic(False)

            # 3. Restore atom formal charges, H counts, and aromaticity
            for pa, props in derx_ops.get("atom_restores", {}).items():
                oa = rmatch[pa]
                atom = rw.GetAtomWithIdx(oa)
                atom.SetIsAromatic(props.get("is_aromatic", False))
                atom.SetFormalCharge(props.get("formal_charge", 0))
                atom.SetNumExplicitHs(props.get("explicit_H", 0))
                atom.SetNoImplicit(True)

            break  # one match per pattern type

    if rw is None or not xlink_pairs:
        return xlink_pairs, mol

    # Sanitize the de-reacted mol (no atom addition/removal, just bond/charge changes)
    try:
        _C.SanitizeMol(rw)
        return xlink_pairs, rw.GetMol()
    except (ValueError, RuntimeError):
        try:
            _C.SanitizeMol(
                rw, _C.SanitizeFlags.SANITIZE_ALL ^ _C.SanitizeFlags.SANITIZE_KEKULIZE
            )
            return xlink_pairs, rw.GetMol()
        except (ValueError, RuntimeError):
            return xlink_pairs, mol  # fallback: original mol


def _s2c_get_ring_crosslink_patterns():
    """Return patterns for ring-forming crosslinks derived from reactions.yaml.

    For each reaction whose product places >1 atom between the two sidechain
    attachment points (e.g. CuAAC 1,2,3-triazole, SPAAC, IEDDA), extract a
    product SMARTS mol suitable for substructure matching, plus the atom indices
    of the two attachment points within that mol and the R-group slot numbers.

    Returns list of (rxn_id, prod_mol, attach_idx_a, attach_idx_b, slot_a, slot_b, derx_ops)
    where derx_ops is a dict of bond-surgery instructions (or None if not de-reactable).
    Cached for the process lifetime (derived from the static reactions.yaml).
    """
    global _s2c_ring_xlink_cache
    with _s2c_lib_lock:
        if _s2c_ring_xlink_cache is not None:
            return _s2c_ring_xlink_cache

    from rdkit import Chem as _C
    from rdkit.Chem import AllChem as _AC
    from rdkit.Chem import rdmolops as _rmo

    from pyPept.interfaces.reaction_library import REACTIONS as _RXNS

    patterns = []
    for rxn_id, entry in _RXNS.items():
        steps = entry.get("steps", [])
        if not steps:
            continue
        slot_a = entry.get("slot_a", 4)
        slot_b = entry.get("slot_b", 4)

        smirks = steps[-1]
        if " >> " not in smirks:
            continue
        try:
            rxn = _AC.ReactionFromSmarts(smirks)
        except (ValueError, RuntimeError):
            continue
        if rxn is None or rxn.GetNumProductTemplates() == 0:
            continue

        # Find which mapped atoms are directly bonded to dummy [n*] in reactants.
        # Those are the sidechain attachment points that trace back to backbone residues.
        attach_mapnums = []
        for ri in range(rxn.GetNumReactantTemplates()):
            tmpl = rxn.GetReactantTemplate(ri)
            for atom in tmpl.GetAtoms():
                if atom.GetAtomicNum() == 0 and atom.GetIsotope() > 0:
                    for nb in atom.GetNeighbors():
                        mn = nb.GetAtomMapNum()
                        if mn > 0:
                            attach_mapnums.append(mn)
                            break
        if len(attach_mapnums) != 2:
            continue

        mapnum_a, mapnum_b = attach_mapnums

        # Convert product template to a SMARTS query mol.
        prod_tmpl = rxn.GetProductTemplate(0)
        try:
            prod_smarts = _C.MolToSmarts(prod_tmpl)
            prod_mol = _C.MolFromSmarts(prod_smarts)
        except (ValueError, RuntimeError):
            continue
        if prod_mol is None:
            continue

        # Locate attachment point atom indices in the product mol.
        idx_a = next(
            (a.GetIdx() for a in prod_mol.GetAtoms() if a.GetAtomMapNum() == mapnum_a),
            None,
        )
        idx_b = next(
            (a.GetIdx() for a in prod_mol.GetAtoms() if a.GetAtomMapNum() == mapnum_b),
            None,
        )
        if idx_a is None or idx_b is None:
            continue

        # Only keep ring-forming crosslinks: attachment points separated by ≥3 bonds
        # (shortest path through the product ≥ 4 atoms).
        # - Single-bond crosslinks (S-S, S-C): path length 2 → excluded
        # - 2-bond links (NHS amide N-CO-C, oxime, hydrazone): path length 3 → excluded
        # - Triazole (CuAAC/SPAAC): path goes through 3 ring atoms → length ≥ 4 → kept
        try:
            path = _rmo.GetShortestPath(prod_mol, idx_a, idx_b)
        except (ValueError, RuntimeError):
            continue
        if len(path) < 4:
            continue

        derx_ops = _s2c_extract_triazole_derx_ops(prod_mol, idx_a, idx_b)
        patterns.append((rxn_id, prod_mol, idx_a, idx_b, slot_a, slot_b, derx_ops))

    with _s2c_lib_lock:
        _s2c_ring_xlink_cache = patterns
    return patterns


def _s2c_r_group_at_atom(orig_mol, orig_atom_idx, frag_atoms, abbr, raw_lib):
    """Find the R-group number (1-6) for orig_atom_idx (a junction atom) in monomer abbr.

    Matches the fragment's core atom graph against the raw SDF mol (with dummy R-atoms),
    then returns the isotope label of the dummy atom adjacent to the matched position.
    Returns None if no match found.
    """
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol as _RW

    raw_mol = raw_lib.get(abbr)
    if raw_mol is None:
        return None

    # Collect: for each heavy atom in raw_mol, note adjacent dummy isotopes
    dummy_nbrs = {}
    for a in raw_mol.GetAtoms():
        if a.GetAtomicNum() == 0 and a.GetIsotope() > 0:
            for nb in a.GetNeighbors():
                dummy_nbrs.setdefault(nb.GetIdx(), []).append(a.GetIsotope())

    # Build core fragment mol from frag_atoms (bonds within fragment only)
    rw = _RW()
    old_to_new = {}
    for oi in sorted(frag_atoms):
        old_to_new[oi] = rw.AddAtom(_C.Atom(orig_mol.GetAtomWithIdx(oi).GetAtomicNum()))
    for bond in orig_mol.GetBonds():
        bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in old_to_new and ei in old_to_new:
            rw.AddBond(old_to_new[bi], old_to_new[ei], bond.GetBondType())
    core_frag = rw.GetMol()
    new_junct_idx = old_to_new.get(orig_atom_idx)
    if new_junct_idx is None:
        return None

    # Build stripped raw mol (remove dummy atoms) and track index mapping
    raw_to_stripped = {}
    stripped_to_raw = {}
    si = 0
    for ri in range(raw_mol.GetNumAtoms()):
        if raw_mol.GetAtomWithIdx(ri).GetAtomicNum() != 0:
            raw_to_stripped[ri] = si
            stripped_to_raw[si] = ri
            si += 1
    rw2 = _RW(raw_mol)
    for di in sorted(
        [a.GetIdx() for a in raw_mol.GetAtoms() if a.GetAtomicNum() == 0], reverse=True
    ):
        rw2.RemoveAtom(di)
    stripped = rw2.GetMol()

    # Match core_frag within stripped (find core_frag as sub-graph of stripped)
    matches = stripped.GetSubstructMatches(core_frag, useChirality=False)
    for match in matches:
        if new_junct_idx < len(match):
            stripped_idx = match[new_junct_idx]
            raw_idx = stripped_to_raw.get(stripped_idx)
            if raw_idx is None:
                continue
            isos = dummy_nbrs.get(raw_idx, [])
            if isos:
                return min(isos)

    # Reverse: find stripped within core_frag.
    # Needed when the assembled segment has extra free-cap atoms (e.g. unused R2 α-OH stays
    # when R4 is the chain continuation) — core_frag > stripped so forward match fails.
    rev_matches = core_frag.GetSubstructMatches(stripped, useChirality=False)
    for rev_match in rev_matches:
        for stripped_i, core_i in enumerate(rev_match):
            if core_i == new_junct_idx:
                raw_idx = stripped_to_raw.get(stripped_i)
                if raw_idx is None:
                    continue
                isos = dummy_nbrs.get(raw_idx, [])
                if isos:
                    return min(isos)
    return None


_BRANCH_SEG_SKIP: frozenset = frozenset()  # no deprecated aliases remain


def _s2c_match_branch_segment(seg_frag, full_lib):
    """Match a branch segment mol against the full library; return (abbr, n_matched).

    Tries two directions:
    1. Forward: lib as subgraph of segment (segment >= lib in size).
    2. Reverse: segment as subgraph of lib (lib has one or two extra cap atoms like -OH
       that were consumed when forming the amide bond in the full peptide).

    R-group slot is determined by _s2c_r_group_at_atom from connectivity.
    """
    seg_norm = _s2c_normalize(seg_frag)
    if seg_norm is None:
        return None, 0
    best_abbr = None
    best_n = 0
    best_lib_n = float(
        "inf"
    )  # prefer smallest library entry for same n_m (tightest fit)
    seg_n = seg_norm.GetNumAtoms()
    for abbr, cm, norm_cm, n_lib in full_lib:
        if abbr in _BRANCH_SEG_SKIP:
            continue
        n_m = 0
        # Forward: find lib monomer as subgraph of our segment
        if n_lib <= seg_n + 2:
            matches = seg_norm.GetSubstructMatches(norm_cm, useChirality=False)
            if matches:
                n_m = len(matches[0])
        # Reverse: find our segment as subgraph of lib (lib has extra cap -OH/-H atoms).
        # Only try reverse if forward failed — protected variants (e.g. Glu_OAll) share
        # the Glu core but carry extra atoms; tightest-fit tie-break below prefers base E.
        if n_m == 0 and seg_n < n_lib <= seg_n + 4:
            matches = norm_cm.GetSubstructMatches(seg_norm, useChirality=False)
            if matches:
                n_m = len(matches[0])  # = seg_n atoms matched within lib
        # Prefer: (1) higher n_m; (2) smaller n_lib for same n_m (tightest match).
        # This ensures E (n_lib=10) beats Glu_OAll (n_lib=13) when both match 9 atoms.
        if n_m > best_n or (n_m == best_n and n_m > 0 and n_lib < best_lib_n):
            best_n = n_m
            best_abbr = abbr
            best_lib_n = n_lib
            if best_n == seg_n and best_lib_n == seg_n:
                break  # perfect exact-size match — stop searching
    return best_abbr, best_n


def _s2c_walk_branch(
    mol,
    glu_atoms,
    glu_outgoing_n,
    orphan_atoms,
    full_lib,
    raw_lib,
    anchor_abbr,
    anchor_n_atom,
    junction_abbr="E_g",
):
    """Walk the branch chain starting from the junction monomer outward.

    Returns list of (abbr, prev_r, cur_r) for each piece in the branch
    (including the junction monomer).
    glu_outgoing_n: junction monomer's alpha-N atom index (outgoing)
    anchor_n_atom: K's epsilon-N atom index (main_at)
    junction_abbr: abbreviation of the junction monomer (detected, not hardcoded)
    """
    from collections import deque as _dq

    from rdkit import Chem as _C
    from rdkit.Chem import RWMol as _RW

    results = []

    # ── Step 1: R-group for K (anchor) → junction monomer connection ───────
    branch_at = None
    for nb in mol.GetAtomWithIdx(anchor_n_atom).GetNeighbors():
        if nb.GetIdx() in glu_atoms:
            branch_at = nb.GetIdx()
            break

    anchor_r = _s2c_r_group_at_atom(
        mol,
        anchor_n_atom,
        set(range(mol.GetNumAtoms())) - glu_atoms - orphan_atoms,
        anchor_abbr,
        raw_lib,
    )
    glu_in_r = (
        _s2c_r_group_at_atom(mol, branch_at, glu_atoms, junction_abbr, raw_lib)
        if branch_at is not None
        else None
    )

    # ── Step 2: Identify junction monomer's outgoing R-group ────────────────
    glu_out_r = _s2c_r_group_at_atom(
        mol, glu_outgoing_n, glu_atoms, junction_abbr, raw_lib
    )

    results.append((junction_abbr, anchor_r or 4, glu_in_r or 4))

    if glu_outgoing_n is None or not orphan_atoms:
        return results

    # ── Step 3: Segment orphan atoms at amide N→C=O bonds ───────────────────
    # Walk from the first orphan atom (bonded to glu_outgoing_n) through orphan_atoms
    # Split at each internal amide bond (orphan-N bonded to orphan-C=O)

    # Find orphan entry atom (bonded to glu_outgoing_n)
    entry_orphan = None
    for nb in mol.GetAtomWithIdx(glu_outgoing_n).GetNeighbors():
        if nb.GetIdx() in orphan_atoms:
            entry_orphan = nb.GetIdx()
            break
    if entry_orphan is None:
        return results

    # BFS to find all connected orphan atoms in chain order
    # Build segments by splitting at each amide-N → C=O bond within orphans
    segments = []
    current_seg = set()
    visited = set()

    def is_co(ai):
        atom = mol.GetAtomWithIdx(ai)
        if atom.GetSymbol() != "C":
            return False
        return any(
            nb.GetSymbol() == "O"
            and mol.GetBondBetweenAtoms(ai, nb.GetIdx()).GetBondTypeAsDouble() == 2.0
            for nb in atom.GetNeighbors()
        )

    # Walk the linear orphan chain
    queue = _dq([entry_orphan])
    junctions = []  # list of (N_atom, CO_atom) internal amide bonds

    while queue:
        ai = queue.popleft()
        if ai in visited:
            continue
        visited.add(ai)
        current_seg.add(ai)
        for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
            ni = nb.GetIdx()
            if ni in visited or ni not in orphan_atoms:
                continue
            # Check if this is an amide bond crossing (orphan-N to orphan-C=O or vice versa)
            atom_ai = mol.GetAtomWithIdx(ai)
            if atom_ai.GetSymbol() == "N" and is_co(ni):
                junctions.append((ai, ni))
            queue.append(ni)

    # Use junctions to split orphan_atoms into segments
    if not junctions:
        segments = [set(orphan_atoms)]
    else:
        # Rebuild segments by connectivity after removing junction bonds
        remaining = set(orphan_atoms)
        for n_atom, co_atom in junctions:
            # Split: everything reachable from entry_orphan without crossing n_atom→co_atom bond
            seg = set()
            q2 = _dq([entry_orphan])
            while q2:
                ai2 = q2.popleft()
                if ai2 in seg:
                    continue
                seg.add(ai2)
                for nb2 in mol.GetAtomWithIdx(ai2).GetNeighbors():
                    ni2 = nb2.GetIdx()
                    if ni2 in seg or ni2 not in remaining:
                        continue
                    if ai2 == n_atom and ni2 == co_atom:
                        continue  # skip the junction bond
                    q2.append(ni2)
            segments.append(seg)
            entry_orphan = co_atom
            remaining -= seg

        if remaining:
            segments.append(remaining)

    # ── Step 4: Match each segment and find R-groups ─────────────────────────
    prev_out_r = glu_out_r or 1  # E_g's outgoing R (R1 for backbone_n/in_n)
    seg_entry_atoms = []  # entry atom for each segment (bonded to previous piece)

    # Reconstruct entry atom for each segment
    seg_entry = glu_outgoing_n
    for seg_idx, seg in enumerate(segments):
        # Find which atom in seg is bonded to the previous piece's outgoing atom
        seg_in_atom = None
        for ai in seg:
            for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                if nb.GetIdx() == seg_entry or (
                    seg_idx > 0 and nb.GetIdx() in segments[seg_idx - 1]
                ):
                    seg_in_atom = ai
                    break
            if seg_in_atom is not None:
                break
        # Hmm, above is fragile. Better: find atom in seg bonded to atom NOT in seg
        # that is either glu_outgoing_n or in previous segment
        if seg_in_atom is None:
            for ai in seg:
                for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                    ni2 = nb.GetIdx()
                    if ni2 not in seg and ni2 in (
                        set()
                        if seg_idx == 0
                        else segments[seg_idx - 1] | {glu_outgoing_n}
                    ):
                        seg_in_atom = ai
                        break
                if seg_in_atom is not None:
                    break
        seg_entry_atoms.append(seg_in_atom)

    # Find entry atom for first segment (bonded to glu_outgoing_n)
    seg0_in = None
    for ai in segments[0]:
        for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
            if nb.GetIdx() == glu_outgoing_n:
                seg0_in = ai
                break
        if seg0_in is not None:
            break
    if seg0_in is not None:
        seg_entry_atoms[0] = seg0_in

    for seg_idx, seg in enumerate(segments):
        # Build segment mol
        rw3 = _RW()
        s2n = {}
        for ai in sorted(seg):
            _src = mol.GetAtomWithIdx(ai)
            _dst = _C.Atom(_src.GetAtomicNum())
            _dst.SetIsAromatic(_src.GetIsAromatic())
            _dst.SetFormalCharge(_src.GetFormalCharge())
            _dst.SetNumExplicitHs(_src.GetNumExplicitHs())
            _dst.SetNoImplicit(_src.GetNoImplicit())
            s2n[ai] = rw3.AddAtom(_dst)
        for bond in mol.GetBonds():
            bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bi in s2n and ei in s2n:
                rw3.AddBond(s2n[bi], s2n[ei], bond.GetBondType())
        seg_mol = _C.RemoveHs(rw3.GetMol())

        abbr, n_matched = _s2c_match_branch_segment(seg_mol, full_lib)

        # 1-2 rule: prefer the alternative whose ENTRY slot uses the lowest R-group.
        # E_g (R2=γ-COOH entry) beats E (R4=γ-COOH entry) for γ-linked chains.
        if abbr:
            _seg_in = (
                seg_entry_atoms[seg_idx] if seg_idx < len(seg_entry_atoms) else None
            )
            if _seg_in is not None:
                _cur_in_r = _s2c_r_group_at_atom(mol, _seg_in, seg, abbr, raw_lib)
                if _cur_in_r is not None and _cur_in_r > 2:
                    _seg_norm = _s2c_normalize(seg_mol)
                    _seg_n = _seg_norm.GetNumAtoms() if _seg_norm else 0
                    for _alt_abbr, _alt_cm, _alt_norm, _alt_n in full_lib:
                        if _alt_abbr == abbr or _alt_abbr in _BRANCH_SEG_SKIP:
                            continue
                        _alt_nm = 0
                        if _seg_norm and _alt_n <= _seg_n + 2:
                            _m = _seg_norm.GetSubstructMatches(
                                _alt_norm, useChirality=False
                            )
                            if _m:
                                _alt_nm = len(_m[0])
                        if _alt_nm == 0 and _seg_n < _alt_n <= _seg_n + 4:
                            _m = _alt_norm.GetSubstructMatches(
                                _seg_norm, useChirality=False
                            )
                            if _m:
                                _alt_nm = _seg_n
                        if _alt_nm < n_matched:
                            continue
                        _alt_in_r = _s2c_r_group_at_atom(
                            mol, _seg_in, seg, _alt_abbr, raw_lib
                        )
                        if _alt_in_r is not None and _alt_in_r < _cur_in_r:
                            abbr = _alt_abbr
                            break

        seg_in_atom = (
            seg_entry_atoms[seg_idx] if seg_idx < len(seg_entry_atoms) else None
        )

        # cur_r: R-group on this segment at incoming junction
        cur_r = None
        if seg_in_atom is not None and abbr:
            cur_r = _s2c_r_group_at_atom(mol, seg_in_atom, seg, abbr, raw_lib)

        # Find outgoing junction atom of this segment (bonded to next segment or terminus)
        seg_out_atom = None
        if seg_idx < len(segments) - 1:
            next_seg = segments[seg_idx + 1]
            for ai in seg:
                for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                    if nb.GetIdx() in next_seg:
                        seg_out_atom = ai
                        break
                if seg_out_atom is not None:
                    break

        # Format: (abbr, prev_out_R, cur_R) — AEEA(1,2) means E_g's R1 connects
        # to AEEA's R2; prev_out_R=1 (E_g.in_n), cur_R=2 (AEEA.out_co).
        results.append((abbr or "?", prev_out_r, cur_r or 2))

        if seg_out_atom is not None and abbr:
            prev_out_r = (
                _s2c_r_group_at_atom(mol, seg_out_atom, seg, abbr, raw_lib) or 1
            )
        else:
            prev_out_r = 1

    return results


def _s2c_cip_score(aa_mol, ref_mol, match_tuple):
    """Score CIP stereo agreement for a match.

    Returns (n_agree, n_ref_stereo) where n_agree counts stereocenters with
    matching CIP codes and n_ref_stereo counts defined stereocenters in ref.
    Disagreements score 0 (not -1) because CIP codes are context-dependent:
    the same absolute configuration can produce different R/S labels when the
    molecular graph changes (e.g. isolated fragment vs full peptide chain).
    """
    from rdkit.Chem import AllChem as _AC

    _AC.AssignStereochemistry(aa_mol, cleanIt=True, force=True)
    _AC.AssignStereochemistry(ref_mol, cleanIt=True, force=True)
    n_agree = 0
    n_ref_stereo = 0
    for ref_idx, aa_idx in enumerate(match_tuple):
        ref_cip = ref_mol.GetAtomWithIdx(ref_idx).GetPropsAsDict().get("_CIPCode")
        aa_cip = aa_mol.GetAtomWithIdx(aa_idx).GetPropsAsDict().get("_CIPCode")
        if ref_cip:
            n_ref_stereo += 1
            if aa_cip and ref_cip == aa_cip:
                n_agree += 1
    return n_agree, n_ref_stereo


def _s2c_match(aa_mol, lib):
    """Match an isolated residue mol against the library.

    Priority: (1) most atoms matched, (2) ring-count agreement (prevents linear
    patterns false-matching cyclic residues like Pip/Aze), (3) CIP stereo score,
    (4) E/Z double-bond stereo score.
    Early break only when all four are perfect (full CIP + E/Z agreement).
    """
    from rdkit.Chem import rdchem as _RC
    from rdkit.Chem import rdMolDescriptors as _RDM

    norm_mol = _s2c_normalize(aa_mol)
    n_q = aa_mol.GetNumAtoms()
    n_rings_aa = _RDM.CalcNumRings(aa_mol)

    # Save E/Z stereo BEFORE _s2c_cip_score (AssignStereochemistry cleanIt=True) corrupts it.
    # Keyed by frozenset({bi, ei}) → stereo enum; used for E/Z scoring in the loop.
    _BS = _RC.BondStereo
    _ez_saved = {}
    for b in aa_mol.GetBonds():
        if b.GetBondTypeAsDouble() == 2.0:
            s = b.GetStereo()
            if s not in (_BS.STEREONONE, _BS.STEREOANY):
                _ez_saved[frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))] = s
    query_has_ez = bool(_ez_saved)
    n_query_stereo = sum(a.HasProp("_CIPCode") for a in aa_mol.GetAtoms())

    # D-form tiebreaker: when CIP scores are equal (e.g. flat/achiral input),
    # prefer L-form entries so T beats dT, A beats dA, etc.
    import re as _re2

    def _is_d_form(tok):
        return bool(_re2.match(r"^(D_|[Dd][A-Z])", tok))

    best_abbr = None
    best_atom_n = 0
    best_stereo = -999
    best_n_ref_stereo = 0
    best_ring_match = False
    best_ez = -2
    best_is_d = True  # so the first L-form candidate always beats the initial state
    best_n_rgroups = -1

    for (
        abbr,
        ref_orig,
        ref_norm,
        n_ref,
        n_rings_ref,
        ref_has_ez,
        _bb_dist,
        n_rgroups_extra,
        *_,
    ) in lib:
        if n_ref > n_q:
            continue
        match = aa_mol.GetSubstructMatches(ref_orig, useChirality=False)
        if match:
            n_m = len(match[0])
            pairs = [_s2c_cip_score(aa_mol, ref_orig, m) for m in match]
            sc = max(a for a, _ in pairs)
            n_rs = max(r for _, r in pairs)
        else:
            match = norm_mol.GetSubstructMatches(ref_norm, useChirality=False)
            if not match:
                continue
            n_m = len(match[0])
            sc = 0
            n_rs = 0

        rings_match = n_rings_ref == n_rings_aa

        # E/Z scoring using saved pre-CIP state: +1 match, 0 neutral, -1 mismatch
        if query_has_ez and ref_has_ez:
            ez_sc = 0
            for b in ref_orig.GetBonds():
                if b.GetBondTypeAsDouble() != 2.0:
                    continue
                ref_s = b.GetStereo()
                if ref_s in (_BS.STEREONONE, _BS.STEREOANY):
                    continue
                ri, rj = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                qi = match[0][ri] if ri < len(match[0]) else None
                qj = match[0][rj] if rj < len(match[0]) else None
                if qi is None or qj is None:
                    continue
                saved = _ez_saved.get(frozenset((qi, qj)))
                if saved is None:
                    ez_sc = 0
                    break
                ez_sc = 1 if saved == ref_s else -1
        else:
            ez_sc = 0

        is_d = _is_d_form(abbr)
        # Tiebreaker: when all stereo+ring+L/D criteria are equal, prefer the
        # monomer with more declared R-group slots (e.g. XlQuat over Ile when
        # both match the chondramide quat-C residue — XlQuat declares R4 for
        # the aryl-ether crosslink, Ile doesn't, so XlQuat's CABILN can
        # actually round-trip).
        prior_eq = (
            n_m == best_atom_n
            and rings_match == best_ring_match
            and sc == best_stereo
            and n_rs == best_n_ref_stereo
            and ez_sc == best_ez
            and is_d == best_is_d
        )
        better = (
            n_m > best_atom_n
            or (n_m == best_atom_n and rings_match and not best_ring_match)
            or (
                n_m == best_atom_n
                and rings_match == best_ring_match
                and sc > best_stereo
            )
            or (
                n_m == best_atom_n
                and rings_match == best_ring_match
                and sc == best_stereo
                and n_rs > best_n_ref_stereo
            )
            or (
                n_m == best_atom_n
                and rings_match == best_ring_match
                and sc == best_stereo
                and n_rs == best_n_ref_stereo
                and ez_sc > best_ez
            )
            or (
                n_m == best_atom_n
                and rings_match == best_ring_match
                and sc == best_stereo
                and n_rs == best_n_ref_stereo
                and ez_sc == best_ez
                and not is_d
                and best_is_d
            )
            or (prior_eq and n_rgroups_extra > best_n_rgroups)
        )
        if better:
            best_atom_n = n_m
            best_abbr = abbr
            best_stereo = sc
            best_n_ref_stereo = n_rs
            best_ring_match = rings_match
            best_ez = ez_sc
            best_is_d = is_d
            best_n_rgroups = n_rgroups_extra
            if n_m == n_q and rings_match:
                if sc == n_query_stereo and (not query_has_ez or ez_sc > 0):
                    break
    return best_abbr, best_atom_n


# ── Coverage library: backbone monomers with tracked in_n / out_co ────────────
_CovEntry = _nt2(
    "_CovEntry",
    [
        "abbr",
        "m_type",
        "orig_mol",
        "norm_mol",
        "smarts_mol",
        "in_n_idx",
        "out_co_idx",  # indices into orig_mol / norm_mol
        "smarts_in_n_idx",
        "smarts_out_co_idx",  # indices into smarts_mol
        "n_atoms",
        "n_rings",
        "has_ez",
        "bb_dist",  # shortest-path bond count N→backbone-C
        # iso SMARTS: carboxyl cap O removed so isopeptide-bonded monomers (E_g→K)
        # match without pulling the backbone N into the match set.  None if no
        # carboxyl R-groups are present.
        "smarts_mol_iso",
        "smarts_in_n_idx_iso",
        "smarts_out_co_idx_iso",
    ],
)
_PlacedNode = _nt2(
    "_PlacedNode",
    [
        "abbr",
        "m_type",
        "mol_atoms",
        "in_n",
        "out_co",
        "entry",
        "match_tuple",
    ],
)

_s2c_cov_lib_cache: list = []
_s2c_cov_lib_stamp: tuple | None = None


def _s2c_get_coverage_lib():
    """Backbone monomer library; each entry tracks in_n / out_co atom indices."""
    global _s2c_cov_lib_cache, _s2c_cov_lib_stamp
    stamp = _s2c_library_stamp()
    with _s2c_lib_lock:
        if _s2c_cov_lib_cache and _s2c_cov_lib_stamp == stamp:
            return _s2c_cov_lib_cache
        lib = _s2c_build_cov_lib()
        _s2c_cov_lib_cache = lib
        _s2c_cov_lib_stamp = stamp
        return lib


def _s2c_build_cov_lib():
    """Build coverage library: backbone monomers with atom-map-tracked terminals."""
    from rdkit import Chem as _C
    from rdkit.Chem import RWMol as _RW
    from rdkit.Chem import rdchem as _RC
    from rdkit.Chem import rdMolDescriptors as _RDM

    _MAP_IN, _MAP_OUT = 9901, 9902
    raw_mols, _ = _load_sdf()
    lib = []
    seen: set = set()
    for raw_mol in raw_mols:
        props = raw_mol.GetPropsAsDict()
        abbr = props.get("m_abbr", "").strip()
        if not abbr or abbr in seen:
            continue
        m_type = props.get("m_type", "aa")
        ct_str = props.get("m_chem_types", "")
        rg_str = props.get("m_Rgroups", "")
        cts: dict = {}
        for item in ct_str.split(","):
            item = item.strip()
            if ":" in item:
                k, v = item.split(":", 1)
                try:
                    cts[int(k)] = v.strip()
                except ValueError:
                    pass
        bn_r = next(
            (r for r, t in cts.items() if t in ("backbone_n", "backbone_o")), None
        )
        bc_r = next(
            (
                r
                for r, t in cts.items()
                if t in ("backbone_c", "backbone_c_red", "sp3_c_anchor")
            ),
            None,
        )
        if bn_r is None or bc_r is None:
            continue
        bn_dummy = next(
            (
                a.GetIdx()
                for a in raw_mol.GetAtoms()
                if a.GetAtomicNum() == 0 and a.GetIsotope() == bn_r
            ),
            None,
        )
        bc_dummy = next(
            (
                a.GetIdx()
                for a in raw_mol.GetAtoms()
                if a.GetAtomicNum() == 0 and a.GetIsotope() == bc_r
            ),
            None,
        )
        if bn_dummy is None or bc_dummy is None:
            continue
        bn_nbrs = [
            nb.GetIdx() for nb in raw_mol.GetAtomWithIdx(bn_dummy).GetNeighbors()
        ]
        bc_nbrs = [
            nb.GetIdx() for nb in raw_mol.GetAtomWithIdx(bc_dummy).GetNeighbors()
        ]
        if not bn_nbrs or not bc_nbrs:
            continue
        in_n_raw = bn_nbrs[0]
        out_co_raw = bc_nbrs[0]
        # Relax H-count SMARTS for atoms adjacent to R-groups that undergo
        # connectivity changes during crosslinking (terminal_alkene → =CH2→=CH-).
        # Restrict to alkene-type R-groups ONLY — other non-backbone R-groups
        # like formamide_c/amide_nh must keep their exact H count or they absorb
        # neighbouring residue fragments (e.g. Lys_For absorbing E_g gamma-CO).
        _MAP_XL = 9904
        # Chem-types whose adjacent host carbon undergoes H-count change on
        # crosslink formation. Listed atoms get [#6] (no H pin) in SMARTS so
        # the unbonded monomer and the bonded peptide both match.
        #   terminal_alkene / crosslink_alkene: =CH2 → =CH-
        #   quat_c_anchor:   sp3 CHn → sp3 C-X  (chondramide aryl-ether bridge)
        _XL_ACTIVE_TYPES = frozenset(
            {
                "terminal_alkene",
                "crosslink_alkene",
                "quat_c_anchor",
                # Reduced-amide C-side: cap H becomes the next residue's N-bond,
                # so the host CH3-cap → CH2-bonded transition needs no H pin.
                "backbone_c_red",
            }
        )
        xl_raw_attach: set = set()
        for _rl_r, _rl_t in cts.items():
            if _rl_t not in _XL_ACTIVE_TYPES:
                continue
            _xl_dum = next(
                (
                    a.GetIdx()
                    for a in raw_mol.GetAtoms()
                    if a.GetAtomicNum() == 0 and a.GetIsotope() == _rl_r
                ),
                None,
            )
            if _xl_dum is not None:
                for _xl_nb in raw_mol.GetAtomWithIdx(_xl_dum).GetNeighbors():
                    _xl_idx = _xl_nb.GetIdx()
                    if _xl_idx in (in_n_raw, out_co_raw):
                        continue
                    xl_raw_attach.add(_xl_idx)
                    for _xl_bd in raw_mol.GetAtomWithIdx(_xl_idx).GetBonds():
                        if _xl_bd.GetBondTypeAsDouble() == 2.0:
                            _pi = (
                                _xl_bd.GetEndAtomIdx()
                                if _xl_bd.GetBeginAtomIdx() == _xl_idx
                                else _xl_bd.GetBeginAtomIdx()
                            )
                            if _pi != _xl_dum and _pi not in (in_n_raw, out_co_raw):
                                xl_raw_attach.add(_pi)
        # For carboxyl-type R-groups (isopeptide bond forming), the R-group cap
        # atom (O from [OH] capping) ends up in the SMARTS as [#8] but in the
        # assembled molecule that position is the amide-N of the isopeptide bond.
        # Track attachment atoms with _MAP_CR so we can find and remove the
        # capping O from the SMARTS (same treatment as bc_r_cap_idx for out_co).
        _MAP_CR = 9906
        _CR_ACTIVE_TYPES = frozenset({"carboxyl"})
        cap_remove_raw_attach: set = set()
        for _rl_r, _rl_t in cts.items():
            if _rl_t not in _CR_ACTIVE_TYPES:
                continue
            _cr_dum = next(
                (
                    a.GetIdx()
                    for a in raw_mol.GetAtoms()
                    if a.GetAtomicNum() == 0 and a.GetIsotope() == _rl_r
                ),
                None,
            )
            if _cr_dum is not None:
                for _cr_nb in raw_mol.GetAtomWithIdx(_cr_dum).GetNeighbors():
                    _cr_idx = _cr_nb.GetIdx()
                    if _cr_idx not in (in_n_raw, out_co_raw):
                        cap_remove_raw_attach.add(_cr_idx)
        rw = _RW(raw_mol)
        rw.GetAtomWithIdx(in_n_raw).SetAtomMapNum(_MAP_IN)
        rw.GetAtomWithIdx(out_co_raw).SetAtomMapNum(_MAP_OUT)
        for _xl_raw in xl_raw_attach:
            rw.GetAtomWithIdx(_xl_raw).SetAtomMapNum(_MAP_XL)
        for _cr_raw in cap_remove_raw_attach:
            rw.GetAtomWithIdx(_cr_raw).SetAtomMapNum(_MAP_CR)
        smi = _C.MolToSmiles(rw.GetMol(), allHsExplicit=False)
        # Cap bn_r with [H] so the query matches in any peptide context.
        # Keep bc_r's SDF default (e.g. [OH] for standard AAs, [H] for _al forms)
        # so that Gly vs Gly_al can be distinguished by bc_r cap presence.
        rg_parts = [rg.strip() for rg in rg_str.split(",")]
        for _ri in range(len(rg_parts)):
            if (_ri + 1) == bn_r:
                rg_parts[_ri] = "[H]"
        cov_rg_str = ", ".join(rg_parts)
        capped_smi = _s2c_cap_smiles(smi, cov_rg_str)
        try:
            cm = _C.MolFromSmiles(capped_smi)
        except (ValueError, RuntimeError):
            continue
        if cm is None:
            continue
        cm = _C.RemoveHs(cm)
        in_n_idx = next(
            (a.GetIdx() for a in cm.GetAtoms() if a.GetAtomMapNum() == _MAP_IN), None
        )
        out_co_idx = next(
            (a.GetIdx() for a in cm.GetAtoms() if a.GetAtomMapNum() == _MAP_OUT), None
        )
        if in_n_idx is None or out_co_idx is None:
            continue
        xl_cm_idxs = frozenset(
            a.GetIdx() for a in cm.GetAtoms() if a.GetAtomMapNum() == _MAP_XL
        )
        cap_remove_cm_idxs = frozenset(
            a.GetIdx() for a in cm.GetAtoms() if a.GetAtomMapNum() == _MAP_CR
        )
        # Degree-1 single-bond neighbors of carboxyl-attachment atoms are the
        # O atoms added when the carboxyl R-group dummy was replaced by [OH].
        # In the assembled molecule that position is the isopeptide amide-N, so
        # remove these capping O atoms from the SMARTS.
        cap_remove_cap_set: set = set()
        for _cr_cm in cap_remove_cm_idxs:
            for _cr_nb in cm.GetAtomWithIdx(_cr_cm).GetNeighbors():
                _bd = cm.GetBondBetweenAtoms(_cr_cm, _cr_nb.GetIdx())
                if _cr_nb.GetDegree() == 1 and _bd.GetBondTypeAsDouble() != 2.0:
                    cap_remove_cap_set.add(_cr_nb.GetIdx())
        # Build H-count-specific SMARTS while cm still carries _MAP_IN/_MAP_OUT.
        # - in_n → [#7:MAP] (any N; matches free, protonated, or amide N-terminus)
        # - out_co → [#{anum}H{nh}:MAP] using ACTUAL H count from capped mol.
        #   Standard AAs (bc_r=[OH] default) → 0H = distinguishes them from _al forms.
        #   Aldehyde forms (bc_r=[H] default → removed by RemoveHs) → 1H = CHO.
        # - bc_r cap atom (degree-1 neighbor of out_co, not =O) → marked 9903,
        #   then removed from smarts_mol so it doesn't constrain C-terminal context.
        # - All other side-chain atoms → [#{anum}H{nh}] (exact H count) to prevent
        #   neighbouring cap atoms from being absorbed into the match.
        import re as _re2

        _MAP_CAP = 9903
        bc_r_cap_idx = None
        for _nb in cm.GetAtomWithIdx(out_co_idx).GetNeighbors():
            _bd = cm.GetBondBetweenAtoms(out_co_idx, _nb.GetIdx())
            if _bd.GetBondTypeAsDouble() != 2.0 and _nb.GetDegree() == 1:
                bc_r_cap_idx = _nb.GetIdx()
                break
        # Use a per-atom map number to track carboxyl-cap O atoms so we can
        # replace their SMARTS with [#8,#7] (O or N) rather than [#8] alone.
        # This lets E_g match in both free-carboxyl (-OH) and isopeptide (-N)
        # forms while still requiring D/E's carboxyl in free-acid contexts.
        _MAP_CR_CAP = 9907
        rw_s = _RW(cm)
        idx_to_mn: dict = {}
        for _sa in rw_s.GetAtoms():
            _si = _sa.GetIdx()
            if _sa.GetAtomMapNum() in (_MAP_IN, _MAP_OUT):
                idx_to_mn[_si] = _sa.GetAtomMapNum()  # keep 9901/9902 for tracking
            elif _si == bc_r_cap_idx:
                _sa.SetAtomMapNum(_MAP_CAP)
                idx_to_mn[_si] = _MAP_CAP
            elif _si in cap_remove_cap_set:
                _sa.SetAtomMapNum(_MAP_CR_CAP)
                idx_to_mn[_si] = _MAP_CR_CAP
            else:
                _mn = _si + 1  # 1-indexed; won't clash with reserved map nums
                _sa.SetAtomMapNum(_mn)
                idx_to_mn[_si] = _mn
        _smarts_str = _C.MolToSmarts(rw_s.GetMol())
        for _sa in cm.GetAtoms():
            _si = _sa.GetIdx()
            _mn = idx_to_mn[_si]
            _nh = _sa.GetTotalNumHs()
            _ch = _sa.GetFormalCharge()
            _cs = f"+{_ch}" if _ch > 0 else (str(_ch) if _ch < 0 else "")
            if _si == in_n_idx:
                # Reduced-amide N-side: the N must be a sp3 secondary amine
                # (NX3 with at least 1 H in the free monomer, or 0 H when bonded
                # to next residue). Exclude amide N (NX3 bonded to C=O), aromatic
                # n, and iminol/sp2 =N — those belong to other residue types.
                if cts.get(bc_r) == "backbone_c_red":
                    _repl = f"[NX3;!$(N-C=O);!$(N=*):{_MAP_IN}]"
                else:
                    _repl = f"[#{_sa.GetAtomicNum()}:{_MAP_IN}]"  # permissive
            elif _si == out_co_idx:
                # Reduced-amide C-side (backbone_c_red): free terminal CH3 (H=3)
                # becomes CH2 (H=2) once bonded — allow H2 or H3.
                # sp3_c_anchor (e.g. ValAryl alpha-C bonded directly to aryl):
                # free CH (H=2 because R2 capped with H) becomes H=1 once bonded
                # to aryl — allow H1 or H2.
                if cts.get(bc_r) == "backbone_c_red":
                    _repl = f"[#{_sa.GetAtomicNum()};H2,H3{_cs}:{_MAP_OUT}]"
                elif cts.get(bc_r) == "sp3_c_anchor":
                    _repl = f"[#{_sa.GetAtomicNum()};H1,H2{_cs}:{_MAP_OUT}]"
                else:
                    _repl = f"[#{_sa.GetAtomicNum()}H{_nh}{_cs}:{_MAP_OUT}]"  # exact H count
            elif _si == bc_r_cap_idx:
                _repl = f"[*:{_MAP_CAP}]"  # temporary; will be removed below
            elif _si in cap_remove_cap_set:
                # Carboxyl R-group cap atom: keep as [#8] in standard SMARTS
                # (correct for D/E backbone detection).  The iso SMARTS removes
                # this atom entirely so E_g matches in isopeptide form.
                _repl = f"[#{_sa.GetAtomicNum()}{_cs}:{_MAP_CR_CAP}]"
            elif _sa.GetAtomicNum() in (7, 8, 16):
                # Side-chain heteroatoms (N/O/S) may form crosslinks/isopeptide
                # bonds that change their H count; match atom type only.
                _repl = f"[#{_sa.GetAtomicNum()}{_cs}]"
            elif _si in xl_cm_idxs:
                # Atom adjacent to a non-backbone R-group (e.g. terminal alkene
                # for olefin-staple crosslinks).  Crosslinking changes H count
                # (=CH2 → =CH-), so don't pin it.
                _repl = f"[#{_sa.GetAtomicNum()}{_cs}]"
            else:
                _repl = f"[#{_sa.GetAtomicNum()}H{_nh}{_cs}]"
            _smarts_str = _re2.sub(
                r"\[[^\]]*:" + str(_mn) + r"\]", _repl, _smarts_str, count=1
            )
        smarts_mol = _C.MolFromSmarts(_smarts_str)
        if smarts_mol is None:
            continue
        # Remove bc_r_cap atom ([*:9903]) so the query imposes no constraint on
        # the C-terminal context (amide bond, am cap, free acid all match).
        _cap_idxs_in_sm = sorted(
            [
                a.GetIdx()
                for a in smarts_mol.GetAtoms()
                if a.GetAtomMapNum() == _MAP_CAP
            ],
            reverse=True,
        )
        if _cap_idxs_in_sm:
            rw_sm0 = _RW(smarts_mol)
            for _cap_idx in _cap_idxs_in_sm:
                rw_sm0.RemoveAtom(_cap_idx)
            smarts_mol = rw_sm0.GetMol()
        # Build iso SMARTS by further removing the carboxyl cap atoms (tagged
        # _MAP_CR_CAP).  This lets E_g match in isopeptide form (no -OH at
        # gamma-carboxyl) without pulling K's epsilon-N into the match set.
        _cr_cap_idxs = sorted(
            [
                a.GetIdx()
                for a in smarts_mol.GetAtoms()
                if a.GetAtomMapNum() == _MAP_CR_CAP
            ],
            reverse=True,
        )
        if _cr_cap_idxs:
            rw_iso = _RW(smarts_mol)
            for _ci in _cr_cap_idxs:
                rw_iso.RemoveAtom(_ci)
            smarts_mol_iso = rw_iso.GetMol()
            smarts_in_n_idx_iso = next(
                (
                    a.GetIdx()
                    for a in smarts_mol_iso.GetAtoms()
                    if a.GetAtomMapNum() == _MAP_IN
                ),
                None,
            )
            smarts_out_co_idx_iso = next(
                (
                    a.GetIdx()
                    for a in smarts_mol_iso.GetAtoms()
                    if a.GetAtomMapNum() == _MAP_OUT
                ),
                None,
            )
            if smarts_in_n_idx_iso is None or smarts_out_co_idx_iso is None:
                smarts_mol_iso = None
                smarts_in_n_idx_iso = None
                smarts_out_co_idx_iso = None
            else:
                rw_iso2 = _RW(smarts_mol_iso)
                for _sa in rw_iso2.GetAtoms():
                    _sa.SetAtomMapNum(0)
                smarts_mol_iso = rw_iso2.GetMol()
        else:
            smarts_mol_iso = None
            smarts_in_n_idx_iso = None
            smarts_out_co_idx_iso = None
        # Strip _MAP_CR_CAP from standard smarts_mol (atoms stay, map num goes).
        for _sa in smarts_mol.GetAtoms():
            if _sa.GetAtomMapNum() == _MAP_CR_CAP:
                _sa.SetAtomMapNum(0)
        # Locate in_n / out_co in smarts_mol by map num, then strip map nums.
        smarts_in_n_idx = next(
            (a.GetIdx() for a in smarts_mol.GetAtoms() if a.GetAtomMapNum() == _MAP_IN),
            None,
        )
        smarts_out_co_idx = next(
            (
                a.GetIdx()
                for a in smarts_mol.GetAtoms()
                if a.GetAtomMapNum() == _MAP_OUT
            ),
            None,
        )
        if smarts_in_n_idx is None or smarts_out_co_idx is None:
            continue
        rw_sm = _RW(smarts_mol)
        for _sa in rw_sm.GetAtoms():
            _sa.SetAtomMapNum(0)
        smarts_mol = rw_sm.GetMol()
        rw2 = _RW(cm)
        for a in rw2.GetAtoms():
            if a.GetAtomMapNum() in (_MAP_IN, _MAP_OUT):
                a.SetAtomMapNum(0)
        orig_cm = rw2.GetMol()
        norm_cm = _s2c_normalize(orig_cm)
        if norm_cm is None:
            continue
        n_rings = _RDM.CalcNumRings(orig_cm)
        _BS = _RC.BondStereo
        has_ez = any(
            b.GetStereo() not in (_BS.STEREONONE, _BS.STEREOANY)
            for b in orig_cm.GetBonds()
            if b.GetBondTypeAsDouble() == 2.0
        )
        _bb_path = _C.GetShortestPath(raw_mol, in_n_raw, out_co_raw)
        _bb_dist = len(_bb_path) - 1 if _bb_path else 0
        lib.append(
            _CovEntry(
                abbr=abbr,
                m_type=m_type,
                orig_mol=orig_cm,
                norm_mol=norm_cm,
                smarts_mol=smarts_mol,
                in_n_idx=in_n_idx,
                out_co_idx=out_co_idx,
                smarts_in_n_idx=smarts_in_n_idx,
                smarts_out_co_idx=smarts_out_co_idx,
                n_atoms=norm_cm.GetNumAtoms(),
                n_rings=n_rings,
                has_ez=has_ez,
                bb_dist=_bb_dist,
                smarts_mol_iso=smarts_mol_iso,
                smarts_in_n_idx_iso=smarts_in_n_idx_iso,
                smarts_out_co_idx_iso=smarts_out_co_idx_iso,
            )
        )
        seen.add(abbr)
    lib.sort(key=lambda e: e.n_atoms, reverse=True)
    return lib


def _s2c_place_monomers(mol, cov_lib):
    """Return all placements of backbone monomers in mol as _PlacedNode objects.

    Match phase uses useChirality=True so stereo-specific library entries only
    match Cα atoms with the corresponding chirality (when the source HAS
    stereo).  Source atoms with unspecified chirality match BOTH L and D
    library SMARTS; the L-preference in the dedup sort below picks L by
    default (per project policy: no-stereo source → default L).
    """
    placements: list = []
    for entry in cov_lib:
        matches = mol.GetSubstructMatches(entry.smarts_mol, useChirality=True)
        for match in matches:
            placements.append(
                _PlacedNode(
                    abbr=entry.abbr,
                    m_type=entry.m_type,
                    mol_atoms=frozenset(match),
                    in_n=match[entry.smarts_in_n_idx],
                    out_co=match[entry.smarts_out_co_idx],
                    entry=entry,
                    match_tuple=match,
                )
            )
    return placements


def _s2c_build_backbone(mol, placements):
    """Find the longest non-overlapping backbone chain via amide-bond DAG.

    A→B when A.out_co bonds to B.in_n.  Among equal-length chains prefers
    those with the most m_type='aa' members (avoids picking lipid branch chains
    of the same length as the backbone, e.g. C20FA-AEEA-E_g vs K-G-K).
    Returns list of _PlacedNode in N→C order.
    """
    from collections import defaultdict as _dd3

    if not placements:
        return []

    # Deduplicate by atom coverage: L and D amino acids share identical heavy-atom
    # L and D amino acids share identical heavy-atom SMARTS (useChirality=False)
    # so they produce the same placement — _s2c_match (step 4) resolves stereo.
    # HOWEVER beta/gamma backbone variants (E vs E_g, D vs D_b) cover identical
    # atoms but use a DIFFERENT out_co (different carbonyl as R2).  Both must
    # survive deduplication so the chain DAG can select whichever one's R2 is the
    # actual amide-bond exit — that is the (1,2) backbone connection.
    # Dedup key: (mol_atoms, out_co) keeps one entry per distinct backbone exit.
    # Sort priority: aa type, carbonyl-backbone (prefer standard amide over reduced-
    # amide/no-carbonyl monomers), larger residue, more R-group slots, abbr.
    def _extra_rgroups(entry):
        # Count dummies with isotope >= 3 (R3, R4, ...) — proxies for "this monomer
        # supports crosslink/side-chain bonds beyond plain backbone amide".
        return sum(
            1
            for a in entry.orig_mol.GetAtoms()
            if a.GetAtomicNum() == 0 and a.GetIsotope() >= 3
        )

    def _has_carbonyl_backbone(entry):
        # The out_co atom of a standard amide-backbone monomer has a =O neighbour.
        # Reduced-amide (backbone_c_red) monomers don't — they shouldn't outrank
        # standard amino acids whose backbone atoms happen to overlap.
        a = entry.orig_mol.GetAtomWithIdx(entry.out_co_idx)
        for nb in a.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                b = entry.orig_mol.GetBondBetweenAtoms(entry.out_co_idx, nb.GetIdx())
                if b and b.GetBondTypeAsDouble() == 2.0:
                    return True
        return False

    # If any placement has a carbonyl backbone, drop reduced-amide (no-carbonyl)
    # placements globally — reduced-amide monomers (redG2/redG3) are only
    # meaningful for pure-polyamine molecules. Otherwise they'd hijack lipid
    # side-chain N-CH2-CH2 segments (e.g. K's epsilon-amino lipidation linker)
    # and produce wrong-output linear chains.
    _any_carbonyl = any(_has_carbonyl_backbone(p.entry) for p in placements)
    if _any_carbonyl:
        placements = [p for p in placements if _has_carbonyl_backbone(p.entry)]
    placements.sort(
        key=lambda p: (
            p.m_type != "aa",
            not _has_carbonyl_backbone(p.entry),  # carbonyl backbone preferred
            -p.entry.n_atoms,
            -_extra_rgroups(p.entry),
            p.abbr.startswith("D_"),  # prefer L_* over D_* (no-stereo source → L)
            p.abbr,
        )
    )
    seen_mol_out: set = set()
    deduped: list = []
    for p in placements:
        key = (p.mol_atoms, p.out_co)
        if key not in seen_mol_out:
            seen_mol_out.add(key)
            deduped.append(p)

    in_n_map = _dd3(list)
    for p in deduped:
        in_n_map[p.in_n].append(p)

    # ── Scaffold-portal pre-pass ──────────────────────────────────────────────
    # Detect multi-arm scaffold instances in the molecule (TBMB-style linkers,
    # PhosOxScaffold, ImzScaffold). Build atom_to_scaffold_anchors so the chain
    # walker can "tunnel through" a scaffold from one anchor to another when no
    # direct out_co→in_n bond exists. Anchor atom = the molecule atom that maps
    # to a scaffold's [*:n] attachment-point slot; other_anchor_atoms = list of
    # the other anchors of the SAME scaffold instance.
    atom_to_scaffold_anchors: dict = (
        {}
    )  # atom_idx → list[(other_anchor_atom_set, scaffold_abbr)]
    _sc_patterns = _s2c_get_scaffold_patterns()
    for _sc_entry in _sc_patterns:
        _sc_abbr = _sc_entry[0]
        _sc_smarts = _sc_entry[1]
        _sc_r_pos = _sc_entry[2]  # smarts_atom_idx → r_group_num
        for _sc_match in mol.GetSubstructMatches(_sc_smarts, useChirality=False):
            _anchor_atoms = {_sc_match[_si] for _si in _sc_r_pos}
            for _ai in _anchor_atoms:
                _others = frozenset(_anchor_atoms - {_ai})
                atom_to_scaffold_anchors.setdefault(_ai, []).append((_others, _sc_abbr))

    # out_cos from placements that don't share atoms with a given node — these
    # are genuine predecessor exits from OTHER residues.  Excludes the node's own
    # out_co and any alternate-backbone out_co for the same atom set (e.g. E and
    # E_g share mol_atoms; neither should count as the other's predecessor).
    def _external_out_cos(node):
        return {p.out_co for p in deduped if not (p.mol_atoms & node.mol_atoms)}

    def has_predecessor(node):
        ext = _external_out_cos(node)
        for nb in mol.GetAtomWithIdx(node.in_n).GetNeighbors():
            if nb.GetIdx() in ext:
                return True
        # Scaffold-portal predecessor: in_n is reachable through a scaffold's
        # anchors from some external out_co.  Check in_n itself AND its
        # neighbors as potential scaffold anchors (in_n can overlap with a
        # scaffold attachment point, e.g. CP01557 where ValAryl.in_n is also
        # an ImzScaffold anchor).
        _pred_check = [node.in_n] + [
            nb.GetIdx() for nb in mol.GetAtomWithIdx(node.in_n).GetNeighbors()
        ]
        for _pc in _pred_check:
            for other_anchors, _ in atom_to_scaffold_anchors.get(_pc, ()):
                for _far_anchor in other_anchors:
                    if _far_anchor in ext:
                        return True
                    for _far_nb in mol.GetAtomWithIdx(_far_anchor).GetNeighbors():
                        if _far_nb.GetIdx() in ext:
                            return True
        return False

    def successors(node, used_atoms):
        result = []
        seen_q = set()  # dedup if same q reachable via direct + scaffold paths
        # Rule 1: direct neighbour chain step (existing).
        for nb in mol.GetAtomWithIdx(node.out_co).GetNeighbors():
            for q in in_n_map[nb.GetIdx()]:
                if q.in_n not in used_atoms and q.out_co not in used_atoms:
                    if id(q) not in seen_q:
                        seen_q.add(id(q))
                        result.append(q)
        # Rule 2: scaffold-portal chain step — out_co itself or a neighbor
        # bonds to a scaffold anchor; check other anchors of that scaffold
        # for downstream in_n matches.  Include out_co itself because it can
        # overlap with a scaffold attachment point (e.g. CP01557 where
        # ValAryl.out_co is also an ImzScaffold anchor).
        _succ_check = [node.out_co] + [
            nb.GetIdx() for nb in mol.GetAtomWithIdx(node.out_co).GetNeighbors()
        ]
        for _sc in _succ_check:
            for other_anchors, _sc_abbr in atom_to_scaffold_anchors.get(_sc, ()):
                if any(a in used_atoms for a in other_anchors):
                    continue
                for _far_anchor in other_anchors:
                    for q in in_n_map.get(_far_anchor, []):
                        if q.in_n in used_atoms or q.out_co in used_atoms:
                            continue
                        if id(q) not in seen_q:
                            seen_q.add(id(q))
                            result.append(q)
                    for _far_nb in mol.GetAtomWithIdx(_far_anchor).GetNeighbors():
                        for q in in_n_map[_far_nb.GetIdx()]:
                            if q.in_n in used_atoms or q.out_co in used_atoms:
                                continue
                            if id(q) in seen_q:
                                continue
                            seen_q.add(id(q))
                            result.append(q)
        result.sort(key=lambda q: (q.m_type != "aa", -q.entry.n_atoms))
        return result

    starts = [p for p in deduped if not has_predecessor(p)]
    is_cyclic_topology = not starts
    if is_cyclic_topology:
        # Pure cyclic peptide: any rotation is valid. Only try one start.
        starts = deduped
        starts.sort(key=lambda p: (p.m_type != "aa", -p.entry.n_atoms))
        starts = starts[:1]
    else:
        # Mixed topology: linear lipid branches co-exist with a cyclic backbone.
        # Residues on the cyclic backbone all have predecessors → not in starts.
        # Explicitly add them so the DFS can find the (longer) cyclic chain.
        cyclic_cands = [p for p in deduped if has_predecessor(p)]
        starts = starts + cyclic_cands
        starts.sort(key=lambda p: (p.m_type != "aa", -p.entry.n_atoms))

    _all_placed = frozenset().union(*(p.mol_atoms for p in deduped))
    _true_starts = frozenset(p for p in deduped if not has_predecessor(p))

    best: list = []
    best_score = (-1, -1, -1, -1, 1, -1)

    def chain_score(c):
        _is_true_start = c[0] in _true_starts
        _in_n_unclaimed = any(
            _nb.GetAtomicNum() > 1 and _nb.GetIdx() not in _all_placed
            for _nb in mol.GetAtomWithIdx(c[0].in_n).GetNeighbors()
        )
        return (
            len(c),
            int(_is_true_start),
            int(_in_n_unclaimed),
            sum(1 for nd in c if nd.m_type == "aa"),
            -sum(nd.entry.bb_dist for nd in c),
            sum(len(nd.mol_atoms) for nd in c),
        )

    def dfs(node, chain, used):
        nonlocal best_score
        sc = chain_score(chain)
        if sc > best_score:
            best[:] = chain[:]
            best_score = sc
        for nxt in successors(node, used)[:3]:  # limit branching factor
            chain.append(nxt)
            dfs(nxt, chain, used | nxt.mol_atoms)
            chain.pop()

    for start in starts:
        dfs(start, [start], start.mol_atoms)

    return best


# ── Cap query-mol library (built once per SDF version) ───────────────────────
_s2c_n_cap_lib_cache: list | None = None
_s2c_c_cap_lib_cache: list | None = None
_s2c_cap_lib_stamp: tuple | None = None
_s2c_cap_lib_lock = threading.Lock()


def _s2c_get_cap_libs():
    """Return (n_cap_lib, c_cap_lib), building and caching on the first call per SDF version.

    n_cap_lib: [(n_atoms, abbr, query_mol), ...] sorted by -n_atoms.
               Query mol has the backbone_c dummy replaced with N, others removed.
    c_cap_lib: [(n_atoms, abbr, query_mol), ...] sorted by (-n_atoms, abbr).
               Query mol has the backbone_n dummy replaced with C, others removed.
    """
    global _s2c_n_cap_lib_cache, _s2c_c_cap_lib_cache, _s2c_cap_lib_stamp
    stamp = _s2c_library_stamp()
    with _s2c_cap_lib_lock:
        if _s2c_n_cap_lib_cache is not None and _s2c_cap_lib_stamp == stamp:
            return _s2c_n_cap_lib_cache, _s2c_c_cap_lib_cache
        from rdkit import Chem as _C
        from rdkit.Chem import Atom, RWMol

        _, by_abbr = _load_sdf()
        n_lib: list = []
        c_lib: list = []
        for abbr, cap_mol in by_abbr.items():
            p = cap_mol.GetPropsAsDict()
            if p.get("m_type", "") != "cap":
                continue
            chem_types = p.get("m_chem_types", "")
            has_bn = "backbone_n" in chem_types
            has_bc = "backbone_c" in chem_types
            if has_bc and not has_bn:
                attach, placeholder, target = "backbone_c", "N", n_lib
            elif has_bn and not has_bc:
                attach, placeholder, target = "backbone_n", "C", c_lib
            else:
                continue
            r_num = None
            for part in chem_types.split(","):
                part = part.strip()
                if ":" not in part:
                    continue
                rn, rtype = part.split(":", 1)
                if rtype.strip() == attach:
                    try:
                        r_num = int(rn.strip())
                        break
                    except ValueError:
                        pass
            if r_num is None:
                continue
            dummy_idx = next(
                (
                    a.GetIdx()
                    for a in cap_mol.GetAtoms()
                    if a.GetAtomicNum() == 0 and a.GetIsotope() == r_num
                ),
                None,
            )
            if dummy_idx is None:
                continue
            rw2 = RWMol(cap_mol)
            rw2.ReplaceAtom(dummy_idx, Atom(placeholder))
            _rg_list = [r.strip() for r in p.get("m_Rgroups", "").split(",")]
            _resolve_dummy_rgroups(rw2, _rg_list, r_num)
            try:
                _C.SanitizeMol(rw2)
            except (ValueError, RuntimeError):
                continue
            target.append((cap_mol.GetNumAtoms(), abbr, rw2.GetMol()))
        n_lib.sort(key=lambda x: -x[0])
        c_lib.sort(key=lambda x: (-x[0], x[1].lstrip("_")))
        _s2c_n_cap_lib_cache = n_lib
        _s2c_c_cap_lib_cache = c_lib
        _s2c_cap_lib_stamp = stamp
        return n_lib, c_lib


def _s2c_identify_cap_from_atom(
    mol, anchor_idx, cap_start_atom, excluded_atoms, placeholder_sym, cap_lib
):
    """Shared BFS + fragment-build + match kernel for both N-cap and C-cap identification.

    anchor_idx:     backbone atom the cap is bonded to (N for N-caps, carbonyl C for C-caps).
    cap_start_atom: first atom on the cap side of the anchor bond.
    placeholder_sym: element symbol of the anchor placeholder atom ('N' or 'C').
    cap_lib:        pre-built [(n_atoms, abbr, query_mol), ...] from _s2c_get_cap_libs().
    """
    from rdkit import Chem as _C
    from rdkit.Chem import Atom, RWMol

    cap_atoms: set = set()
    queue = [cap_start_atom]
    while queue:
        ai = queue.pop()
        if ai in cap_atoms or ai == anchor_idx or ai in excluded_atoms:
            continue
        cap_atoms.add(ai)
        for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
            ni = nb.GetIdx()
            if ni not in cap_atoms and ni != anchor_idx and ni not in excluded_atoms:
                queue.append(ni)
    if not cap_atoms:
        return None
    rw = RWMol()
    amap: dict = {}
    for ai in cap_atoms:
        amap[ai] = rw.AddAtom(mol.GetAtomWithIdx(ai))
    anchor_placeholder = rw.AddAtom(Atom(placeholder_sym))
    for bond in mol.GetBonds():
        bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if bi in cap_atoms and ei in cap_atoms:
            rw.AddBond(amap[bi], amap[ei], bond.GetBondType())
        elif bi in cap_atoms and ei == anchor_idx:
            rw.AddBond(amap[bi], anchor_placeholder, bond.GetBondType())
        elif ei in cap_atoms and bi == anchor_idx:
            rw.AddBond(amap[ei], anchor_placeholder, bond.GetBondType())
    try:
        _C.SanitizeMol(rw)
    except (ValueError, RuntimeError):
        return None
    cap_frag = rw.GetMol()
    # Require EXACT coverage: the cap query must cover every heavy atom of
    # cap_frag.  HasSubstructMatch alone is a substructure test, so e.g.
    # NHEt (4 atoms) would falsely match the first 4 atoms of a 24-atom
    # macrocyclic linker that wraps back to the N-terminus — the dominant
    # rt-mismatch bug in cyclicpepedia (CP00011 etc., 155/383 mismatches
    # as of 2026-05-16).  Insisting on full coverage rejects those false
    # positives while leaving true linear caps unaffected.
    cap_frag_size = cap_frag.GetNumAtoms()
    for _, abbr, qmol in cap_lib:
        m = cap_frag.GetSubstructMatch(qmol, useChirality=False)
        if m and len(m) == cap_frag_size:
            return abbr
    return None


def _s2c_identify_n_cap_from_atom(mol, n_idx, cap_start_atom, excluded_atoms):
    """Identify N-cap; returns abbreviation string or None."""
    n_lib, _ = _s2c_get_cap_libs()
    return _s2c_identify_cap_from_atom(
        mol, n_idx, cap_start_atom, excluded_atoms, "N", n_lib
    )


def _s2c_identify_c_cap_from_atom(mol, co_idx, cap_start_atom, excluded_atoms):
    """Identify C-cap; returns abbreviation string or None."""
    _, c_lib = _s2c_get_cap_libs()
    return _s2c_identify_cap_from_atom(
        mol, co_idx, cap_start_atom, excluded_atoms, "C", c_lib
    )


def _s2c_cap_standalone(smi, rgroups_str):
    """Replace [n*] dummies with their leaving-group atom for standalone matching.

    Delegates to _s2c_cap_smiles, which now handles the full leaving-group table
    ([OH]→O, [Br]→Br, [Cl]→Cl, [I]→I, [SH]→S, else→[H]).
    """
    return _s2c_cap_smiles(smi, rgroups_str)


def _s2c_find_amide_bonds(mol):
    """Return list of (n_idx, co_idx) for every N-C(=O) amide bond in mol."""
    amide_bonds = []
    for bond in mol.GetBonds():
        bi, ei = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        for n_idx, c_idx in ((bi, ei), (ei, bi)):
            a_n = mol.GetAtomWithIdx(n_idx)
            a_c = mol.GetAtomWithIdx(c_idx)
            if a_n.GetAtomicNum() != 7 or a_c.GetAtomicNum() != 6:
                continue
            if any(
                mol.GetAtomWithIdx(nb.GetIdx()).GetAtomicNum() == 8
                and mol.GetBondBetweenAtoms(c_idx, nb.GetIdx()).GetBondTypeAsDouble()
                == 2.0
                for nb in mol.GetAtomWithIdx(c_idx).GetNeighbors()
            ):
                amide_bonds.append((n_idx, c_idx))
                break
    return amide_bonds


def _s2c_multi_amide_chain(mol, amide_pairs):
    """Walk a chain of ≥2 amide bonds when normal backbone detection failed.

    Finds the start amide (CO-side not reachable from any other amide N), then
    walks segment by segment toward the C-terminus.  Each middle segment's atoms
    — from the amide N (inclusive) to the next amide CO (inclusive, but excluding
    the next amide's N) — are isolated and matched against the backbone library.
    The N-cap and C-cap are identified using the existing cap-helper functions.

    Robustness: when multiple amide COs are reachable from the current N (e.g.
    sidechain amides in Asn/Gln), the one found first by BFS (shortest bond path)
    is chosen as the backbone CO.  This works for most standard residues; unusual
    branching may produce incorrect results, so the function returns None on any
    matching failure rather than raising.

    Returns (cabiln_str, details) or None.
    """
    from collections import deque as _dq2

    all_amide_cos = {co for _, co in amide_pairs}
    by_co = {co: (n, co) for n, co in amide_pairs}

    # Find the start amide: its CO-side is NOT reachable from any other amide N.
    def _co_reached_from_other_n(n_idx, co_idx):
        for other_n, other_co in amide_pairs:
            if other_n == n_idx:
                continue
            visited = {other_co}
            q = _dq2([other_n])
            while q:
                ai = q.popleft()
                if ai in visited:
                    continue
                visited.add(ai)
                if ai == co_idx:
                    return True
                for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                    if nb.GetIdx() not in visited:
                        q.append(nb.GetIdx())
        return False

    start = next(
        ((n, co) for n, co in amide_pairs if not _co_reached_from_other_n(n, co)),
        None,
    )
    if start is None:
        return None  # cyclic or complex topology

    start_n, start_co = start
    n_cap = _s2c_identify_n_cap_from_atom(mol, start_n, start_co, set())
    if n_cap is None:
        return None

    lib = _s2c_get_lib()
    middle_abbrs = []
    current_n = start_n
    prev_co = start_co

    while True:
        # BFS from current_n (excluding prev_co) to find the nearest amide CO.
        visited = {prev_co}
        q = _dq2([current_n])
        next_co = None
        while q:
            ai = q.popleft()
            if ai in visited:
                continue
            visited.add(ai)
            if ai != current_n and ai in all_amide_cos:
                next_co = ai
                break
            for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                if nb.GetIdx() not in visited:
                    q.append(nb.GetIdx())

        if next_co is None:
            # No more amide COs reachable → current_n starts the C-cap.
            c_cap = _s2c_identify_c_cap_from_atom(mol, prev_co, current_n, set())
            break

        next_amide = by_co.get(next_co)
        if next_amide is None:
            return None
        next_n = next_amide[0]

        # Collect middle-segment atoms: BFS from current_n, excluding prev_co
        # and next_n.  This includes next_co and its =O so isolation can cap it
        # with -OH (giving the correct free amino acid form for _s2c_match).
        seg_atoms: set = set()
        vis2 = {prev_co, next_n}
        sq = _dq2([current_n])
        while sq:
            ai = sq.popleft()
            if ai in vis2:
                continue
            vis2.add(ai)
            seg_atoms.add(ai)
            for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                if nb.GetIdx() not in vis2:
                    sq.append(nb.GetIdx())

        if not seg_atoms:
            return None

        aas, _ = _s2c_isolate_residues(mol, [list(seg_atoms)])
        aa_mol = aas[0] if aas else None
        if aa_mol is None:
            return None

        abbr, _ = _s2c_match(aa_mol, lib)
        if abbr is None:
            return None
        middle_abbrs.append(abbr)

        current_n = next_n
        prev_co = next_co

    if c_cap is None:
        return None

    all_abbrs = [n_cap] + middle_abbrs + [c_cap]
    return "-".join(all_abbrs), [(a, 1.0, 0) for a in all_abbrs]


def _s2c_canonicalize_guanidinium(mol):
    """Flip Cδ-N=C(N)N guanidinium tautomers to the Cδ-N-C(=N)N form.

    Arg's library SMARTS bakes in the Cδ-Nε single-bond / Cζ=Nη double-bond
    tautomer (the form stored in monomers.sdf).  PubChem-derived SMILES
    frequently encode the alternate Cδ-N=Cζ(N)N tautomer, which has the same
    heavy-atom connectivity but different bond orders.  RDKit's substructure
    matcher is bond-order sensitive, so Arg's SMARTS misses the alternate
    tautomer; the Cδ-N then matches Orn's relaxed sidechain-N, and the walker
    silently mis-types Arg as Orn (dropping the 3 guanidinium atoms).  Run
    this once on input, before placement, to keep the library SMARTS valid.
    """
    from rdkit import Chem as _Chem
    from rdkit.Chem import RWMol as _RWMol

    pat = _Chem.MolFromSmarts("[CX4][NX2;!R]=[CX3]([NX3])[NX3]")
    matches = mol.GetSubstructMatches(pat)
    if not matches:
        return mol
    rw = _RWMol(mol)
    for _cd, ne, cz, nh1, _nh2 in matches:
        b_ne_cz = rw.GetBondBetweenAtoms(ne, cz)
        b_cz_nh1 = rw.GetBondBetweenAtoms(cz, nh1)
        if b_ne_cz is None or b_cz_nh1 is None:
            continue
        b_ne_cz.SetBondType(_Chem.BondType.SINGLE)
        b_cz_nh1.SetBondType(_Chem.BondType.DOUBLE)
        for _ai in (ne, nh1):
            _a = rw.GetAtomWithIdx(_ai)
            _a.SetNumExplicitHs(0)
            _a.SetNoImplicit(False)
    try:
        _Chem.SanitizeMol(rw)
    except (ValueError, RuntimeError):
        return mol
    return rw.GetMol()


def _s2c_normalize_iminol(mol):
    """Return a copy of mol with any iminol tautomers flipped to the amide form.

    Iminol: N=C(OH)  →  amide: NH-C(=O)
    This is a bond-order surgery (same connectivity, different bond orders) that
    sometimes appears in database SMILES for peptide backbones.  Applied only as
    a rescue pass — never mutates the original mol.

    Returns the modified molecule if any iminol bonds were found, otherwise
    returns the original mol object unchanged (caller can test ``is mol``).
    """
    from rdkit import Chem as _Chem
    from rdkit.Chem import RWMol as _RWMol

    # strict: sp2 N=sp2 C bearing a hydroxyl — the genuine amide-iminol tautomer.
    # Excludes aromatic n (won't carry an explicit '=') and amidoximes/oximes
    # (OH on N, not C), so it is a safe no-op on ordinary amide-form peptides.
    _IMINOL_SMARTS = _Chem.MolFromSmarts("[NX2:1]=[CX3:2][OX2H1:3]")
    matches = mol.GetSubstructMatches(_IMINOL_SMARTS, useChirality=False)
    if not matches:
        return mol
    rw = _RWMol(mol)
    for n_idx, c_idx, o_idx in matches:
        n_atom = rw.GetAtomWithIdx(n_idx)
        o_atom = rw.GetAtomWithIdx(o_idx)
        nc_bond = rw.GetBondBetweenAtoms(n_idx, c_idx)
        co_bond = rw.GetBondBetweenAtoms(c_idx, o_idx)
        if nc_bond is None or co_bond is None:
            continue
        nc_bond.SetBondType(_Chem.BondType.SINGLE)
        co_bond.SetBondType(_Chem.BondType.DOUBLE)
        n_atom.SetNumExplicitHs(0)
        n_atom.SetNoImplicit(False)
        o_atom.SetNumExplicitHs(0)
        o_atom.SetNoImplicit(False)
    try:
        _Chem.SanitizeMol(rw)
    except (ValueError, RuntimeError):
        return mol
    return rw.GetMol()


def _s2c_backbone_less_fallback(mol):
    """Identify backbone-less molecules: single entities, two-cap, or multi-amide chains.

    Called when smiles_to_cabiln_core finds no backbone monomers (no residue
    with both R1 and R2).  Handles in order:

    1. Whole-molecule single entity — SMILES matches one library monomer exactly
       after capping dummies with the actual leaving-group atom.  Catches
       standalone scaffolds (TBMB, TATA, …) and single cap molecules.

    2. Single-amide two-cap chain — one N-C(=O) joins an N-cap to a C-cap
       (e.g. ac-am, fmoc-am, boc-am).

    3. Multi-amide chain (≥2 amide bonds) — N-cap + zero or more middle residues
       + C-cap, identified by walking the amide-bond graph.  Middle segments are
       isolated and matched against the backbone library.  Handles molecules whose
       coverage-library SMARTS failed to match despite the residues being known.

    Returns (cabiln_str, match_details) or (None, None).
    """
    from rdkit import Chem as _C

    input_can = _C.MolToSmiles(mol)

    # ── 1. Whole-molecule match ───────────────────────────────────────────────
    rlib = _s2c_get_raw_lib()
    for abbr, raw_mol in rlib.items():
        p = raw_mol.GetPropsAsDict()
        rg_str = p.get("m_Rgroups", "")
        raw_smi = _C.MolToSmiles(raw_mol)
        try:
            standalone_smi = _s2c_cap_standalone(raw_smi, rg_str)
            sm = _C.MolFromSmiles(standalone_smi)
            if sm is not None and _C.MolToSmiles(sm) == input_can:
                return abbr, [(abbr, 1.0, mol.GetNumAtoms())]
        except (ValueError, RuntimeError):
            pass

    # ── 2 & 3. Amide-bond chain walk ─────────────────────────────────────────
    amide_bonds = _s2c_find_amide_bonds(mol)

    if len(amide_bonds) == 1:
        n_idx, co_idx = amide_bonds[0]
        n_cap = _s2c_identify_n_cap_from_atom(mol, n_idx, co_idx, set())
        c_cap = _s2c_identify_c_cap_from_atom(mol, co_idx, n_idx, set())
        if n_cap and c_cap:
            return f"{n_cap}-{c_cap}", [(n_cap, 1.0, 0), (c_cap, 1.0, 0)]

    if len(amide_bonds) >= 2:
        result = _s2c_multi_amide_chain(mol, amide_bonds)
        if result is not None:
            return result

    return None, None


def _s2c_local_residue_template(molecule, node):
    """Preserve a placed residue with its actual backbone attachment atoms.

    This deliberately narrow rescue accepts only a residue whose boundaries
    are single bonds at its identified backbone N and carbonyl C. Sidechain
    attachments and partially covered sidechains are declined. FragmentOnBonds
    retains source atom ordering and stereo, including explicit imine H atoms.
    """
    from rdkit import Chem

    def carbonyl(atom):
        return atom.GetAtomicNum() == 6 and any(
            other.GetAtomicNum() == 8
            and molecule.GetBondBetweenAtoms(
                atom.GetIdx(), other.GetIdx()
            ).GetBondType()
            == Chem.BondType.DOUBLE
            for other in atom.GetNeighbors()
        )

    if molecule.GetAtomWithIdx(node.in_n).GetAtomicNum() != 7 or not carbonyl(
        molecule.GetAtomWithIdx(node.out_co)
    ):
        return None
    atoms = set(node.mol_atoms)
    for atom_index in tuple(atoms):
        atoms.update(
            atom.GetIdx()
            for atom in molecule.GetAtomWithIdx(atom_index).GetNeighbors()
            if atom.GetAtomicNum() == 1
        )
    # A free terminal hydroxyl is the R2 leaving group. Most coverage patterns
    # omit it already, but an entry may include it.
    for atom in molecule.GetAtomWithIdx(node.out_co).GetNeighbors():
        if (
            atom.GetAtomicNum() == 8
            and atom.GetDegree() == 1
            and atom.GetTotalNumHs() == 1
        ):
            atoms.discard(atom.GetIdx())

    cuts, slots = [], []
    for bond in molecule.GetBonds():
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if (begin in atoms) == (end in atoms):
            continue
        inside, outside = (begin, end) if begin in atoms else (end, begin)
        outside_atom = molecule.GetAtomWithIdx(outside)
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return None
        if inside == node.in_n and carbonyl(outside_atom):
            slot = 1
        elif inside == node.out_co and outside_atom.GetAtomicNum() in (7, 8):
            slot = 2
        else:
            return None
        if slot in slots:
            return None
        cuts.append(bond.GetIdx())
        slots.append(slot)
    if 2 not in slots:
        return None

    opened = Chem.FragmentOnBonds(
        molecule, cuts, dummyLabels=[(slot, slot) for slot in slots]
    )
    mappings = []
    fragments = Chem.GetMolFrags(opened, asMols=True, fragsMolAtomMapping=mappings)
    original_atoms = set(range(molecule.GetNumAtoms()))
    for fragment, mapping in zip(fragments, mappings):
        if set(mapping) & original_atoms != atoms:
            continue
        if 1 not in slots:
            rw = Chem.RWMol(fragment)
            nitrogen_index = mapping.index(node.in_n)
            nitrogen = rw.GetAtomWithIdx(nitrogen_index)
            if not nitrogen.GetTotalNumHs():
                return None
            if nitrogen.GetNumExplicitHs():
                nitrogen.SetNumExplicitHs(nitrogen.GetNumExplicitHs() - 1)
            dummy = Chem.Atom(0)
            dummy.SetIsotope(1)
            dummy_index = rw.AddAtom(dummy)
            rw.AddBond(nitrogen_index, dummy_index, Chem.BondType.SINGLE)
            fragment = rw.GetMol()
            Chem.SanitizeMol(fragment)
        return f"<{Chem.MolToSmiles(fragment)}>"
    return None


def _s2c_decompose(mol, *, source_molecule=None):
    """Convert a peptide SMILES string to CABILN notation.

    Returns (cabiln_str, match_details) where match_details is a list of
    (abbr, n_matched, n_total_atoms) per residue in chain order.

    Uses library-first backbone detection: matches all backbone monomers against
    the full input mol, builds an amide-bond DAG, and finds the longest chain.
    This handles β-AAs, peptoids, and any future SDF-registered monomer types.
    """
    from collections import deque as _dq

    from rdkit import Chem as _C

    # ── 0. Pre-detect ring crosslinks; build de-reacted mol for backbone detection
    # Ring-forming reactions (e.g. CuAAC triazole) consume reactive sidechains
    # (alkyne, azide) into the ring, making coverage-library detection fail.
    # Bond surgery restores the pre-reaction groups without changing atom indices,
    # so the de-reacted mol's backbone atom indices are valid for the original mol.
    _raw_xlinks: list = []  # (matom_a, matom_b, slot_a, slot_b) – original indices
    _rxlink_pats = _s2c_get_ring_crosslink_patterns()
    _derx_mol = mol
    if _rxlink_pats:
        _raw_xlinks, _derx_mol = _s2c_pre_detect_ring_xlinks(mol, _rxlink_pats)

    # ── 1. Library-first backbone detection ──────────────────────────────────
    cov_lib = _s2c_get_coverage_lib()
    placements = _s2c_place_monomers(_derx_mol, cov_lib)
    # iso_placements: like placements but uses the iso SMARTS (carboxyl cap
    # removed) for monomers that have one.  Used only for branch detection so
    # that isopeptide-bonded monomers (E_g→K) are found without pulling the
    # backbone amide-N into the match set.
    iso_placements: list = []
    for _ie in cov_lib:
        if _ie.smarts_mol_iso is not None:
            _ism, _ini, _ioi = (
                _ie.smarts_mol_iso,
                _ie.smarts_in_n_idx_iso,
                _ie.smarts_out_co_idx_iso,
            )
        else:
            _ism, _ini, _ioi = (
                _ie.smarts_mol,
                _ie.smarts_in_n_idx,
                _ie.smarts_out_co_idx,
            )
        for _imatch in _derx_mol.GetSubstructMatches(_ism, useChirality=False):
            iso_placements.append(
                _PlacedNode(
                    abbr=_ie.abbr,
                    m_type=_ie.m_type,
                    mol_atoms=frozenset(_imatch),
                    in_n=_imatch[_ini],
                    out_co=_imatch[_ioi],
                    entry=_ie,
                    match_tuple=_imatch,
                )
            )
    backbone = _s2c_build_backbone(_derx_mol, placements)
    if not backbone:
        _fb_cabiln, _fb_details = _s2c_backbone_less_fallback(mol)
        if _fb_cabiln is not None:
            return _fb_cabiln, _fb_details
        raise ValueError("Cannot detect peptide backbone (no library monomers matched)")

    n = len(backbone)
    backbone_atoms_all: set = set()
    for node in backbone:
        backbone_atoms_all.update(node.mol_atoms)

    # ── 2. Cyclic detection ───────────────────────────────────────────────────
    cyclic = mol.GetBondBetweenAtoms(backbone[-1].out_co, backbone[0].in_n) is not None

    # ── 3. Cap detection ──────────────────────────────────────────────────────
    n_cap: str | None = None
    c_cap: str | None = None
    if not cyclic:
        # N-cap: carbonyl-bearing C bonded to backbone[0].in_n, outside backbone
        for nb in mol.GetAtomWithIdx(backbone[0].in_n).GetNeighbors():
            nbi = nb.GetIdx()
            if nbi not in backbone_atoms_all and nb.GetSymbol() == "C":
                has_co = any(
                    x.GetSymbol() == "O"
                    and mol.GetBondBetweenAtoms(nbi, x.GetIdx()).GetBondTypeAsDouble()
                    == 2.0
                    for x in nb.GetNeighbors()
                )
                if has_co:
                    # Identify the N-cap explicitly; do NOT fall back to 'ac'.
                    # The old `or 'ac'` fallback was masking macrocyclic
                    # ring-closures (where the in_n's external carbonyl is the
                    # NEXT residue around the ring, not a real cap), producing
                    # 246+ spurious 'ac' rt-mismatches in cyclicpepedia.  If the
                    # cap can't be identified, leave n_cap unset.
                    identified = _s2c_identify_n_cap_from_atom(
                        mol, backbone[0].in_n, nbi, backbone_atoms_all
                    )
                    if identified:
                        n_cap = identified
                        break
                    # else: keep looking at other neighbours of in_n.
                    continue
        # C-cap: any non-backbone substituent on out_co (am, NHEt, OEt, OtBu, …)
        for nb in mol.GetAtomWithIdx(backbone[-1].out_co).GetNeighbors():
            nbi = nb.GetIdx()
            if nbi in backbone_atoms_all:
                continue
            # Guard: reject peptide-bond Ns (alpha-N of the next residue whose
            # Cα is itself bonded to a carbonyl C).  A real cap N has no such
            # carbonyl-bearing C neighbour on its non-out_co side.
            if nb.GetSymbol() == "N":
                is_peptide_n = any(
                    any(
                        x.GetSymbol() == "O"
                        and mol.GetBondBetweenAtoms(
                            c.GetIdx(), x.GetIdx()
                        ).GetBondTypeAsDouble()
                        == 2.0
                        for x in c.GetNeighbors()
                    )
                    for c in nb.GetNeighbors()
                    if c.GetIdx() != backbone[-1].out_co and c.GetSymbol() == "C"
                )
                if is_peptide_n:
                    continue
            identified = _s2c_identify_c_cap_from_atom(
                mol, backbone[-1].out_co, nbi, backbone_atoms_all
            )
            if identified:
                c_cap = identified
                break
            # Fallback: plain am (NH2 with no additional C)
            if nb.GetSymbol() == "N":
                c_nbrs = [x for x in nb.GetNeighbors() if x.GetSymbol() == "C"]
                if len(c_nbrs) <= 1:
                    c_cap = "am"
                    break

    # ── 4. Match each backbone residue with stereo-aware matching ────────────
    from rdkit.Chem import AllChem as _AC3

    _AC3.AssignStereochemistry(mol, cleanIt=True, force=True)
    lib = _s2c_get_lib()
    details = []
    abbrs = []
    # Use de-reacted mol for residue isolation/matching when ring crosslinks
    # consumed reactive sidechains into a ring (e.g. CuAAC triazole consumes
    # Pra's alkyne and AzK's azide). _derx_mol has identical atom indices but
    # correct pre-reaction bond types so library matching works correctly.
    _match_mol = _derx_mol
    for pos, node in enumerate(backbone):
        atom_list = list(node.mol_atoms)
        # For free-acid C-terminus (no am cap): include the carboxyl OH in the
        # atom list so isolation gives COOH not CHO (which would match Gly_al).
        if not cyclic and pos == n - 1 and c_cap is None:
            for _nb in _match_mol.GetAtomWithIdx(node.out_co).GetNeighbors():
                _bd = _match_mol.GetBondBetweenAtoms(node.out_co, _nb.GetIdx())
                if (
                    _nb.GetIdx() not in backbone_atoms_all
                    and _nb.GetSymbol() == "O"
                    and _nb.GetTotalNumHs() > 0
                    and _bd.GetBondTypeAsDouble() == 1.0
                ):
                    atom_list.append(_nb.GetIdx())
                    break
        aas, _ = _s2c_isolate_residues(_match_mol, [atom_list])
        aa_mol = aas[0]
        # Terminal residues: strip cap artifact introduced by isolate_residues
        # (the cap atoms are outside node.mol_atoms so the isolated mol already
        # has the correct N—H / C=O—OH termini; only strip if isolation added
        # a spurious cap atom from the coverage-lib match itself).
        if not cyclic:
            if pos == 0 and n_cap is not None:
                stripped, did = _s2c_strip_n_cap(aa_mol)
                if did:
                    aa_mol = stripped
            if pos == n - 1 and c_cap is not None:
                stripped, did = _s2c_strip_c_cap(aa_mol)
                if did:
                    aa_mol = stripped
        # Filter lib to entries with the same backbone path length as the
        # backbone-chain placement.  This ensures E_g is preferred over E when
        # the gamma-COOH was the actual chain exit (bb_dist=4 vs 2).  Falls
        # back to the full lib if nothing matches (e.g. novel unregistered dist).
        _node_bb_dist = node.entry.bb_dist
        _lib_filtered = [e for e in lib if e[6] == _node_bb_dist]
        abbr, score = _s2c_match(aa_mol, _lib_filtered or lib)
        if abbr is None:
            raise ValueError(
                f"Unrecognised monomer at backbone position {pos} "
                f"({aa_mol.GetNumAtoms()} heavy atoms); "
                f"register it in the monomer library first"
            )
        details.append((abbr, score, aa_mol.GetNumAtoms()))
        abbrs.append(abbr)

    # ── 5. Branch detection and formatting ────────────────────────────────────
    # Non-backbone placements whose out_co bonds to a backbone atom are branch
    # anchors (e.g. E_g attached to K's epsilon-N via isopeptide bond).
    # Deduplicate branch placements by atom set: L/D pairs and similar monomers
    # (DGlu/E_g/E) cover identical heavy atoms.  When multiple abbrs match
    # the same atoms, prefer standard E/D (whose R4=γ/β-COOH gives E(4,4) etc.)
    # over gamma-form E_g at the isopeptide junction; E_g is preferred only in
    # orphan-segment matching (step 4 of _s2c_walk_branch) via the 1-2 rule.
    # Priority 0 = most preferred; unlisted abbrs get 999.
    _ISOPEP_PRIORITY: dict = {
        "E": 0,
        "dE": 1,  # standard Glu: slot annotation e.g. E(4,4)
        "D": 2,
        "dD": 3,  # standard Asp
        "DGlu": 4,
        "DGlu_g": 5,  # D-Glu variants
        "E_g": 6,  # gamma-backbone Glu (lower — E preferred at isopeptide junction)
    }
    _bp_by_key: dict = {}  # frozenset(mol_atoms) → best _PlacedNode
    for p in iso_placements:
        if p.mol_atoms & backbone_atoms_all:
            continue  # overlaps backbone — skip
        key = frozenset(p.mol_atoms)
        if key not in _bp_by_key:
            _bp_by_key[key] = p
        else:
            # Keep the placement with higher isopeptide priority
            cur_pri = _ISOPEP_PRIORITY.get(_bp_by_key[key].abbr, 999)
            new_pri = _ISOPEP_PRIORITY.get(p.abbr, 999)
            if new_pri < cur_pri:
                _bp_by_key[key] = p
    branch_placements = list(_bp_by_key.values())
    branch_junctions: dict = {}  # backbone_pos → list of {'main_at', 'glu_node'}

    for bi, bb_node in enumerate(backbone):
        _seen_main_ats: set = set()  # one branch per unique junction atom
        for bp in branch_placements:
            main_at = -1
            # Check ALL atoms of bp for a bond to any backbone atom — not just
            # bp.out_co, because E_g's isopeptide attachment is its gamma-C=O
            # (in mol_atoms) rather than its alpha-C=O (out_co).
            for bpa in bp.mol_atoms:
                for ba in bb_node.mol_atoms:
                    if mol.GetBondBetweenAtoms(bpa, ba) is not None:
                        main_at = ba
                        break
                if main_at != -1:
                    break
            if main_at != -1 and main_at not in _seen_main_ats:
                _seen_main_ats.add(main_at)
                branch_junctions.setdefault(bi, []).append(
                    {"main_at": main_at, "glu_node": bp}
                )

    if branch_junctions:
        raw_lib = _s2c_get_raw_lib()
        full_lib = _s2c_get_full_lib()
        branch_brackets: dict = {}

        for bi, junction_list in branch_junctions.items():
            anchor_abbr = abbrs[bi]
            bracket_acc = ""
            for jinfo in junction_list:
                main_at = jinfo["main_at"]
                glu_node = jinfo["glu_node"]
                glu_atoms = set(glu_node.mol_atoms)
                # in_n (alpha-N) is the junction from E_g toward AEEA/C20FA —
                # AEEA's out_co forms an amide bond with E_g's in_n.
                # out_co (alpha-C=O) is a free carboxyl terminus in the assembled
                # molecule when E_g is used as a branch (isopeptide via R4).
                glu_outgoing_n = glu_node.in_n

                # Orphan atoms: reachable from glu_outgoing_n, outside backbone+glu
                orphan_atoms: set = set()
                vis_o = backbone_atoms_all | glu_atoms
                q_o = _dq()
                for nb in mol.GetAtomWithIdx(glu_outgoing_n).GetNeighbors():
                    ni = nb.GetIdx()
                    if ni not in vis_o:
                        q_o.append(ni)
                while q_o:
                    ai = q_o.popleft()
                    if ai in vis_o:
                        continue
                    vis_o.add(ai)
                    orphan_atoms.add(ai)
                    for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                        if nb.GetIdx() not in vis_o:
                            q_o.append(nb.GetIdx())

                chain = _s2c_walk_branch(
                    mol,
                    glu_atoms,
                    glu_outgoing_n,
                    orphan_atoms,
                    full_lib,
                    raw_lib,
                    anchor_abbr,
                    main_at,
                    junction_abbr=glu_node.abbr,
                )
                if chain:
                    parts = ".".join(f"{a}({p},{c})" for a, p, c in chain)
                    bracket_acc += f".[{parts}]"
            if bracket_acc:
                branch_brackets[bi] = bracket_acc

        abbrs = [a + branch_brackets.get(i, "") for i, a in enumerate(abbrs)]

    # ── 5.3 Generic sidechain cap detection ──────────────────────────────────
    # After multi-monomer branch detection, atoms not in any backbone residue
    # or branch anchor (e.g. trt on Cys thiol, boc on Lys ε-N) are matched
    # against the full library as single-monomer sidechain caps.
    #
    # _sc_claimed tracks atoms already attributed so we don't double-annotate.
    _sc_claimed: set = set(backbone_atoms_all)
    for _jl in branch_junctions.values():
        for _j in _jl:
            _sc_claimed.update(_j["glu_node"].mol_atoms)

    _sc_unclaimed: set = set(range(mol.GetNumAtoms())) - _sc_claimed
    if _sc_unclaimed:
        from rdkit.Chem import RWMol as _RW2

        _sc_raw = _s2c_get_raw_lib()
        _sc_full = _s2c_get_full_lib()

        for _bi, _node in enumerate(backbone):
            _base_abbr = details[_bi][0]
            # Find all junction points: backbone atoms bonded to unclaimed atoms
            for _ai in list(_node.mol_atoms):
                for _nb in mol.GetAtomWithIdx(_ai).GetNeighbors():
                    _ni = _nb.GetIdx()
                    if _ni not in _sc_unclaimed:
                        continue
                    # BFS to collect the full connected unclaimed cluster
                    _cluster: set = set()
                    _q2 = _dq([_ni])
                    while _q2:
                        _ci = _q2.popleft()
                        if _ci in _sc_claimed or _ci in _cluster:
                            continue
                        _cluster.add(_ci)
                        for _cnb in mol.GetAtomWithIdx(_ci).GetNeighbors():
                            if _cnb.GetIdx() not in _sc_claimed:
                                _q2.append(_cnb.GetIdx())
                    if not _cluster:
                        continue
                    # Build a fragment mol from the cluster
                    _rw_cap = _RW2()
                    _cap_s2n: dict = {}
                    for _cai in sorted(_cluster):
                        _cap_s2n[_cai] = _rw_cap.AddAtom(mol.GetAtomWithIdx(_cai))
                    for _bond in mol.GetBonds():
                        _bai, _eai = _bond.GetBeginAtomIdx(), _bond.GetEndAtomIdx()
                        if _bai in _cap_s2n and _eai in _cap_s2n:
                            _rw_cap.AddBond(
                                _cap_s2n[_bai], _cap_s2n[_eai], _bond.GetBondType()
                            )
                    try:
                        _cap_frag = _C.RemoveHs(_rw_cap.GetMol())
                    except (ValueError, RuntimeError):
                        for _aat in _rw_cap.GetAtoms():
                            _aat.SetIsAromatic(False)
                        for _abd in _rw_cap.GetBonds():
                            if _abd.GetBondTypeAsDouble() == 1.5:
                                _abd.SetBondType(_C.BondType.SINGLE)
                        _cap_frag = _C.RemoveHs(_rw_cap.GetMol(), sanitize=False)
                    if _cap_frag is None or _cap_frag.GetNumAtoms() == 0:
                        continue
                    # Match against library
                    _cap_abbr, _ = _s2c_match_branch_segment(_cap_frag, _sc_full)
                    if not _cap_abbr or _cap_abbr == "?":
                        continue
                    # Skip scaffold monomers — they're handled in step 5.5
                    if _cap_abbr in {a for a, *_ in _s2c_get_scaffold_patterns()}:
                        continue
                    # Determine R-groups at the junction
                    _res_r = _s2c_r_group_at_atom(
                        mol, _ai, _node.mol_atoms, _base_abbr, _sc_raw
                    )
                    _cap_r = _s2c_r_group_at_atom(
                        mol, _ni, _cluster, _cap_abbr, _sc_raw
                    )
                    # Only accept sidechain attachments (R4+).
                    # R1/R2/R3 are backbone N/C connections — those are terminal
                    # caps (ac, am, fmoc, boc) already handled by cap stripping.
                    if _res_r is None or _res_r < 4 or _cap_r is None:
                        continue
                    abbrs[_bi] += f".{_cap_abbr}({_res_r},{_cap_r})"
                    _sc_claimed.update(_cluster)
                    _sc_unclaimed -= _cluster
                    break  # one cap per junction atom; next junction

    # ── 5.5 Scaffold crosslink detection ─────────────────────────────────────
    # Detect multi-arm non-amide crosslink scaffolds (TBMB, etc.) library-
    # driven: any SDF monomer with ≥2 non-backbone dummy attachment points.
    # For each pattern, [n*] in its CHUCKLES was converted to [*:n] SMARTS;
    # matched atoms at those positions must be backbone atoms.  R-group labels
    # are assigned in backbone chain order (!1, !2, …) for deterministic output.
    scaffold_suffix = ""
    _scaffold_xlinks = 0
    _scaffold_patterns = _s2c_get_scaffold_patterns()

    if _scaffold_patterns:
        _rlib = _s2c_get_raw_lib()
        for (
            _sc_abbr,
            _sc_smarts,
            _sc_r_pos,
            _sc_min_slot,
            _sc_arm_types,
            _sc_min_required,
            *_sc_flags,
        ) in _scaffold_patterns:
            _sc_is_partial = bool(_sc_flags[0]) if _sc_flags else False
            _found = False
            for _match in mol.GetSubstructMatches(_sc_smarts, useChirality=False):
                # Collect R-group positions that map to backbone atoms.
                # Allow partial matches: require ≥2 backbone hits so a scaffold
                # with one unreacted arm is still detected.
                _attach: dict = {}  # smarts_idx → (backbone_pos, r_num)
                for _si, _rnum in _sc_r_pos.items():
                    _mat_idx = _match[_si]
                    _bp = next(
                        (
                            pos
                            for pos, node in enumerate(backbone)
                            if _mat_idx in node.mol_atoms
                        ),
                        None,
                    )
                    if _bp is not None:
                        _attach[_si] = (_bp, _rnum)

                if len(_attach) < _sc_min_required:
                    # Partial patterns already encode only the reacted arms in their SMARTS;
                    # every remaining [*:n] slot must hit backbone — no further slack allowed.
                    if _sc_is_partial:
                        continue
                    # Relax when unreacted arms end at an atom consistent with
                    # that arm's chem_type (scaffold-local, not a generic halogen check).
                    # alkyl_halide_c → unreacted = Br/Cl/I leaving group still present
                    # thia_michael_c → unreacted = terminal alkene C (has a C=C double bond)
                    if len(_attach) >= 1:
                        _unmatched_slots = [
                            _si for _si in _sc_r_pos if _si not in _attach
                        ]
                        _accept = True
                        for _usi in _unmatched_slots:
                            _uatom = mol.GetAtomWithIdx(_match[_usi])
                            _ct = _sc_arm_types.get(_sc_r_pos[_usi], "")
                            if _ct == "alkyl_halide_c":
                                if _uatom.GetAtomicNum() not in {35, 17, 53}:
                                    _accept = False
                                    break
                            elif _ct == "thia_michael_c":
                                # Unreacted arm is H-capped by the mol builder (propanoyl
                                # terminus CH3, not vinyl =CH2), so just check for carbon.
                                if _uatom.GetAtomicNum() != 6:
                                    _accept = False
                                    break
                            else:
                                _accept = False
                                break
                        if not _accept:
                            continue

                # Sort attachments by backbone chain order.  Renumber R-groups
                # consecutively from the scaffold's minimum slot (e.g. R4,R5,R6
                # for TBMB) so partial matches (2-of-3 arms) don't produce
                # non-consecutive slot numbers like (4,6) when (4,5) is canonical.
                _sorted_bpos = sorted(_attach.values(), key=lambda x: x[0])
                _sorted_r = list(range(_sc_min_slot, _sc_min_slot + len(_sorted_bpos)))
                # Scaffold XL IDs must not reuse !1 when the backbone ring
                # closure already claims it. Start from 2 for cyclic peptides.
                _xl_id = 2 if cyclic else 1
                _tags = []
                for (_bp, _), _scaffold_r in zip(_sorted_bpos, _sorted_r):
                    _res_r = _s2c_crosslink_r(details[_bp][0], _rlib)
                    _tags.append(f"!{_xl_id}")
                    abbrs[_bp] += f".!{_xl_id}({_res_r},{_scaffold_r})"
                    _xl_id += 1
                scaffold_suffix = f"%{_sc_abbr}." + ".".join(_tags)
                # Track the highest XL ID used so sidechain crosslinks start
                # after it (the sidechain counter uses _scaffold_xlinks + 1).
                _scaffold_xlinks = _xl_id - 1
                _found = True
                break

            if _found:
                break  # one scaffold per molecule

    # ── 5.7 Ring-forming crosslink detection (YAML-driven) ────────────────────
    # Reactions that produce a ring between two sidechain attachment points
    # (e.g. CuAAC 1,4-triazole) cannot be caught by the inter-residue bond
    # scan below — the ring atoms belong to neither residue's mol_atoms.
    # When _raw_xlinks were pre-detected (step 0), the attachment atom indices
    # are already known; we just look them up in the backbone.  Fallback to
    # substructure matching on the original mol for non-de-reactable patterns.
    _ring_xlinks: list = []  # (bi, bj, slot_a, slot_b)
    _ring_xlink_pairs: set = (
        set()
    )  # frozenset({bi, bj}) — to suppress sc_pairs duplication
    if _raw_xlinks:
        for _matom_a, _matom_b, _rslot_a, _rslot_b in _raw_xlinks:
            _rbp_a = next(
                (pos for pos, nd in enumerate(backbone) if _matom_a in nd.mol_atoms),
                None,
            )
            _rbp_b = next(
                (pos for pos, nd in enumerate(backbone) if _matom_b in nd.mol_atoms),
                None,
            )
            if _rbp_a is None or _rbp_b is None or _rbp_a == _rbp_b:
                continue
            _pair_key = frozenset((_rbp_a, _rbp_b))
            if _pair_key in _ring_xlink_pairs:
                continue  # duplicate: same residue pair from another matching pattern
            _ring_xlink_pairs.add(_pair_key)
            _ring_xlinks.append((_rbp_a, _rbp_b, _rslot_a, _rslot_b))
    elif _rxlink_pats:
        # Fallback: no de-reacted backbone was built, try direct substructure match
        if not branch_junctions:
            _rlib = _s2c_get_raw_lib()
        for (
            _rxn_id,
            _prod_mol,
            _ridx_a,
            _ridx_b,
            _rslot_a,
            _rslot_b,
            *_,
        ) in _rxlink_pats:
            for _rmatch in mol.GetSubstructMatches(_prod_mol, useChirality=False):
                _matom_a = _rmatch[_ridx_a]
                _matom_b = _rmatch[_ridx_b]
                _rbp_a = next(
                    (
                        pos
                        for pos, nd in enumerate(backbone)
                        if _matom_a in nd.mol_atoms
                    ),
                    None,
                )
                _rbp_b = next(
                    (
                        pos
                        for pos, nd in enumerate(backbone)
                        if _matom_b in nd.mol_atoms
                    ),
                    None,
                )
                if _rbp_a is None or _rbp_b is None or _rbp_a == _rbp_b:
                    continue
                _ring_xlinks.append((_rbp_a, _rbp_b, _rslot_a, _rslot_b))
                break  # one match per reaction type

    # ── 6. Crosslink detection ────────────────────────────────────────────────
    # Side-chain bonds between pairs of backbone residues that are not the
    # adjacent backbone amide bonds (out_co[i] → in_n[i+1]).
    branch_node_atoms: set = set()
    for jl in branch_junctions.values():
        for jinfo in jl:
            branch_node_atoms.update(jinfo["glu_node"].mol_atoms)

    sc_pairs: list = []  # (bi, bj, atom_in_bi, atom_in_bj)
    backbone_amide_pairs = {
        frozenset((backbone[i].out_co, backbone[i + 1].in_n)) for i in range(n - 1)
    }
    if cyclic:
        backbone_amide_pairs.add(frozenset((backbone[-1].out_co, backbone[0].in_n)))
    for bi in range(n):
        for bj in range(bi + 1, n):
            if frozenset((bi, bj)) in _ring_xlink_pairs:
                continue  # already annotated as ring crosslink; ring atoms aren't
                # direct mol bonds between backbone nodes anyway
            ni_atoms = backbone[bi].mol_atoms
            nj_atoms = backbone[bj].mol_atoms
            for bond in mol.GetBonds():
                ba, ea = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                if ba in ni_atoms and ea in nj_atoms:
                    if frozenset((ba, ea)) not in backbone_amide_pairs:
                        sc_pairs.append((bi, bj, ba, ea))
                        break
                elif ba in nj_atoms and ea in ni_atoms:
                    if frozenset((ba, ea)) not in backbone_amide_pairs:
                        sc_pairs.append((bi, bj, ea, ba))
                        break

    if sc_pairs or _ring_xlinks:
        if not branch_junctions:
            raw_lib = _s2c_get_raw_lib()
        # !1 is reserved for the backbone ring-closure marker when cyclic.
        # Sidechain crosslinks must not reuse it, so start from at least 2.
        xlink_ctr = max(1 if cyclic else 0, _scaffold_xlinks) + 1

        # Ring crosslinks (e.g. CuAAC triazole) — R-group slots from YAML.
        for bi, bj, slot_i, slot_j in _ring_xlinks:
            tag = f"!{xlink_ctr}"
            xlink_ctr += 1
            abbrs[bi] += f".{tag}({slot_i},{slot_j})"
            abbrs[bj] += f".{tag}({slot_j},{slot_i})"

        for bi, bj, atom_i, atom_j in sc_pairs:
            base_i = details[bi][0]
            base_j = details[bj][0]
            # Use backbone slot numbers when the crosslink touches a backbone atom:
            # slot 2 = backbone_c (out_co), slot 1 = backbone_n (in_n).
            # This distinguishes macrolactam amide bonds (sidechain-N to backbone-C=O,
            # giving (4,2)) from sidechain-to-sidechain crosslinks (4,4).
            if atom_i == backbone[bi].out_co:
                r_i = 2
            elif atom_i == backbone[bi].in_n:
                r_i = 1
            else:
                r_i = _s2c_crosslink_r(base_i, raw_lib)
            if atom_j == backbone[bj].out_co:
                r_j = 2
            elif atom_j == backbone[bj].in_n:
                r_j = 1
            else:
                r_j = _s2c_crosslink_r(base_j, raw_lib)
            tag = f"!{xlink_ctr}"
            xlink_ctr += 1
            abbrs[bi] += f".{tag}({r_i},{r_j})"
            abbrs[bj] += f".{tag}({r_j},{r_i})"

    # ── 6.5 Secondary-chain detection (macrocycles bridged via sidechain) ─────
    # Many natural-product peptides (NRPS lipopeptides, polymyxin/colistin
    # family, dual-chain disulfide knots) contain MULTIPLE backbone chains not
    # captured by _s2c_build_backbone's primary pick.  After the primary chain
    # is built, iteratively re-run chain detection on placements outside the
    # already-claimed atoms — each pass claims one more chain until no more
    # eligible candidates remain.  Bridges are then detected pairwise across
    # all chains and emitted as .!N(r_a, r_b) crosslinks.
    extra_chains: list = []  # list of chain lists (each = list of nodes)
    extra_chain_abbrs: list = []  # parallel list of abbrev lists
    # Skip multi-chain detection entirely when the primary chain already has
    # branch annotations (lipid arms via bracket notation), scaffold crosslinks,
    # or ring crosslinks — the branch-arm atoms would otherwise be misread as
    # a phantom second chain (e.g. K's gGlu-AEEA-C20FA arm), regressing 5+
    # GLP-1 drug round-trip tests.
    _skip_multichain = (
        bool(branch_junctions) or bool(_ring_xlinks) or _scaffold_xlinks > 0
    )
    _all_atom_set = set(range(mol.GetNumAtoms())) if not _skip_multichain else set()
    _claimed_atoms = set(backbone_atoms_all)
    _remaining = _all_atom_set - _claimed_atoms
    while len(_remaining) >= 5:
        _candidates = [p for p in placements if p.mol_atoms.issubset(_remaining)]
        if not _candidates:
            break
        try:
            _new_chain = _s2c_build_backbone(_derx_mol, _candidates)
        except (ValueError, RuntimeError):
            _new_chain = []
        if not _new_chain:
            break
        _new_abbrs = []
        for _node in _new_chain:
            _atom_list = list(_node.mol_atoms)
            try:
                _aas, _ = _s2c_isolate_residues(_match_mol, [_atom_list])
                _aa_mol = _aas[0] if _aas else None
            except (ValueError, RuntimeError):
                _aa_mol = None
            if _aa_mol is None:
                _new_abbrs.append(_node.abbr)
                continue
            _bbd = _node.entry.bb_dist
            _lf = [e for e in lib if e[6] == _bbd]
            try:
                _ab, _ = _s2c_match(_aa_mol, _lf or lib)
            except (ValueError, RuntimeError):
                _ab = None
            _new_abbrs.append(_ab or _node.abbr)
        extra_chains.append(_new_chain)
        extra_chain_abbrs.append(_new_abbrs)
        for _node in _new_chain:
            _claimed_atoms.update(_node.mol_atoms)
        _remaining = _all_atom_set - _claimed_atoms

    if extra_chains:
        # Build index maps for ALL chains (primary at index 0, extras 1..N)
        _all_chains = [backbone] + extra_chains
        _all_abbrs = [abbrs] + extra_chain_abbrs  # mutated in-place
        _chain_idx = []  # per-chain {atom_idx: residue_pos}
        _chain_sets = []  # per-chain set of atom indices
        _chain_anchors = []  # per-chain {in_n, out_co} backbone-anchor atoms
        for _ch in _all_chains:
            _idx = {ai: pos for pos, nd in enumerate(_ch) for ai in nd.mol_atoms}
            _chain_idx.append(_idx)
            _chain_sets.append(set(_idx.keys()))
            _chain_anchors.append({n.in_n for n in _ch} | {n.out_co for n in _ch})
        _all_chain_atoms = set().union(*_chain_sets)

        # ── Detect cross-chain bridges (direct + through-linker) ──────────────
        # Same two flavours as the single-secondary-chain case, but iterated
        # across all chain pairs.  Only sidechain atoms qualify as anchors
        # (in_n/out_co would yield invalid .!N(1,2) annotations).
        _bridge_records: list = []  # (i, j, pos_i, pos_j, atom_i, atom_j)
        _seen_global: set = set()

        def _record(i, j, ai, aj):
            ci, cj = (i, j) if i < j else (j, i)
            xi, xj = (ai, aj) if i < j else (aj, ai)
            pi, pj = _chain_idx[ci][xi], _chain_idx[cj][xj]
            k = (ci, cj, pi, pj)
            if k in _seen_global:
                return False
            _seen_global.add(k)
            _bridge_records.append((ci, cj, pi, pj, xi, xj))
            return True

        # Direct-bond bridges: any inter-chain bond
        for _bd in mol.GetBonds():
            _ba, _ea = _bd.GetBeginAtomIdx(), _bd.GetEndAtomIdx()
            _ci = next((i for i, s in enumerate(_chain_sets) if _ba in s), None)
            _cj = next((i for i, s in enumerate(_chain_sets) if _ea in s), None)
            if _ci is None or _cj is None or _ci == _cj:
                continue
            _record(_ci, _cj, _ba, _ea)

        # Through-linker bridges: BFS each chain's sidechain atoms through
        # unplaced atoms (atoms NOT in any chain), recording the first
        # other-chain sidechain reached.  Hop cap 4 covers disulfide (S-S),
        # thioether (S-CH2), and short aliphatic linkers.
        from collections import deque as _dq_b

        _MAX_LINKER_HOPS = 4
        for i in range(len(_all_chains)):
            for _start in sorted(_chain_sets[i] - _chain_anchors[i]):
                if not any(
                    nb.GetIdx() not in _all_chain_atoms
                    for nb in mol.GetAtomWithIdx(_start).GetNeighbors()
                ):
                    continue
                _q = _dq_b([(_start, 0)])
                _bfs_seen = {_start}
                _found = False
                while _q and not _found:
                    _cur, _hops = _q.popleft()
                    if _hops > _MAX_LINKER_HOPS:
                        continue
                    for _nb in mol.GetAtomWithIdx(_cur).GetNeighbors():
                        _ni = _nb.GetIdx()
                        if _ni in _bfs_seen:
                            continue
                        _target = next(
                            (
                                k
                                for k, s in enumerate(_chain_sets)
                                if k != i and _ni in s
                            ),
                            None,
                        )
                        if _target is not None:
                            if _ni in _chain_anchors[_target]:
                                continue
                            if _record(i, _target, _start, _ni):
                                _found = True
                            break
                        if _ni in _all_chain_atoms:
                            continue
                        _bfs_seen.add(_ni)
                        _q.append((_ni, _hops + 1))

        if _bridge_records:
            try:
                raw_lib
            except NameError:
                raw_lib = _s2c_get_raw_lib()
            try:
                xlink_ctr
            except NameError:
                xlink_ctr = max(1 if cyclic else 0, _scaffold_xlinks) + 1
            for _i, _j, _pi, _pj, _ai, _aj in _bridge_records:
                _ch_i, _ch_j = _all_chains[_i], _all_chains[_j]
                _ab_i, _ab_j = _all_abbrs[_i], _all_abbrs[_j]
                if _ai == _ch_i[_pi].out_co:
                    _r1 = 2
                elif _ai == _ch_i[_pi].in_n:
                    _r1 = 1
                else:
                    _r1 = _s2c_crosslink_r(_ab_i[_pi].split(".")[0], raw_lib) or 4
                if _aj == _ch_j[_pj].out_co:
                    _r2 = 2
                elif _aj == _ch_j[_pj].in_n:
                    _r2 = 1
                else:
                    _r2 = _s2c_crosslink_r(_ab_j[_pj].split(".")[0], raw_lib) or 4
                _tag = f"!{xlink_ctr}"
                xlink_ctr += 1
                _ab_i[_pi] += f".{_tag}({_r1},{_r2})"
                _ab_j[_pj] += f".{_tag}({_r2},{_r1})"

    if source_molecule is not None:
        # Only replace undecorated residues. The earlier stages still use the
        # library's names and slot metadata for branches and reaction products.
        # A substituted source template declares only its verified R1/R2 ends.
        references = {entry[0]: entry[1] for entry in lib}
        for position, node in enumerate(backbone):
            abbreviation = details[position][0]
            if abbrs[position] != abbreviation:
                continue
            token = _s2c_local_residue_template(source_molecule, node)
            if token is None:
                continue
            # This is only a local candidate-selection comparison. R1/R2 were
            # just established as N/carbonyl ends, so their free forms are H/OH.
            # The entire emitted peptide is verified by actual assembly later.
            source_residue = _C.MolFromSmiles(_s2c_cap_smiles(token[1:-1], "[H],[OH]"))
            if source_residue is None:
                continue
            if _s2c_preserves_structure(source_residue, references.get(abbreviation)):
                continue
            abbrs[position] = token
            details[position] = (token, 0, source_residue.GetNumAtoms())

    # ── 7. Build CABILN string ────────────────────────────────────────────────
    if cyclic:
        cabiln = "!1-" + "-".join(abbrs) + "-!1"
    else:
        cabiln = "-".join(abbrs)
        if c_cap:
            cabiln += "-" + c_cap
    if scaffold_suffix:
        cabiln += scaffold_suffix
    if n_cap:
        cabiln = n_cap + "-" + cabiln
    for _ec_abbrs in extra_chain_abbrs:
        if _ec_abbrs:
            cabiln += "%" + "-".join(_ec_abbrs)

    return cabiln, details


def _s2c_roundtrip_molecule(cabiln, *, warning_sink=None):
    """Assemble a candidate; invalid chemistry is a rejected candidate."""
    from pyPept.molecule import Molecule
    from pyPept.sequence import Sequence

    if warning_sink is None:
        warning_sink = lambda _message: None
    try:
        return Molecule(Sequence(cabiln, warning_sink=warning_sink)).get_molecule(
            fmt="ROMol"
        )
    except (ValueError, RuntimeError):
        return None


def _s2c_preserves_structure(source, rebuilt):
    """Require the complete graph and every stereo assignment in the input.

    Library monomers may fill unspecified stereo centers with their default
    configuration. This is the converter's existing inference policy; it never
    permits changing or removing an explicitly supplied configuration.
    """
    return compare_structures(source, rebuilt).compatible


def _s2c_one_fragment_candidates(fragment, ring_tag=1, max_breaks=16):
    """Yield full-structure synthetic candidates without guessing residues.

    Try replacing the termini of an amino acid first. Otherwise open an amide:
    a ring becomes one synthetic residue with a terminal crosslink, while an
    acyclic molecule becomes two synthetic caps joined by the original amide.
    Each candidate must still be checked by assembly against the source.
    """
    from rdkit import Chem

    def carbonyl(atom):
        return atom.GetAtomicNum() == 6 and any(
            nb.GetAtomicNum() == 8
            and fragment.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx()).GetBondType()
            == Chem.BondType.DOUBLE
            for nb in atom.GetNeighbors()
        )

    free_nitrogens = [
        atom.GetIdx()
        for atom in fragment.GetAtoms()
        if atom.GetAtomicNum() == 7
        and atom.GetTotalNumHs() >= 1
        and not any(carbonyl(nb) for nb in atom.GetNeighbors())
    ]
    acid_oxygens = [
        atom.GetIdx()
        for atom in fragment.GetAtoms()
        if atom.GetAtomicNum() == 8
        and atom.GetDegree() == 1
        and atom.GetTotalNumHs() == 1
        and carbonyl(atom.GetNeighbors()[0])
    ]
    if free_nitrogens and acid_oxygens:
        # Replace the acid OH; adding another substituent to its carbonyl C
        # would produce pentavalent carbon (the former fallback did exactly that).
        rw = Chem.RWMol(fragment)
        dummy = Chem.Atom(0)
        dummy.SetIsotope(2)
        rw.ReplaceAtom(acid_oxygens[0], dummy)
        n_idx = free_nitrogens[0]
        nitrogen = rw.GetAtomWithIdx(n_idx)
        if nitrogen.GetNumExplicitHs():
            nitrogen.SetNumExplicitHs(nitrogen.GetNumExplicitHs() - 1)
        dummy = Chem.Atom(0)
        dummy.SetIsotope(1)
        d_idx = rw.AddAtom(dummy)
        rw.AddBond(n_idx, d_idx, Chem.BondType.SINGLE)
        try:
            molecule = rw.GetMol()
            Chem.SanitizeMol(molecule)
            yield f"<{Chem.MolToSmiles(molecule)}>", False
        except (ValueError, RuntimeError):
            pass

    n_tried = 0
    for bond in fragment.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetAtomicNum() == 7 and carbonyl(b):
            nitrogen, carbon = a, b
        elif b.GetAtomicNum() == 7 and carbonyl(a):
            nitrogen, carbon = b, a
        else:
            continue
        if n_tried >= max_breaks:
            break
        n_tried += 1
        rw = Chem.RWMol(fragment)
        rw.RemoveBond(nitrogen.GetIdx(), carbon.GetIdx())
        for slot, atom in ((1, nitrogen), (2, carbon)):
            dummy = Chem.Atom(0)
            dummy.SetIsotope(slot)
            d_idx = rw.AddAtom(dummy)
            rw.AddBond(atom.GetIdx(), d_idx, Chem.BondType.SINGLE)
        try:
            opened = rw.GetMol()
            Chem.SanitizeMol(opened)
            parts = Chem.GetMolFrags(opened, asMols=True)
            if len(parts) == 1:
                token = Chem.MolToSmiles(parts[0])
                yield f"!{ring_tag}-<{token}>-!{ring_tag}", True
            elif len(parts) == 2:
                caps = {}
                for part in parts:
                    slot = next(
                        atom.GetIsotope()
                        for atom in part.GetAtoms()
                        if atom.GetAtomicNum() == 0
                    )
                    caps[slot] = Chem.MolToSmiles(part)
                yield f"<{caps[2]}>_-_<{caps[1]}>", False
        except (ValueError, RuntimeError):
            continue


def _s2c_verified_component(molecule):
    """Prefer a verified library decomposition; otherwise use a synthetic one."""
    # These existing normalizers can help recognize a library decomposition,
    # but candidate acceptance and fallback always use the original structure.
    matching_molecule = _s2c_normalize_iminol(_s2c_canonicalize_guanidinium(molecule))
    for preserve_local_residues in (False, True):
        try:
            cabiln, details = _s2c_decompose(
                matching_molecule,
                source_molecule=molecule if preserve_local_residues else None,
            )
        except (ValueError, RuntimeError) as error:
            decomposition_error = error
            continue
        diagnostics = []
        rebuilt = _s2c_roundtrip_molecule(cabiln, warning_sink=diagnostics.append)
        if _s2c_preserves_structure(molecule, rebuilt):
            synthetic_count = sum(name.startswith("<") for name, _, _ in details)
            if synthetic_count:
                noun = "residue" if synthetic_count == 1 else "residues"
                diagnostics.append(
                    f"Preserved {synthetic_count} source {noun} as local synthetic "
                    "CABILN because the library representation would change "
                    "the input structure. Backbone attachment atoms were retained."
                )
            return cabiln, details, rebuilt, bool(synthetic_count), diagnostics
        decomposition_error = ValueError(
            "The monomer decomposition did not preserve the input structure"
        )

    for candidate, _ in _s2c_one_fragment_candidates(molecule):
        diagnostics = []
        rebuilt = _s2c_roundtrip_molecule(candidate, warning_sink=diagnostics.append)
        if _s2c_preserves_structure(molecule, rebuilt):
            diagnostics.append(
                "Could not identify every library residue; preserved the component "
                "as coarse synthetic CABILN. Residue match details are unavailable; "
                "the synthetic slots do not establish peptide residue boundaries."
            )
            return candidate, [], rebuilt, True, diagnostics
    raise ValueError(
        "Cannot represent this SMILES component in CABILN without changing its "
        "structure; register the missing monomer or use another input structure"
    ) from decomposition_error


def convert_smiles(smiles: str) -> ConversionResult:
    """Convert SMILES with explicit status and no Python warning side effects.

    Every connected component is retained. The assembled result must preserve
    connectivity, charge, isotopes and all specified stereochemistry. Known
    monomers supply default configurations for unspecified stereo centers.
    Alternative tautomer forms may be tried during matching, but the returned
    notation must preserve the input form.
    Local synthetic residues retain actual backbone attachment atoms. Coarse
    synthetic fallbacks preserve structure but have no residue match details.
    Inputs that cannot be represented raise ``ValueError``.
    """
    from rdkit import Chem

    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("SMILES must be a non-empty string")
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or not molecule.GetNumAtoms():
        raise ValueError(f"Invalid SMILES: {smiles[:80]}")
    require_supported_stereo(molecule)
    fragments = Chem.GetMolFrags(molecule, asMols=True)
    parts, details, diagnostics, synthetic_components = [], [], [], []
    next_bond_id = 1
    for component_index, fragment in enumerate(fragments):
        part, fragment_details, rebuilt, synthetic, messages = _s2c_verified_component(
            fragment
        )
        if synthetic:
            synthetic_components.append(component_index)
        diagnostics.extend(
            (
                f"Component {component_index + 1}: {message}"
                if len(fragments) > 1
                else message
            )
            for message in messages
        )
        if len(fragments) == 1:
            parts.append(part)
            details.extend(fragment_details)
            break
        # Each independently converted component starts crosslink numbering at
        # one. A global namespace prevents accidental bonds between components.
        bond_ids = {}

        def renumber(match):
            nonlocal next_bond_id
            key = match.group(0)
            if key not in bond_ids:
                bond_ids[key] = f"!{next_bond_id}"
                next_bond_id += 1
            return bond_ids[key]

        parts.append(_re.sub(r"![A-Za-z0-9_]+", renumber, part))
        details.extend(fragment_details)
    cabiln = "%".join(parts)
    if len(parts) > 1:
        rebuilt = _s2c_roundtrip_molecule(cabiln)
    comparison = compare_structures(molecule, rebuilt)
    if not comparison.compatible:
        raise ValueError("Cannot preserve all SMILES components in one CABILN string")
    inferred_stereo = not comparison.exact
    if inferred_stereo:
        diagnostics.append(_STEREO_INFERENCE_MESSAGE)
    return ConversionResult(
        cabiln=cabiln,
        details=details,
        inferred_stereo=inferred_stereo,
        synthetic_components=tuple(synthetic_components),
        warnings=tuple(dict.fromkeys(diagnostics)),
    )


def smiles_to_cabiln_core(smiles: str):
    """Compatibility API returning ``(notation, details)`` and emitting warnings.

    Applications should use :func:`convert_smiles` for request-local status.
    """
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
