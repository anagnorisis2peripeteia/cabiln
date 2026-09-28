"""
A class to parse and manipulate BILN sequences.

From publication: pyPept: a python library to generate atomistic 2D and 3D representations of peptides
Journal of Cheminformatics, 2023
"""

########################################################################################
# Authorship
########################################################################################

__credits__ = ["Rodrigo Ochoa", "J.B. Brown", "Thomas Fox"]
__license__ = "MIT"


########################################################################################
# Modules
########################################################################################
# pylint: disable=E1101

# System libraries
import sys
import copy
import re
import os
import warnings

from pyPept.source import (
    SourceText, bond_marker as _source_bond_marker, group as _source_group,
    join as _source_join, record as _source_record, sub as _source_sub,
    synthetic as _source_synthetic,
)

import string
from collections import defaultdict
from dataclasses import dataclass, field
from importlib.resources import files

# Third-party libraries
import hashlib
import numpy as np
from rdkit import Chem


##########################################################################
# Functions and classes
##########################################################################

# Regex for inline synthetic tokens:
#   `<sidechain>`       — alpha-AA template (sidechain only or full [1*]/[2*] form)
#   `<smi>_`            — N-cap (single R2 slot connecting to next residue's R1 via amide)
#   `_<smi>`            — C-cap (single R1 slot connecting to previous residue's R2)
# Order matters in the regex pre-replacer: caps must be tried BEFORE plain
# alpha-AA tokens so the trailing/leading `_` is consumed with the bracket.
_SYN_NCAP_RE = re.compile(r'<([^<>]*)>_')
_SYN_CCAP_RE = re.compile(r'_<([^<>]*)>')
_SYN_AA_RE   = re.compile(r'<([^<>]*)>')


def _build_synthetic_aa(sidechain_smi):
    """Build an alpha-AA ROMol template from a sidechain SMILES (or a full
    residue template).

    Two input forms accepted:
      1. Sidechain-only (no `[1*]` / `[2*]`): wrapped as
         `[1*]NC(<sc>)C(=O)[2*]`.  Cα stereo is unset.
      2. Full template (contains both `[1*]` and `[2*]`): used as-is.
         This lets the walker preserve Cα stereo by emitting e.g.
         `<[1*]N[C@@H](C)C([2*])=O>` for L-alanine-like residues.

    Layout convention: R1 (isotope 1) on backbone N, R2 (isotope 2) on
    carbonyl C, matching the SDF library.  Returns None on parse error.
    Empty / whitespace `sidechain_smi` produces a glycine equivalent.
    """
    import re as _re
    if sidechain_smi is None:
        sidechain_smi = ''
    sc = sidechain_smi.strip()
    if not sc:
        smi = '[1*]NCC(=O)[2*]'
    elif _re.search(r'\[1\*[:\]]', sc) and _re.search(r'\[2\*[:\]]', sc):
        # Full template form: accept [1*], [1*:n], [2*], [2*:n]. Atom-map digits
        # after the isotope are in-band flags (e.g. [2*:99] = aldehyde terminus).
        smi = sc
    else:
        smi = f'[1*]NC({sc})C(=O)[2*]'
    try:
        return Chem.MolFromSmiles(smi)
    except Exception:
        return None


def _build_synthetic_cap(cap_smi, side):
    """Build a synthetic cap ROMol from a SMILES.

    side='n':  N-cap — the cap's C(=O) bonds to the next residue's R1
               via standard backbone_amide.  Output has R2 only.
    side='c':  C-cap — the cap's N bonds to the previous residue's R2.
               Output has R1 only.

    Two input forms accepted (same as _build_synthetic_aa):
      - Bare SMILES (no `[1*]` / `[2*]`): walker is responsible for
        placing the dummy atom; we wrap minimally.  For N-cap, append
        `[2*]` to the LAST atom; for C-cap, prepend `[1*]` to the first.
      - Full template containing the appropriate `[N*]`: used as-is.
    """
    if cap_smi is None:
        cap_smi = ''
    s = cap_smi.strip()
    if not s:
        return None
    if side == 'n':
        if '[2*]' in s:
            smi = s
        else:
            smi = f'{s}[2*]'
    elif side == 'c':
        if '[1*]' in s:
            smi = s
        else:
            smi = f'[1*]{s}'
    else:
        return None
    try:
        return Chem.MolFromSmiles(smi)
    except Exception:
        return None


def _synthetic_aa_symbol(sidechain_smi):
    """Stable synthetic-AA monomer symbol from sidechain SMILES."""
    h = hashlib.md5((sidechain_smi or '').encode()).hexdigest()[:10]
    return f'__SYN_{h}'


def _synthetic_cap_symbol(cap_smi, side):
    """Stable synthetic cap monomer symbol from cap SMILES + side ('n' / 'c')."""
    h = hashlib.md5(f'{side}:{cap_smi or ""}'.encode()).hexdigest()[:10]
    return f'__SYN{side.upper()}CAP_{h}'


def _infer_synth_chem_type(mol, slot):
    """Return (chem_type, leaving_group) for an R-group dummy in a synthetic
    monomer template.  Slot 1 → backbone_n, slot 2 → backbone_c, slot 3+ →
    inferred from the dummy's neighbour atom (sidechain anchor).

    The default backbone-C leaving group is [OH] (carboxylic acid terminus).
    To emit a peptide-aldehyde terminus instead (C(=O)H), the caller writes the
    slot-2 dummy as [2*:99] in SMILES — atom-map 99 is an in-band marker that
    overrides _lg to [H] without changing the CABILN grammar."""
    if slot == 1:
        return ('backbone_n', '[H]')
    if slot == 2:
        for atom in mol.GetAtoms():
            if (atom.GetAtomicNum() == 0 and atom.GetIsotope() == 2
                    and atom.GetAtomMapNum() == 99):
                return ('backbone_c', '[H]')
        return ('backbone_c', '[OH]')
    # Find the dummy atom and its neighbour
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == slot:
            nbs = atom.GetNeighbors()
            if not nbs:
                return (None, None)
            anchor = nbs[0]
            sym = anchor.GetSymbol()
            # Quick chem_type inference by neighbour atom type and context
            if sym == 'S':
                return ('thiol', '[H]')
            if sym == 'Se':
                return ('selenol', '[H]')
            if sym == 'N':
                # Could be amine or guanidinium etc. Default to amine_primary.
                return ('amine_primary', '[H]')
            if sym == 'O':
                # Phenolic O if anchor is bonded to aromatic C; else aliphatic OH
                for nb2 in anchor.GetNeighbors():
                    if nb2.GetIsAromatic():
                        return ('aryl_phenol_o', '[H]')
                return ('hydroxyl', '[H]')
            if sym == 'C':
                # carbonyl C → carboxyl; alkyl halide neighbour; alkene; etc.
                # Look for =O neighbour → carboxyl
                for nb2 in anchor.GetNeighbors():
                    if nb2.GetAtomicNum() == 8:
                        b = mol.GetBondBetweenAtoms(anchor.GetIdx(), nb2.GetIdx())
                        if b and b.GetBondTypeAsDouble() == 2.0:
                            return ('carboxyl', '[OH]')
                # Alkyl halide
                for nb2 in anchor.GetNeighbors():
                    if nb2.GetSymbol() in ('Cl', 'Br', 'I'):
                        return ('alkyl_halide_c', None)
                # Default carbon anchor
                return ('carbon', '[H]')
            return (None, None)
    return (None, None)

def _attachment_idx(mol, slot):
    """Return attachment atom index for R-group slot (1-based), or None.

    Locates the dummy atom whose isotope equals slot (CHUCKLES convention:
    slot 1 → [1*], slot 2 → [2*], …) and returns its neighbour's index —
    i.e. the heavy atom where the inter-monomer bond will form.
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == slot:
            nb = atom.GetNeighbors()
            return nb[0].GetIdx() if nb else None
    return None


def _rgroup_atom_idx(mol, slot):
    """Return the dummy atom index for R-group slot (1-based), or None."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and atom.GetIsotope() == slot:
            return atom.GetIdx()
    return None


def _slot_for_attachment(mol, atom_idx):
    """Return R-group slot (1-based) whose dummy neighbours atom_idx, or None."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            for nb in atom.GetNeighbors():
                if nb.GetIdx() == atom_idx:
                    return atom.GetIsotope()
    return None


# Matches .!marker with optional (host_r,cap_r) — crosslink endpoint.
# Parens required on first occurrence; may be omitted on second (inverse inferred).
_INLINE_BOND_RE = re.compile(r'\.(!\w+)(?:\((\d+),(\d+)\))?')

# Matches .CapToken(host_r,cap_r) — named cap attachment.  The token may start
# with a letter OR underscore: the monomer library uses `_OMe`, `_Bn`, `_NHBn`,
# etc. for sidechain fragments whose token deliberately begins with `_` to
# signal "non-canonical residue".  Without the underscore in the leading-char
# class, ._OMe(4,1) falls through to the crosslink-branch handler in
# Sequence._read_bonds and trips the odd-bond parity check (28+ rt_fail cases
# in cyclicpepedia, 2026-05-16).
_INLINE_CAP_RE  = re.compile(r'\.([A-Za-z_]\w*)\((\d+),(\d+)\)')

# Detects old BILN bare-integer crosslink annotations: Token(bid,rg) not preceded by '.'.
_OLD_BILN_RE = re.compile(r'(?<![.!\w\[{])([A-Za-z]\w*)\((\d+),(\d+)\)')

# Matches .[...] sequential reaction bracket on a residue: .[Cap1(x,y).Cap2(z,w)...]
# The dot sits outside the bracket (consistent with inline .Cap(r,r) notation).
# Each entry inside uses Fragment(prev_r, cap_r); entries separated by '.'.
# Sub-bracket arms [.Entry(r,r)] inside .[...] are supported: each arm connects from
# the current chain pointer without advancing it, enabling multi-arm hubs.
_BRACKET_RE = re.compile(
    r'\.\[([^\[\]]*(?:\[[^\[\]]*\])*[^\[\]]*)\]'  # .[...] with optional inner [.arm]s
    r'|\.\{([^{}]*)\}'                              # .{...} legacy BILN form
)


def _flatten_nested_brackets(seg):
    """Convert old ``[[A.B].C]`` nested notation to new ``[A[.B][.C]]`` form.

    Old notation used nesting to express multi-arm hubs: ``[[TBMB.C].!3]``
    meant both C and !3 bond to TBMB.  New notation uses explicit sub-bracket
    arms: ``[TBMB[.C][.!3]]``.  This function converts the old form on the way
    into the parser so both notations assemble identically.
    """
    out = []
    i = 0
    while i < len(seg):
        if seg[i:i + 2] == '.[' and i + 2 < len(seg) and seg[i + 2] == '[':
            depth = 0
            j = i + 1
            while j < len(seg):
                if seg[j] == '[':
                    depth += 1
                elif seg[j] == ']':
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            whole = seg[i + 1:j + 1]
            flat = _flatten_one_nested(whole)
            out.append(seg[i:i+1] + flat)
            i = j + 1
        else:
            out.append(seg[i])
            i += 1
    return _source_join('', out)


_ENTRY_ANY = re.compile(r'((?:[A-Za-z]\w*|!\w+)\(\d+,\d+\))')


def _flatten_one_nested(s):
    """Convert a single nested bracket to new sub-bracket arm notation.

    Input ``s`` is the bracket expression starting at the outer ``[``.
    Returns ``[anchor[.inner_rest...][.outer...]]`` — every arm from the
    anchor is an explicit ``[.Entry(r,r)]`` sub-bracket so the parser's
    pointer semantics are unambiguous.
    """
    inner_start = s.index('[', 1)
    depth = 0
    inner_end = inner_start
    while inner_end < len(s):
        if s[inner_end] == '[':
            depth += 1
        elif s[inner_end] == ']':
            depth -= 1
            if depth == 0:
                break
        inner_end += 1

    inner_bracket = s[inner_start:inner_end + 1]
    inner_content = inner_bracket[1:-1]
    outer_part = s[inner_end + 1:-1]
    outer_entries = [_source_group(m, 1) for m in _ENTRY_ANY.finditer(outer_part)]

    # If inner content already uses sub-bracket notation, preserve it and just
    # append outer entries as additional arms.
    if '[' in inner_content:
        outer_arms = _source_join('', ('[.' + e + ']' for e in outer_entries))
        return s[:1] + inner_content + outer_arms + s[-1:]

    # Legacy flat inner content: split anchor from remaining chain entries.
    inner_entries = [_source_group(m, 1) for m in _ENTRY_ANY.finditer(inner_content)]
    anchor = inner_entries[0] if inner_entries else ''
    inner_rest = inner_entries[1:]

    inner_arms = _source_join('', ('[.' + e + ']' for e in inner_rest))
    outer_arms = _source_join('', ('[.' + e + ']' for e in outer_entries))
    return s[:1] + anchor + inner_arms + outer_arms + s[-1:]

# Matches entries inside brackets: both monomer .Token(r,r) and crosslink .!n(r,r).
_BRACKET_ENTRY_RE = re.compile(r'\.((?:[A-Za-z]\w*|!\w+))\((\d+),(\d+)\)')

# Mixed tokeniser for bracket content: matches either a sub-bracket arm [.sub] or a
# flat entry .Entry(r,r).  Used by _sub_bracket to process new-grammar notation.
#   group(1) set  → sub-bracket arm content (without surrounding [. and ])
#   group(2,3,4)  → flat entry tok, prev_r, cur_r
_MIXED_ENTRY_RE = re.compile(
    r'\[\.([^\[\]]*)\]'
    r'|\.((?:[A-Za-z]\w*|!\w+))\((\d+),(\d+)\)'
)


def biln_to_cabiln(biln):
    """Convert old BILN crosslink notation ``Token(bid,rg)`` to CABILN ``.!n(y,z)`` form.

    Old BILN embeds bond IDs as bare integers inside residue names::

        C(1,3)-A-A-A-C(1,3)   # bond 1, R3 on each Cys

    CABILN uses inline dot notation::

        C.!1(4,4)-A-A-A-C.!1  # slot 3 in HELM/old-BILN = slot 4 in pyPept (sidechain)

    Old BILN/HELM R3 means "first sidechain".  pyPept inserts backbone_n_mod at
    slot 3, pushing sidechains to slot 4+, so R3 maps to pyPept slot 4.

    :param biln: BILN string, possibly containing old crosslink annotations.
    :returns: equivalent CABILN string.
    """
    def _remap_slot(rg_str):
        """HELM/BILN R3 (sidechain) → pyPept slot 4 (backbone_n_mod occupies slot 3)."""
        return '4' if rg_str == '3' else rg_str

    matches = list(_OLD_BILN_RE.finditer(biln))
    if not matches:
        return biln

    bond_groups = {}
    for m in matches:
        bid = m.group(2)
        bond_groups.setdefault(bid, []).append(m)

    replacements = {}
    for bid, endpoints in bond_groups.items():
        if len(endpoints) != 2:
            continue
        m1, m2 = endpoints
        tok1, rg1 = m1.group(1), _remap_slot(m1.group(3))
        tok2, rg2 = m2.group(1), _remap_slot(m2.group(3))
        replacements[m1.start()] = (m1, f'{tok1}.!{bid}({rg1},{rg2})')
        replacements[m2.start()] = (m2, f'{tok2}.!{bid}')

    result = biln
    for pos in sorted(replacements, reverse=True):
        m, replacement = replacements[pos]
        result = result[:m.start()] + replacement + result[m.end():]
    return result


def _preprocess_cabiln(biln):
    """Normalise newline- or %-separated CABILN into a single %-delimited string."""
    biln = biln.strip()
    biln = _source_sub(r'[ \t]*[\n%][ \t]*', '%', biln)
    biln = _source_sub(r'%+', '%', biln)
    return biln.strip('%')


def _handle_terminal_bond_markers(seg, seen, count):
    """Replace bare !n at chain start/end with explicit bond annotation.

    !n-A-B-C  → A(!n,1)-B-C   (N-terminal: attaches via first residue's R1)
    A-B-C-!n  → A-B-C(!n,2)   (C-terminal: attaches via last residue's R2)

    Validates the implied R-group against seen[tok][1] (partner rgroup from the
    first explicit .!n(y,z) occurrence).  Raises ValueError on mismatch or > 2
    endpoints for the same bond ID.
    """
    parts = seg.split('-')

    def _is_bare_marker(tok):
        return bool(re.match(r'^!\w+$', tok)) and '(' not in tok

    def _register(tok, rgroup, position_label):
        if tok not in seen:
            # Terminal marker with no explicit .!n(y,z) anywhere — self-register.
            # N-terminal always uses R1 and C-terminal always uses R2, so the
            # partner rgroup is unambiguous from position alone.
            partner_rgroup = 2 if rgroup == 1 else 1
            seen[tok] = (str(rgroup), str(partner_rgroup))
        else:
            partner = seen[tok][1]
            if partner != str(rgroup):
                raise ValueError(
                    f"Branch {position_label} marker {tok!r} implies R{rgroup} "
                    f"attachment, but .{tok}(y,z) declared partner rgroup "
                    f"R{partner} (not R{rgroup}). Fix R-group assignments or reorder.")
        if count.get(tok, 0) >= 2:
            raise ValueError(f"Bond {tok!r}: 3rd endpoint found (max 2 allowed).")
        count[tok] = count.get(tok, 0) + 1

    # N-terminal: first part is bare !n
    if len(parts) >= 2 and _is_bare_marker(parts[0]):
        tok = parts[0]
        _register(tok, 1, 'N-terminal (R1)')
        _source_bond_marker(tok, 1, owner=parts[1], terminal=True)
        parts[1] = parts[1] + f'({tok},1)'
        del parts[0]

    # C-terminal: last part is bare !n (checked after N-terminal, handles both ends)
    if len(parts) >= 2 and _is_bare_marker(parts[-1]):
        tok = parts[-1]
        _register(tok, 2, 'C-terminal (R2)')
        _source_bond_marker(tok, 2, owner=parts[-2], terminal=True)
        parts[-2] = parts[-2] + f'({tok},2)'
        del parts[-1]

    return _source_join('-', parts)


def _emit_warning(message, category=UserWarning, stacklevel=1, *, warning_sink=None):
    """Keep request diagnostics local while retaining the library warning API."""
    if warning_sink is None:
        warnings.warn(message, category, stacklevel=stacklevel + 1)
    else:
        warning_sink(str(message))


def _expand_inline_caps(biln, peptide_branch_threshold=2, warning_sink=None):
    """Pre-process CABILN notation before chain/residue splitting.

    Supported forms:
    - '%' or newline as segment separator (main chain first, branches follow)
    - .boc(3,1)            named cap — auto-assigns bond ID >= 100, appended as pendant
    - .!1(4,2)             crosslink first endpoint — explicit R-groups, defines inverse
    - .!1                  crosslink second endpoint — parens omitted, inverse inferred
    - !1-A-B-C             branch N-terminal marker — attaches via R1 of first residue
    - A-B-C-!1             branch C-terminal marker — attaches via R2 of last residue
    - [.A(r,s).B(t,u)...]  sequential reaction bracket — left-to-right chain of named
                            cap attachments.  For each step after the first, the first
                            R-group number refers to the PRECEDING fragment (not the host
                            residue).  The host receives only the first bond annotation.
                            Auto bond IDs are assigned in left-to-right order.

    Returns (expanded_biln, branch_rgroup).  expanded_biln uses '.' as chain
    separator; branch_rgroup maps each !x to its partner rgroup from first occurrence.
    """
    biln = _preprocess_cabiln(biln)
    segments = [s for s in biln.split('%') if s.strip()]

    # Preliminary scan: collect _seen from all explicit .!n(y,z) occurrences so
    # terminal markers and no-parens second occurrences can look up the inverse
    # regardless of segment order.
    _seen = {}
    _SCAN_RE = re.compile(r'\.(!\w+)\((\d+),(\d+)\)')
    for seg in segments:
        for m in _SCAN_RE.finditer(seg):
            tok, hr, cr = m.group(1), m.group(2), m.group(3)
            if tok not in _seen:
                _seen[tok] = (hr, cr)
            else:
                fhr, fcr = _seen[tok]
                if hr != fcr or cr != fhr:
                    raise ValueError(
                        f"Bond {tok!r}: R-group conflict. "
                        f"First endpoint declared R{fhr}→partner R{fcr}; "
                        f"second endpoint declares R{hr}→partner R{cr} "
                        f"(expected inverse R{fcr}→partner R{fhr}).")

    appended = []
    _ctr = [100]
    _count = {}

    def _sub_bond(m):
        tok, hr, cr = m.group(1), m.group(2), m.group(3)
        if _count.get(tok, 0) >= 2:
            raise ValueError(
                f"Bond {tok!r}: 3rd endpoint found (max 2 allowed).")
        if hr is None:
            # No parens — second occurrence; infer inverse from first
            if tok not in _seen:
                raise ValueError(
                    f"Bond marker {tok!r} used without parens but no prior "
                    f".{tok}(y,z) first-occurrence found.")
            fhr, fcr = _seen[tok]
            hr, cr = fcr, fhr
        _count[tok] = _count.get(tok, 0) + 1
        _source_bond_marker(_source_group(m), int(hr))
        return f'({tok},{hr})'

    def _sub_cap(m):
        tok, hr, cr = _source_group(m, 1), m.group(2), m.group(3)
        _source_record(tok, _source_group(m), 'inline')
        bid = _ctr[0]; _ctr[0] += 1
        appended.append(tok + f'({bid},{cr})')
        return f'({bid},{hr})'

    # Shared container so _sub_bracket can read surrounding segment context
    # set just before each call: (host_name, prefix_before_host, rest_of_seg)
    _bracket_ctx = [None]

    def _sub_bracket(m):
        # dot lives outside bracket in notation; restore for _MIXED_ENTRY_RE parsing
        raw = _source_group(m, 1) if m.group(1) is not None else _source_group(m, 2)
        content = '.' + raw
        tokens = list(_MIXED_ENTRY_RE.finditer(content))
        if not tokens:
            raise ValueError(
                f"Sequential bracket {m.group(0)!r} contains no valid "
                f".Fragment(host_r,cap_r) entries.")
        reconstructed = ''.join(t.group(0) for t in tokens)
        if reconstructed != content:
            raise ValueError(
                f"Sequential bracket {m.group(0)!r} has unrecognised content; "
                f"expected only .Fragment(r,r), .!n(r,r), or [.arm(r,r)] entries.")

        first = tokens[0]
        if first.group(1) is not None:
            raise ValueError(
                f"Sequential bracket {m.group(0)!r}: first entry must be a flat "
                f".Entry(r,r) or .!n(r,r), not a sub-bracket arm.")

        tok1, hr1, cr1 = _source_group(first, 2), first.group(3), first.group(4)
        _source_record(tok1, _source_group(first), 'bracket', _source_group(m), terminal=len(tokens) == 1)

        if tok1.startswith('!'):
            # Pure crosslink bracket: .[!1(4,4)] or .[!1(4,4).!2(5,3)]
            # Each entry is a crosslink bond on the host monomer — no pendant fragment.
            host_bonds = []
            for t in tokens:
                if t.group(1) is not None:
                    raise ValueError(
                        f"Sequential bracket {m.group(0)!r}: sub-bracket arm cannot "
                        f"appear in a pure-crosslink bracket.")
                tok = t.group(2)
                if not tok.startswith('!'):
                    raise ValueError(
                        f"Sequential bracket {m.group(0)!r}: monomer entry {tok!r} "
                        f"cannot follow crosslink-first entries.")
                prev_r, cur_r = t.group(3), t.group(4)
                if tok not in _seen:
                    _seen[tok] = (prev_r, cur_r)
                _count[tok] = _count.get(tok, 0) + 1
                _source_bond_marker(_source_group(t), int(prev_r), bracket=_source_group(m))
                host_bonds.append(f'({tok},{prev_r})')
            return ''.join(host_bonds)

        bid1 = _ctr[0]; _ctr[0] += 1
        frag_tokens = [tok1]
        frag_bond_parts = [f'({bid1},{cr1})']
        _backbone_steps = []
        pointer_idx = 0  # index into frag_bond_parts of the current chain tail

        for token_position, t in enumerate(tokens[1:], 1):
            if t.group(1) is not None:
                # Sub-bracket arm [.sub_content]: process its entries as a chain
                # branching FROM the current pointer.  The arm does NOT advance
                # pointer_idx — the next sibling arm or flat entry still bonds from
                # the same host fragment.
                sub_content = '.' + _source_group(t, 1)
                sub_steps = list(_BRACKET_ENTRY_RE.finditer(sub_content))
                if not sub_steps or ''.join(step.group(0) for step in sub_steps) != sub_content:
                    raise ValueError(
                        f"Sequential sub-bracket {t.group(0)!r} has unrecognised "
                        "content; expected only .Fragment(r,r) or .!n(r,r) entries.")
                sub_ptr = pointer_idx
                for step_position, step in enumerate(sub_steps):
                    tok, prev_r, cur_r = _source_group(step, 1), step.group(2), step.group(3)
                    _source_record(tok, _source_group(step), 'bracket', _source_group(m), _source_group(t), step_position == len(sub_steps) - 1)
                    if tok.startswith('!'):
                        _source_bond_marker(_source_group(step), int(prev_r),
                                            owner=frag_tokens[sub_ptr], bracket=_source_group(m), arm=_source_group(t))
                        frag_bond_parts[sub_ptr] += f'({tok},{prev_r})'
                        if tok not in _seen:
                            _seen[tok] = (prev_r, cur_r)
                        _count[tok] = _count.get(tok, 0) + 1
                    else:
                        bid = _ctr[0]; _ctr[0] += 1
                        frag_bond_parts[sub_ptr] += f'({bid},{prev_r})'
                        frag_tokens.append(tok)
                        frag_bond_parts.append(f'({bid},{cur_r})')
                        sub_ptr = len(frag_bond_parts) - 1
                # pointer_idx intentionally NOT updated: arm is a branch, not chain
            else:
                # Flat entry: crosslink annotates current host; monomer advances chain.
                tok, prev_r, cur_r = _source_group(t, 2), t.group(3), t.group(4)
                _source_record(tok, _source_group(t), 'bracket', _source_group(m), terminal=token_position == len(tokens) - 1)
                if tok.startswith('!'):
                    _source_bond_marker(_source_group(t), int(prev_r),
                                        owner=frag_tokens[pointer_idx], bracket=_source_group(m))
                    frag_bond_parts[pointer_idx] += f'({tok},{prev_r})'
                    if tok not in _seen:
                        _seen[tok] = (prev_r, cur_r)
                    _count[tok] = _count.get(tok, 0) + 1
                else:
                    bid = _ctr[0]; _ctr[0] += 1
                    if prev_r == '2' and cur_r == '1':
                        _backbone_steps.append(tok)
                    frag_bond_parts[pointer_idx] += f'({bid},{prev_r})'
                    frag_tokens.append(tok)
                    frag_bond_parts.append(f'({bid},{cur_r})')
                    pointer_idx = len(frag_bond_parts) - 1

        if (peptide_branch_threshold is not None
                and len(_backbone_steps) >= peptide_branch_threshold):
            flat_chain_toks = [
                t.group(2) for t in tokens
                if t.group(1) is None and not t.group(2).startswith('!')
            ]
            frag_chain = '-'.join(flat_chain_toks)
            branch_seg = f'!n-{frag_chain}' if cr1 == '1' else f'{frag_chain}-!n'
            ctx = _bracket_ctx[0]
            if ctx:
                host_name, pre, rest = ctx
                suggestion = f"{pre}{host_name}.!n({hr1},{cr1}){rest}%%{branch_seg}"
            else:
                suggestion = f"host.!n({hr1},{cr1})-[chain]%%{branch_seg}"
            _emit_warning(
                f"Sequential bracket {m.group(0)!r} contains {len(_backbone_steps)} "
                f"R2->R1 connections: backbone amide pattern detected — this creates "
                f"a peptide branch inside [...].  Did you mean:\n"
                f"  {suggestion}",
                UserWarning, stacklevel=6, warning_sink=warning_sink)

        for tok, bonds in zip(frag_tokens, frag_bond_parts):
            appended.append(tok + bonds)

        return f'({bid1},{hr1})'

    processed_segments = []
    for segment_number, seg in enumerate(segments):
        if isinstance(seg, SourceText):
            seg.tracker.segment = segment_number
            for root in split_outside(seg, '-', '[]'):
                seg.tracker.root(root, segment_number)
        seg = _handle_terminal_bond_markers(seg, _seen, _count)
        seg = _flatten_nested_brackets(seg)
        # Brackets first: they may contain .!n(r,r) crosslink entries that
        # _INLINE_BOND_RE would otherwise grab prematurely.
        parts = []
        last_end = 0
        for bm in _BRACKET_RE.finditer(seg):
            prefix = seg[last_end:bm.start()]
            parts.append(prefix)
            if '-' in prefix:
                pre_host, host_token = prefix.rsplit('-', 1)
                pre_host += '-'
            else:
                pre_host, host_token = '', prefix
            host_name = re.split(r'[(\[{]', host_token)[0].strip()
            _bracket_ctx[0] = (host_name, pre_host, seg[bm.end():])
            parts.append(_sub_bracket(bm))
            last_end = bm.end()
        parts.append(seg[last_end:])
        seg = _source_join('', parts)
        seg = _source_sub(_INLINE_BOND_RE, _sub_bond, seg)
        seg = _source_sub(_INLINE_CAP_RE, _sub_cap, seg)
        processed_segments.append(seg)

    all_bond_ids = set(_seen) | set(_count)
    for tok in all_bond_ids:
        cnt = _count.get(tok, 0)
        if cnt != 2:
            raise ValueError(
                f"Bond {tok!r}: {cnt} endpoint(s) declared; expected exactly 2.")

    result = _source_join('.', processed_segments)
    if appended:
        result += '.' + _source_join('.', appended)
    if isinstance(result, SourceText):
        result.tracker.labels = set(_seen) | set(_count)
    branch_rgroup = {tok: int(cr) for tok, (hr, cr) in _seen.items()}
    return result, branch_rgroup


# ANSI colour palette for bracket pairs — cycles like Excel's rainbow parentheses.
_BRACKET_COLOURS = [
    '\033[93m',  # yellow
    '\033[96m',  # cyan
    '\033[92m',  # green
    '\033[95m',  # magenta
    '\033[91m',  # red
]
_ANSI_RESET = '\033[0m'


def colorize_cabiln(biln, use_ansi=True):
    """Return a CABILN string with matching [...] bracket pairs colour-coded.

    Each bracket pair gets a distinct ANSI colour (cycling through 5 colours),
    so matching open/close brackets are visually obvious at a glance — the same
    visual aid Excel provides for nested parentheses.

    :param biln: a CABILN notation string.
    :param use_ansi: if False, return the plain string unmodified (useful for
                     environments that don't support ANSI escape codes).
    :return: str with ANSI colour codes injected around [...] pairs.
    """
    if not use_ansi:
        return biln
    result = []
    bracket_count = 0
    in_bracket = False
    for ch in biln:
        if ch in ('[', '{'):
            colour = _BRACKET_COLOURS[bracket_count % len(_BRACKET_COLOURS)]
            result.append(colour + ch)
            in_bracket = True
        elif ch in (']', '}'):
            result.append(ch + _ANSI_RESET)
            bracket_count += 1
            in_bracket = False
        else:
            result.append(ch)
    return ''.join(result)


import re as _re


def _parse_bracket_items(bracket_content):
    """Parse bracket content into list of (abbr, r_prev, r_this) tuples."""
    parts = bracket_content.split('.')
    items = []
    for p in parts:
        pm = _re.fullmatch(r'([^(]+)\((\d+),(\d+)\)', p)
        if pm:
            items.append((pm.group(1), pm.group(2), pm.group(3)))
        else:
            items.append((p, None, None))
    return items


def cabiln_to_branch(cabiln):
    """Convert CABILN bracket notation to branch (%) notation.

    All (2,1) continuations — positional branch, anchor at N-term::

        D.[G(4,1).A(2,1).am(2,1)]-G-am  →  D.(4,1)-G-am%G-A-am

    All (1,2) continuations — reversed chain, anchor at C-term with crosslink::

        K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-am
        →  K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1

    A flat sequential bracket cannot reverse direction without reusing an
    attachment slot. Such mixed continuations are preserved for validation,
    rather than reinterpreted as independent arms from the anchor.

    Single-monomer brackets stay inline. ``{}`` brackets are immune.
    """
    result = cabiln
    branches = []
    search_start = 0
    existing_tags = set(int(x) for x in _re.findall(r'!\s*(\d+)', cabiln))
    _xlink_ctr = [max(existing_tags, default=0) + 1]

    def _next_bracket(s, start):
        """Return (dot_pos, open_pos, close_pos+1) of the next .[...] or .{...}
        starting at or after `start`, using depth-aware scanning.  Returns None
        if not found."""
        i = start
        while i < len(s) - 1:
            if s[i] == '.' and s[i + 1] in ('[', '{'):
                open_ch = s[i + 1]
                close_ch = ']' if open_ch == '[' else '}'
                depth = 0
                j = i + 1
                while j < len(s):
                    if s[j] == open_ch:
                        depth += 1
                    elif s[j] == close_ch:
                        depth -= 1
                        if depth == 0:
                            return (i, i + 1, j + 1)
                    j += 1
            i += 1
        return None

    while True:
        span = _next_bracket(result, search_start)
        if span is None:
            break
        dot_pos, open_pos, close_pos = span
        m_start, m_end = dot_pos, close_pos
        open_ch = result[open_pos]

        if open_ch == '{':
            search_start = m_end
            continue

        bracket_content = result[open_pos + 1:close_pos - 1]

        # Sub-bracket arms present: check if they are all pure crosslink annotations
        # (e.g., .[TBMB(4,4)[.!2(5,4)][.!3(6,4)]]).  If so, convert to branch.
        # If arms have residue tokens, there is no clean % equivalent — skip.
        if '[' in bracket_content:
            flat_end = bracket_content.index('[')
            flat_part = bracket_content[:flat_end]
            hub_m = _re.fullmatch(r'([^(]+)\((\d+),(\d+)\)', flat_part.strip())
            if not hub_m:
                search_start = m_end
                continue
            # Extract sub-bracket arm contents
            sub_arms = []
            j = flat_end
            while j < len(bracket_content):
                if bracket_content[j] == '[':
                    depth = 0; k = j
                    while k < len(bracket_content):
                        if bracket_content[k] == '[': depth += 1
                        elif bracket_content[k] == ']':
                            depth -= 1
                            if depth == 0: break
                        k += 1
                    sub_arms.append(bracket_content[j + 1:k])
                    j = k + 1
                else:
                    j += 1
            # Each arm must be a single pure crosslink entry .!n(r,r)
            _xlink_pat = _re.compile(r'^\.?(!\d+)\((\d+),(\d+)\)$')
            arm_info = [_xlink_pat.match(a.strip()) for a in sub_arms]
            if (not sub_arms or any(m is None for m in arm_info)
                    or flat_part + ''.join(f'[{arm}]' for arm in sub_arms) != bracket_content):
                search_start = m_end
                continue
            # Build the conversion
            hub_name = hub_m.group(1)
            r_host_hub = hub_m.group(2)   # slot on main-chain host
            r_hub_host = hub_m.group(3)   # slot on hub connecting to host
            all_used = set(int(x) for x in _re.findall(r'!(\d+)', result))
            new_n = next(i for i in range(1, len(all_used) + 2) if i not in all_used)
            new_tag = f'!{new_n}'
            _xlink_ctr[0] = new_n + 1
            new_result = result[:m_start] + f'.{new_tag}({r_host_hub},{r_hub_host})' + result[m_end:]
            arm_tags = []
            for am in arm_info:
                tag_arm  = am.group(1)   # e.g., !2
                r_hub_arm = am.group(2)  # slot on hub side
                r_partner = am.group(3)  # slot on main-chain partner
                # Promote no-parens second endpoint in main chain to first occurrence
                new_result = _re.sub(
                    r'\.' + _re.escape(tag_arm) + r'(?![\w(])',
                    f'.{tag_arm}({r_partner},{r_hub_arm})',
                    new_result, count=1)
                arm_tags.append(tag_arm)
            result = new_result
            branch_str = hub_name + '.' + '.'.join([new_tag] + arm_tags)
            branches.append(branch_str)
            continue

        items = _parse_bracket_items(bracket_content)
        if (not items or any(rp is None or rt is None for _, rp, rt in items)
                or items[0][0].startswith('!')):
            search_start = m_end
            continue
        anchor_abbr, r_host, r_branch = items[0]

        if len(items) < 2:
            search_start = m_end
            continue

        # Flat multi-crosslink form: [TBMB(4,4).!2(5,4).!3(6,4)]
        # All non-anchor items are bare !n tags — same processing as sub-bracket form.
        _xpat = _re.compile(r'^!\d+$')
        if all(_xpat.match(abbr) for abbr, rp, rt in items[1:]):
            all_used = set(int(x) for x in _re.findall(r'!(\d+)', result))
            new_n = next(i for i in range(1, len(all_used) + 2) if i not in all_used)
            new_tag = f'!{new_n}'
            _xlink_ctr[0] = new_n + 1
            new_result = result[:m_start] + f'.{new_tag}({r_host},{r_branch})' + result[m_end:]
            arm_tags = []
            for abbr, rp, rt in items[1:]:
                # rp = hub-side slot, rt = partner-side slot
                new_result = _re.sub(
                    r'\.' + _re.escape(abbr) + r'(?![\w(])',
                    f'.{abbr}({rt},{rp})',
                    new_result, count=1)
                arm_tags.append(abbr)
            result = new_result
            branches.append(anchor_abbr + '.' + '.'.join([new_tag] + arm_tags))
            continue

        cont = [(rp, rt) for _, rp, rt in items[1:] if rp and rt]
        all_21 = cont and all(rp == '2' and rt == '1' for rp, rt in cont)
        # all_1x: all continuations enter via R1; out-slot may vary
        all_1x = cont and all(rp == '1' for rp, rt in cont)

        if all_21:
            tag = f'!{_xlink_ctr[0]}'
            _xlink_ctr[0] += 1
            host_marker = f'.{tag}({r_host},{r_branch})'
            branch_parts = [f'{anchor_abbr}.{tag}'] + [abbr for abbr, _, _ in items[1:]]
            branch_str = '-'.join(branch_parts)
            result = result[:m_start] + host_marker + result[m_end:]
            branches.append(branch_str)
        elif all_1x:
            tag = f'!{_xlink_ctr[0]}'
            _xlink_ctr[0] += 1
            host_marker = f'.{tag}({r_host},{r_branch})'
            result = result[:m_start] + host_marker + result[m_end:]
            all_12 = all(rt == '2' for _, rt in cont)
            if all_12:
                # Standard (1,2) continuations: reversed chain representation.
                reversed_cont = [abbr for abbr, _, _ in items[1:][::-1]]
                if r_branch == '2':
                    branch_str = '-'.join(reversed_cont + [anchor_abbr, tag])
                else:
                    branch_str = '-'.join(reversed_cont + [anchor_abbr]) + f'.{tag}'
                branches.append(branch_str)
            else:
                # Non-standard out-slot (e.g. E(1,4)): each monomer becomes a
                # standalone % segment, with one dedicated crosslink per bond.
                # K.[AEEA(4,2).E(1,4).C20FA(1,2)] →
                #   K.!1(4,2)-am%AEEA.!1.!2(1,4)%E.!2.!3(1,2)%C20FA.!3
                conts = items[1:]
                n = len(conts)
                cont_tags = [f'!{_xlink_ctr[0] + i}' for i in range(n)]
                _xlink_ctr[0] += n
                rp0, rt0 = conts[0][1], conts[0][2]
                branches.append(f'{anchor_abbr}.{tag}.{cont_tags[0]}({rp0},{rt0})')
                for k in range(n - 1):
                    rp_k, rt_k = conts[k + 1][1], conts[k + 1][2]
                    branches.append(
                        f'{conts[k][0]}.{cont_tags[k]}.{cont_tags[k + 1]}({rp_k},{rt_k})')
                branches.append(f'{conts[-1][0]}.{cont_tags[-1]}')
        else:
            search_start = m_end
            continue

    if branches:
        result += '%' + '%'.join(branches)
    return result


def cabiln_to_bracket(cabiln):
    """Convert CABILN branch (%) notation to bracket notation.

    Simple positional branches (``.(r,r)`` marker, first monomer is anchor)::

        K.(4,4)-G-am%gGlu-AEEA(2,1)-C20FA(2,1)
        →  K.[gGlu(4,4).AEEA(2,1).C20FA(2,1)]-G-am

    Crosslink branches (``.!n(r,r)`` marker, anchor tagged ``.!n`` in the
    branch — may be at the end or midpoint)::

        K.!1(4,4)-G-am%gGlu-AEEA(2,1)-C20FA(2,1).!1
        →  K.[C20FA(4,4).AEEA(1,2).gGlu(1,2)]-G-am

    Unannotated continuation monomers default to ``(2,1)`` (N→C).
    """
    if '%' not in cabiln and '\n' not in cabiln:
        return cabiln

    segments = _re.split(r'[%\n]', cabiln)
    segments = [s.strip() for s in segments if s.strip()]
    if len(segments) < 2:
        return cabiln

    main_seg = segments[0]
    branch_segs = segments[1:]
    deferred_branches = [
        branch for branch in branch_segs
        if _re.search(r'\.(?:[A-Za-z_]|\[|\{)|<', branch)
    ]
    branch_segs = [branch for branch in branch_segs if branch not in deferred_branches]

    def _normalize_terminal_marker(bs):
        """Convert standalone terminal !n tokens to inline .!n form.

        'G-G-!1'    -> 'G-G.!1'    (C-terminal marker on last monomer)
        '!1-G-G-am' -> 'G.!1-G-am' (N-terminal marker on first monomer)
        """
        parts = _re.split(r'(?<!\()[-](?!\))', bs)
        parts = [p.strip() for p in parts]
        if len(parts) >= 2 and _re.match(r'^!\d+$', parts[-1]):
            parts[-2] = parts[-2] + '.' + parts[-1]
            parts = parts[:-1]
        elif len(parts) >= 2 and _re.match(r'^!\d+$', parts[0]):
            parts[1] = parts[1] + '.' + parts[0]
            parts = parts[1:]
        return '-'.join(parts)

    branch_segs = [_normalize_terminal_marker(bs) for bs in branch_segs]

    # Separate branches: crosslink-connected vs positional
    positional = []
    crosslink = []
    for bs in branch_segs:
        if _re.search(r'\.!\d+', bs):
            crosslink.append(bs)
        else:
            positional.append(bs)

    # --- Phase 0: single-monomer chain segments connected by crosslinks ---
    # Pattern: main has .!n(r_host, r_branch); branch segments form a linear
    # chain MONO.!n.!m(a,b) → MONO.!m.!p(c,d) → MONO.!p (terminal).
    # These are emitted by cabiln_to_branch for non-standard (rt≠2) continuations.
    def _parse_chain_seg(bs):
        raw = [p.strip() for p in _re.split(r'(?<!\()[-](?!\))', bs) if p.strip()]
        if len(raw) != 1:
            return None
        seg = raw[0]
        monomer = _re.split(r'\.!', seg)[0].strip()
        outgoing = {}
        for m in _re.finditer(r'\.(!\d+)\((\d+),(\d+)\)', seg):
            outgoing[m.group(1)] = (m.group(2), m.group(3))
        all_tags = _re.findall(r'\.(!\d+)', seg)
        incoming = [t for t in all_tags if t not in outgoing]
        return monomer, incoming, outgoing

    seg_parse = {}
    for bs in crosslink:
        p = _parse_chain_seg(bs)
        if p is not None:
            seg_parse[bs] = p

    chain_processed = set()
    new_crosslink = []
    for bs in crosslink:
        if bs in chain_processed:
            continue
        if bs not in seg_parse:
            new_crosslink.append(bs)
            continue
        monomer, incoming, outgoing = seg_parse[bs]
        if not outgoing:
            new_crosslink.append(bs)
            continue
        anchor_tag = None
        for t in incoming:
            if _re.search(_re.escape(f'.{t}') + r'\((\d+),(\d+)\)', main_seg):
                anchor_tag = t
                break
        if anchor_tag is None:
            new_crosslink.append(bs)
            continue
        host_m = _re.search(_re.escape(f'.{anchor_tag}') + r'\((\d+),(\d+)\)', main_seg)
        r_host, r_branch = host_m.group(1), host_m.group(2)
        prev_processed = set(chain_processed)
        chain = [(monomer, r_host, r_branch)]
        chain_processed.add(bs)
        cur_out = outgoing
        ok = True
        while cur_out:
            if len(cur_out) != 1:
                ok = False
                break
            out_tag, (rp, rt) = next(iter(cur_out.items()))
            next_bs = None
            for bs2, (m2, inc2, out2) in seg_parse.items():
                if bs2 not in chain_processed and out_tag in inc2:
                    next_bs = bs2
                    break
            if next_bs is None:
                ok = False
                break
            m2, inc2, out2 = seg_parse[next_bs]
            chain.append((m2, rp, rt))
            chain_processed.add(next_bs)
            cur_out = out2
        if not ok:
            chain_processed = prev_processed
            new_crosslink.append(bs)
            continue
        bracket_items = [f'{chain[0][0]}({chain[0][1]},{chain[0][2]})']
        for mono, rp, rt in chain[1:]:
            bracket_items.append(f'{mono}({rp},{rt})')
        bracket_str = '.[' + '.'.join(bracket_items) + ']'
        main_seg = main_seg[:host_m.start()] + bracket_str + main_seg[host_m.end():]
    crosslink = new_crosslink

    # --- Phase 1: crosslink branches (.!n anchor in the branch) ---
    unconverted_crosslink = []
    for branch_seg in crosslink:
        all_tags = _re.findall(r'\.(!\d+)', branch_seg)
        unique_tags = list(dict.fromkeys(all_tags))

        if len(unique_tags) > 1:
            # Multi-tag hub (e.g., TBMB.!1.!2.!3) — single monomer only
            raw_parts = _re.split(r'(?<!\()[-](?!\))', branch_seg)
            raw_parts = [p.strip() for p in raw_parts if p.strip()]
            if len(raw_parts) != 1:
                unconverted_crosslink.append(branch_seg)
                continue
            hub_name = _re.split(r'\.!', branch_seg)[0].strip()
            tag_info = {}
            for tag in unique_tags:
                host_pat = _re.escape(f'.{tag}') + r'\((\d+),(\d+)\)'
                hm = _re.search(host_pat, main_seg)
                if hm:
                    tag_info[tag] = (hm.group(1), hm.group(2))
            if len(tag_info) < len(unique_tags):
                unconverted_crosslink.append(branch_seg)
                continue
            anchor_tag = unique_tags[0]
            r_host, r_branch = tag_info[anchor_tag]
            remaining = []
            for tag in unique_tags[1:]:
                th, tb = tag_info[tag]
                remaining.append(f'{tag}({tb},{th})')
            # Simplify remaining host tags BEFORE inserting bracket
            for tag in unique_tags[1:]:
                main_seg = _re.sub(
                    _re.escape(f'.{tag}') + r'\(\d+,\d+\)',
                    f'.{tag}', main_seg, count=1)
            inner = f'{hub_name}({r_host},{r_branch})'
            if not remaining:
                bracket_str = f'.[{inner}]'
            else:
                # Flat dot notation: all additional arms are pure crosslinks,
                # so sub-brackets are unnecessary — [hub.!2(5,4).!3(6,4)]
                flat_arms = ''.join(f'.{r}' for r in remaining)
                bracket_str = f'.[{inner}{flat_arms}]'
            host_pat = _re.escape(f'.{anchor_tag}') + r'\(\d+,\d+\)'
            main_seg = _re.sub(host_pat, bracket_str, main_seg, count=1)
            continue
        tag_m = _re.search(r'\.(!\d+)', branch_seg)
        if not tag_m:
            unconverted_crosslink.append(branch_seg)
            continue
        tag = tag_m.group(1)
        host_pat = _re.escape(f'.{tag}') + r'\((\d+),(\d+)\)'
        host_m = _re.search(host_pat, main_seg)
        if not host_m:
            unconverted_crosslink.append(branch_seg)
            continue
        r_host, r_branch = host_m.group(1), host_m.group(2)

        raw_parts = _re.split(r'(?<!\()[-](?!\))', branch_seg)
        raw_parts = [p.strip() for p in raw_parts if p.strip()]

        # Find which monomer is the anchor (has .!n)
        anchor_idx = None
        parsed = []
        for j, rp in enumerate(raw_parts):
            if f'.{tag}' in rp:
                anchor_idx = j
            clean = _re.sub(r'\.' + _re.escape(tag), '', rp)
            rg_m = _re.search(r'\((\d+),(\d+)\)', clean)
            abbr = _re.sub(r'\([^)]*\)', '', clean).strip()
            if rg_m:
                parsed.append((abbr, rg_m.group(1), rg_m.group(2)))
            else:
                parsed.append((abbr, None, None))
        if anchor_idx is None:
            unconverted_crosslink.append(branch_seg)
            continue

        # Build bracket: anchor first, then C-terminal side, then N-terminal
        bracket_items = [f'{parsed[anchor_idx][0]}({r_host},{r_branch})']

        # After anchor (C-terminal side in N→C branch) — keep R-groups
        for j in range(anchor_idx + 1, len(parsed)):
            abbr, rp, rt = parsed[j]
            bracket_items.append(
                f'{abbr}({rp},{rt})' if rp and rt else f'{abbr}(2,1)')

        if 0 < anchor_idx < len(parsed) - 1:
            after_arm = '.'.join(bracket_items[1:])
            bracket_items = [bracket_items[0] + f'[.{after_arm}]']

        # Before anchor (N-terminal side) — reverse order.
        # Branch items use swapped (rt,rp) perspective; swap back to get bracket (rp,rt).
        for j in range(anchor_idx - 1, -1, -1):
            abbr, cur_rp, cur_rt = parsed[j]
            if cur_rp and cur_rt:
                bracket_items.append(f'{abbr}({cur_rt},{cur_rp})')
            else:
                bracket_items.append(f'{abbr}(1,2)')

        bracket_str = '.[' + '.'.join(bracket_items) + ']'
        main_seg = (main_seg[:host_m.start()] + bracket_str
                    + main_seg[host_m.end():])

    # --- Phase 2: positional branches (.(r,r) marker) ---
    unconverted_positional = []
    for branch_seg in positional:
        host_m = _re.search(r'\.\((\d+),(\d+)\)', main_seg)
        if not host_m:
            unconverted_positional.append(branch_seg)
            continue
        r_host, r_branch = host_m.group(1), host_m.group(2)

        branch_parts = _re.split(r'(?<!\()[-](?!\))', branch_seg)
        branch_parts = [p.strip() for p in branch_parts if p.strip()]

        bracket_items = []
        for i, bp in enumerate(branch_parts):
            if i == 0:
                abbr = _re.sub(r'\([^)]*\)', '', bp).strip()
                bracket_items.append(f'{abbr}({r_host},{r_branch})')
            else:
                if '(' in bp:
                    bracket_items.append(bp)
                else:
                    bracket_items.append(f'{bp}(2,1)')
        bracket_str = '.[' + '.'.join(bracket_items) + ']'

        main_seg = (main_seg[:host_m.start()] + bracket_str
                    + main_seg[host_m.end():])

    remaining = unconverted_crosslink + unconverted_positional + deferred_branches
    if remaining:
        main_seg += '%' + '%'.join(remaining)
    return main_seg


def _check_bond_chemistry(mol1, at1, mol2, at2, bond_label='', warning_sink=None):
    """
    Validate the proposed inter-monomer bond.

    Standard peptide chemistry bonds pass silently:
      N–C(=O)  amide / isopeptide
      S–S      disulfide
      O–C(=O)  ester
      Se–Se    diselenide

    Bonds that are exotic but have documented peptide chemistry uses emit a
    UserWarning but are not blocked (thioester, sulfenamide, N–N hydrazide,
    non-carbonyl N–C, non-carbonyl O–C ether).

    Bonds with no plausible inter-monomer chemistry (e.g. C–C) raise
    ValueError — these almost always indicate wrong R-group numbers.
    """
    def _is_carbonyl_carbon(mol, atom):
        return (atom.GetAtomicNum() == 6
                and any(nb.GetAtomicNum() == 8
                        and mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                               .GetBondTypeAsDouble() == 2.0
                        for nb in atom.GetNeighbors()))

    a1 = mol1.GetAtomWithIdx(at1)
    a2 = mol2.GetAtomWithIdx(at2)
    sym1, sym2 = a1.GetAtomicNum(), a2.GetAtomicNum()
    pair = frozenset([sym1, sym2])

    # N(7)–C(6): amide/isopeptide — C must be carbonyl
    if pair == frozenset([7, 6]):
        c_mol, c_atom = (mol1, a1) if sym1 == 6 else (mol2, a2)
        if _is_carbonyl_carbon(c_mol, c_atom):
            return  # standard amide bond
        _emit_warning(
            f"Bond {bond_label}: N–C join where C is not a carbonyl carbon. "
            "This forms a C–N bond without amide character (e.g. reductive amination "
            "product). Intentional?",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # S(16)–S(16): disulfide
    if pair == frozenset([16]):
        return  # standard disulfide

    # O(8)–C(6): ester — C must be carbonyl
    if pair == frozenset([8, 6]):
        c_mol, c_atom = (mol1, a1) if sym1 == 6 else (mol2, a2)
        if _is_carbonyl_carbon(c_mol, c_atom):
            return  # standard ester
        _emit_warning(
            f"Bond {bond_label}: O–C join where C is not a carbonyl carbon. "
            "This forms an ether, not an ester. Intentional?",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # Se(34)–Se(34): diselenide
    if pair == frozenset([34]):
        return

    # N(7)–N(7): hydrazide / hydrazone
    if pair == frozenset([7]):
        _emit_warning(
            f"Bond {bond_label}: N–N join (hydrazide/hydrazone). "
            "Unusual in peptide chemistry — check R-group assignments.",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # S–C(=O): thioester; S–C(aliphatic): thioether
    if pair == frozenset([16, 6]):
        c_mol, c_atom = (mol1, a1) if sym1 == 6 else (mol2, a2)
        if _is_carbonyl_carbon(c_mol, c_atom):
            _emit_warning(
                f"Bond {bond_label}: S–C(=O) thioester bond. "
                "Valid for native chemical ligation but unusual for standard assembly. "
                "Intentional?",
                UserWarning, stacklevel=4, warning_sink=warning_sink,
            )
        else:
            _emit_warning(
                f"Bond {bond_label}: S–C(aliphatic) thioether bond. "
                "Valid for thioether-linked protecting groups (e.g. trt, acm) "
                "and related ligation chemistry. Intentional?",
                UserWarning, stacklevel=4, warning_sink=warning_sink,
            )
        return

    # S–N: sulfenamide
    if pair == frozenset([16, 7]):
        _emit_warning(
            f"Bond {bond_label}: S–N sulfenamide bond. "
            "Valid for Cys PTMs, oxidative macrolactamisation, bioconjugation, "
            "and NCL variant chemistry, but unusual in standard assembly. "
            "Intentional?",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # C(6)–C(6): RCM olefin staple (alkene–alkene metathesis) or other C–C bond
    if pair == frozenset([6]):
        def _is_vinyl(mol, atom):
            return any(
                mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble() == 2.0
                and nb.GetAtomicNum() == 6
                for nb in atom.GetNeighbors()
            )
        if _is_vinyl(mol1, a1) and _is_vinyl(mol2, a2):
            return  # all-hydrocarbon alkene staple (RCM)
        _emit_warning(
            f"Bond {bond_label}: C–C inter-monomer bond (non-vinyl). "
            "Unusual — check R-group assignments unless this is intentional "
            "bioconjugation chemistry.",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # O(8)–N(7): hydroxamic acid / hydroxylamine / isoxazole linkage
    if pair == frozenset([8, 7]):
        _emit_warning(
            f"Bond {bond_label}: O–N hydroxylamine/hydroxamic acid bond. "
            "Valid for O-amino acid caps (OBn_, OMe_) and hydroxamate "
            "bioconjugation, but unusual in standard peptide assembly. "
            "Intentional?",
            UserWarning, stacklevel=4, warning_sink=warning_sink,
        )
        return

    # Anything else — no plausible inter-monomer chemistry; raise so the caller
    # gets a clear message rather than a silent bad molecule or a raw RDKit crash.
    sym_names = {6: 'C', 7: 'N', 8: 'O', 16: 'S', 34: 'Se'}
    s1 = sym_names.get(sym1, str(sym1))
    s2 = sym_names.get(sym2, str(sym2))
    raise ValueError(
        f"Bond {bond_label}: {s1}–{s2} inter-monomer bond has no recognised "
        "peptide chemistry context. Check that the correct R-group numbers "
        "were specified for both monomers."
    )


class SequenceConstants:
    """
    A class to hold defaults values related to pyPept.Sequence objects.
    """
    def_path = "pyPept.data"
    def_lib_filename = "monomers.sdf"
    monomer_join = "-"
    chain_separator = "."
    csv_separator = ","
    helm_polymer = '|'
    max_rgroups = 4


# End of Sequence class-related constants definition.
############################################################

@dataclass
class ValidationReport:
    """Result of Sequence.validate(). Check .ok before calling .build()."""
    biln: str
    ok: bool
    bonds: list    # (m1_idx, slot1, m2_idx, slot2) tuples
    warnings: list # chemistry warning strings
    errors: list   # parse/validation error strings

    def __post_init__(self):
        self._seq = None  # set by Sequence.validate()

    def build(self):
        """Return the constructed Sequence, or raise ValueError if validation failed."""
        if not self.ok:
            msg = self.errors[0] if self.errors else "validation failed"
            raise ValueError(f"Cannot build invalid sequence: {msg}")
        return self._seq


class Sequence:
    """
    This class holds the information for a given sequence in BILN format,
    parsing monomer and bond information.
    """

    ############################################################
    def __init__(self, input_biln, path=SequenceConstants.def_path,
                 monomer_lib=SequenceConstants.def_lib_filename, fmt=None, *, warning_sink=None, track_source=False):
        """
        Instantiate a pyPept.Sequence object with the input BILN/CABILN sequence.

        :param input_biln: CABILN sequence string (default), or old BILN when fmt='biln'.
        :type input_biln: str
        :param path: an override option to specify a new location of a monomer
                     library.
        :type path: str
        :param monomer_lib: an override option to specific a different monomer
                           library file name.
        :type monomer_lib: str
        :param fmt: pass ``'biln'`` to accept old BILN crosslink notation
            (``Token(bondID,Rgroup)``), which is auto-converted to CABILN.
            Without this flag, old BILN notation raises ``ValueError``.
        :type fmt: str or None
        :param track_source: retain original token locations in ``s_sources``.
        """
        if not isinstance(input_biln, str) or not input_biln.strip():
            raise ValueError("CABILN must be a non-empty string.")
        if fmt not in (None, 'cabiln', 'biln'):
            raise ValueError(f"Unknown sequence format: {fmt!r}")

        # Variables to store the monomers and bonds
        self.s_inputbiln = input_biln
        self._warning_sink = warning_sink
        self.s_mid = -1
        self.s_bonds = []
        self.s_nbonds = 0
        self._used_slots = set()
        self.s_monomers = []
        self.s_nmonomers = 0
        self.__is_valid = True

        if track_source and fmt == 'biln':
            raise ValueError('Source tracking requires CABILN; convert old BILN first.')

        # Old BILN bare-integer crosslink notation requires explicit opt-in via fmt='biln'.
        if _OLD_BILN_RE.search(input_biln):
            m = _OLD_BILN_RE.search(input_biln)
            tok, bid, rg = m.group(1), m.group(2), m.group(3)
            if fmt != 'biln':
                raise ValueError(
                    f"Old BILN crosslink notation detected: {tok}({bid},{rg}). "
                    f"Pass fmt='biln' to auto-convert, or call biln_to_cabiln() explicitly."
                )
            input_biln = biln_to_cabiln(input_biln)

        # Protect synthetic SMILES before interpreting CABILN separators.
        expanded = SourceText.original(input_biln) if track_source else input_biln
        self.s_sources = []

        # Pre-expand synthetic tokens.  N-caps (`<smi>_`) and C-caps (`_<smi>`)
        # are matched BEFORE bare alpha-AA `<smi>` to ensure the trailing or
        # leading underscore is consumed with the bracket.  Each token is
        # replaced by a deterministic synthetic monomer symbol; the SMILES is
        # remembered so the synthetic monomer rows are appended to monomer_df
        # below.
        self._synthetic_aa_smiles = {}   # symbol -> sidechain SMILES (alpha-AA)
        self._synthetic_caps = {}        # symbol -> (cap_smi, side)
        def _ncap_repl(_m):
            cs = _m.group(1)
            sym = _synthetic_cap_symbol(cs, 'n') + '_'  # trailing _ matches Ac_ / Bz_ naming
            self._synthetic_caps[sym] = (cs, 'n')
            return _source_synthetic(_m, sym)
        def _ccap_repl(_m):
            cs = _m.group(1)
            sym = '_' + _synthetic_cap_symbol(cs, 'c')  # leading _ matches _NH2 / _OBn naming
            self._synthetic_caps[sym] = (cs, 'c')
            return _source_synthetic(_m, sym)
        def _syn_repl(_m):
            sc = _m.group(1)
            sym = _synthetic_aa_symbol(sc)
            self._synthetic_aa_smiles[sym] = sc
            return _source_synthetic(_m, sym)
        expanded = _source_sub(_SYN_NCAP_RE, _ncap_repl, expanded)
        expanded = _source_sub(_SYN_CCAP_RE, _ccap_repl, expanded)
        expanded = _source_sub(_SYN_AA_RE, _syn_repl, expanded)
        expanded, _branch_rgroup = _expand_inline_caps(expanded, warning_sink=warning_sink)

        seq = split_outside(expanded,
                            by_element=SequenceConstants.monomer_join,
                            outside='[]')
        if any(not token.strip() for token in seq):
            raise ValueError("Empty monomer between backbone separators.")
        seq = _source_join(SequenceConstants.monomer_join, seq)
        # Remove dangling '-' around chain breaks
        self.s_biln = _source_sub(r'[-]*\.[-]*',
                             SequenceConstants.chain_separator, seq)

        # Read the monomer dictionary
        default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
        if (path == SequenceConstants.def_path
                and monomer_lib == SequenceConstants.def_lib_filename):
            from pyPept.monomer_store import library_path
            monomer_df_filepath = library_path()
        else:
            monomer_df_filepath = files(path).joinpath(monomer_lib)

        if monomer_df_filepath.is_file() is False:
            monomer_df_filepath = default_monomer_df_filepath

        self.monomer_df = get_monomer_info(str(monomer_df_filepath))

        # Register synthetic alpha-AA monomers built from `<sidechain>` tokens.
        # Malformed SMILES must never silently change the requested structure.
        for _sym, _sc in self._synthetic_aa_smiles.items():
            if _sym in self.monomer_df.index:
                continue
            _romol = _build_synthetic_aa(_sc)
            _name = f'synthetic_aa<{_sc}>'
            if _romol is None:
                raise ValueError(f"Synthetic AA token '<{_sc}>' contains invalid SMILES.")
            # Scan the synthetic mol for ALL R-group dummies (not just R1/R2).
            # Build m_Rgroups and m_chem_types dynamically: R1 = backbone_n,
            # R2 = backbone_c, R3+ = inferred from the atom the dummy bonds to.
            _r_slots = []
            for _a in _romol.GetAtoms():
                if _a.GetAtomicNum() == 0 and _a.GetIsotope() >= 1:
                    _r_slots.append(_a.GetIsotope())
            _r_slots = sorted(set(_r_slots))
            _max_slot = max(_r_slots) if _r_slots else 2
            _rgroups_list = [None] * _max_slot
            _cts_parts = []
            for _slot in _r_slots:
                # Determine chem_type & leaving group
                _ct, _lg = _infer_synth_chem_type(_romol, _slot)
                if _ct is None:
                    continue
                _rgroups_list[_slot - 1] = _lg
                _cts_parts.append(f'{_slot}:{_ct}')
            # Fill any missing slots with None
            _rgroups_list = [r if r is not None else '' for r in _rgroups_list]
            _cts = ','.join(_cts_parts)
            _romol.SetProp('m_chem_types', _cts)
            self.monomer_df.loc[_sym] = {
                'm_Rgroups':    _rgroups_list,
                'm_abbr':       _sym,
                'm_name':       _name,
                'm_type':       'aa',
                'm_subtype':    'synthetic',
                'm_chem_types': _cts,
                'ID':           _sym,
                'm_romol':      _romol,
            }

        # Register synthetic N-cap / C-cap monomers.
        # N-cap: R2 only (carbonyl C bonds to next residue's R1 via backbone_amide)
        # C-cap: R1 only (N bonds to previous residue's R2)
        for _sym, (_cs, _side) in self._synthetic_caps.items():
            if _sym in self.monomer_df.index:
                continue
            _romol = _build_synthetic_cap(_cs, _side)
            if _romol is None:
                raise ValueError(
                    f"Synthetic {_side}-cap token contains invalid SMILES: {_cs!r}")
            if _side == 'n':
                _rgroups = ['', '[OH]']
                _cts = '2:backbone_c'
            else:
                _rgroups = ['[H]', '']
                _cts = '1:backbone_n'
            _romol.SetProp('m_chem_types', _cts)
            self.monomer_df.loc[_sym] = {
                'm_Rgroups':    _rgroups,
                'm_abbr':       _sym,
                'm_name':       f'synthetic_{_side}cap<{_cs}>',
                'm_type':       'cap',
                'm_subtype':    'synthetic',
                'm_chem_types': _cts,
                'ID':           _sym,
                'm_romol':      _romol,
            }

        self.__parse_biln()
        self.__parse_biln_bonds()

    ########################################################################################
    def __parse_biln(self):
        """Split BILN string into individual monomers, do some basic checks,
        then construct the monomers and store the chemical representations
        of the monomers in the Sequence object.
        """

        biln = self.s_biln
        chains = split_outside(biln,
                               by_element=SequenceConstants.chain_separator,
                               outside='[]')

        self.s_chains = {}
        self.s_chains["s_nChains"] = len(chains)
        self.s_chains["s_cType"] = [None] * len(chains)
        self.s_chains["s_monomerIDs"] = [None] * len(chains)

        m_idx = 0
        # Iterate over the chains
        for num_chain, chain in enumerate(chains):
            residues = split_outside(chain,
                                     by_element=SequenceConstants.monomer_join,
                                     outside='[]')
            monomer_types = set()
            monomer_ids = []
            for _res_pos, res in enumerate(residues):
                # Extract the residue information
                # Strip bond annotations: (n,m) and (!n,m) crosslink IDs
                resname = _source_sub(r'\((?:!\w+|\d+),\d+\)', '', res)
                if isinstance(resname, SourceText):
                    self.s_sources.append(resname.tracker.occurrence(resname))
                resname = re.sub(r'^[\[{](.*?)[\]}]$', '\\1', resname)

                # Check if name exists in MonomerDic — with degeneracy resolution

                if resname not in self.monomer_df.index:
                    if resname in self.monomer_df.attrs.get('_ambiguous_aliases', {}):
                        raise ValueError(
                            f"Ambiguous monomer alias '{resname}' names multiple library entries; "
                            "use a canonical symbol instead.")
                    # 1) Check stored abbreviations and unambiguous CSV aliases.
                    synonyms = self.monomer_df.attrs.get('_synonyms', {})
                    if resname in synonyms:
                        resname = synonyms[resname]
                    # 2) Check degenerate cap aliases (R1/R2 structural pairs)
                    elif resname in self.monomer_df.attrs.get('_degen_aliases', {}):
                        degen = self.monomer_df.attrs['_degen_aliases']
                        variants = degen[resname]
                        res_idx = _res_pos
                        if res_idx == 0:
                            chosen = [v for v in variants
                                      if v.endswith('_') and not v.startswith('_')]
                            if not chosen:
                                chosen = [v for v in variants if v.endswith('_')]
                        elif res_idx == len(residues) - 1:
                            chosen = [v for v in variants
                                      if v.startswith('_') and not v.endswith('_')]
                            if not chosen:
                                chosen = [v for v in variants if v.startswith('_')]
                        else:
                            raise ValueError(
                                f"Degenerate cap '{resname}' used mid-chain "
                                f"(position {res_idx}); cannot resolve variant.")
                        if chosen:
                            resname = chosen[0]
                        else:
                            raise ValueError(
                                f"Monomer {resname!r} has no terminal variant for this position.")
                    else:
                        raise ValueError(
                            f"Monomer {resname!r} is not in the monomer library.")

                # Add additional information of the monomers
                mm_info = self.monomer_df.loc[resname, :].to_dict()
                mm_value = {'m_name': resname,
                            'm_name_in_biln': res,
                            'm_chainID': num_chain}
                mm_value.update(mm_info)

                # Append the monomer in the sequence object
                self.__append_monomer(mm_value)

                monomer_ids.append(m_idx)
                monomer_types.add(mm_info['m_type'])
                m_idx += 1

            # Save the type of monomer
            if len(monomer_types) == 1:
                if monomer_types == {'aa'}:
                    self.s_chains["s_cType"][num_chain] = 'peptide'
                else:
                    self.s_chains["s_cType"][num_chain] = 'chem'
            else:
                self.s_chains["s_cType"][num_chain] = 'mixed'

            self.s_chains["s_monomerIDs"][num_chain] = monomer_ids

            # end of looping over each residue in a chain
        # end of looping over all chains in the peptide

    # end of internal-use method for parsing BILN.

    ########################################################################################
    def __append_monomer(self, monomer):
        """Place a new monomer at the end of the sequence, i.e.
        append the current monomer to the monomer list.

        :param monomer: monomer dictionary of added monomer
        :type monomer: dictionary
        """
        self.s_mid += 1

        # Fields required for the addition of the monomer
        keys_needed = ['m_name', 'm_abbr', 'm_name_in_biln', 'm_type',
                       'm_subtype', 'm_chainID', 'm_Rgroups', 'm_romol',
                       'm_chem_types']

        # check if all necessary keys are available in the monomer definition
        for key in keys_needed:
            try:
                monomer[key]
            except KeyError:
                raise ValueError(f'Key {key} missing in monomer description')

        # Construct the monomer dictionary
        m_dict = {}
        m_dict.update({key: monomer[key] for key in keys_needed})

        # Update class variables
        self.s_nmonomers += 1
        self.s_monomers.append(m_dict)

    ########################################################################################
    def __parse_biln_bonds(self):
        """
        Function to parse bond information from the BILN string
        """

        # Local variable to store bond information
        bond_info = []
        bond_info_helm = []

        residues = split_outside(self.s_biln,
                                 SequenceConstants.chain_separator +
                                 SequenceConstants.monomer_join, '[]')
        nres = len(residues)
        # Iterate over residues
        for num_res, res in enumerate(residues):

            if num_res < nres - 1:
                if self.__get_monomer_prop('m_chainID', num_res) == \
                        self.__get_monomer_prop('m_chainID', num_res + 1):
                    # We have a bond between the two monomers.
                    # Read attachment points directly from the mol via isotope labels
                    # (CHUCKLES convention) rather than the sidecar index property.
                    mol1 = self.__get_monomer_prop('m_romol', num_res)
                    mol2 = self.__get_monomer_prop('m_romol', num_res + 1)
                    r1_attach = _attachment_idx(mol1, 2)  # slot 2 = R2 = C-term exit
                    r2_attach = _attachment_idx(mol2, 1)  # slot 1 = R1 = N-term entry
                    if r1_attach is None or r2_attach is None:
                        raise ValueError(
                            f"Cannot form backbone R2→R1 bond between residues "
                            f"{num_res} and {num_res + 1}: one or both monomers "
                            f"lack the required attachment point. Check that "
                            f"terminal-only monomers (e.g. ac, am, fmoc) are "
                            f"not placed mid-chain.")
                    else:
                        _check_bond_chemistry(
                            mol1, r1_attach, mol2, r2_attach,
                            bond_label=f'backbone {num_res}→{num_res+1}',
                            warning_sink=self._warning_sink)
                        self.__add_bond(
                            num_res, r1_attach, num_res + 1, r2_attach, slot1=2, slot2=1)

            # Find connections with other chains.
            # Bond IDs may be plain integers ("1") or !n markers ("!1").
            match = re.findall(r"\((!\w+|\d+),(\d+)\)", res)
            for i, m_value in enumerate(match):
                b_idx, rgrp_idx = m_value[0], int(m_value[1])
                mol = self.__get_monomer_prop('m_romol', num_res)
                attach = _attachment_idx(mol, rgrp_idx)
                if attach is None:
                    raise ValueError(
                        f"Residue {num_res} has no R{rgrp_idx} attachment "
                        f"point but BILN specifies a bond using R{rgrp_idx}.")
                bond_info.append([b_idx, num_res, rgrp_idx, attach])
                bond_info_helm.append([b_idx, num_res, rgrp_idx])

            # For multiple branching (two crosslinks on the same residue)
            match = re.findall(r"\((!\w+|\d+),(\d+);(!\w+|\d+),(\d+)\)", res)
            for i, m_value in enumerate(match):
                b1_idx, rgrp1_idx = m_value[0], int(m_value[1])
                b2_idx, rgrp2_idx = m_value[2], int(m_value[3])
                mol = self.__get_monomer_prop('m_romol', num_res)
                attach1 = _attachment_idx(mol, rgrp1_idx)
                attach2 = _attachment_idx(mol, rgrp2_idx)
                if attach1 is None:
                    raise ValueError(
                        f"Residue {num_res} has no R{rgrp1_idx} attachment "
                        f"point but BILN specifies a bond using R{rgrp1_idx}.")
                if attach2 is None:
                    raise ValueError(
                        f"Residue {num_res} has no R{rgrp2_idx} attachment "
                        f"point but BILN specifies a bond using R{rgrp2_idx}.")
                bond_info.append([b1_idx, num_res, rgrp1_idx, int(attach1)])
                bond_info.append([b2_idx, num_res, rgrp2_idx, int(attach2)])
                bond_info_helm.append([b1_idx, num_res, rgrp1_idx])
                bond_info_helm.append([b2_idx, num_res, rgrp2_idx])

        # Now that we have the bond info we collect additional information
        nbonds = int(len(bond_info))

        if (nbonds % 2) == 1:
            raise ValueError("Must be an even number of extra bonds.")

        if nbonds > 0:
            # Group bond_info entries by bond identifier (str — may be "1" or "!1").
            grouped: dict = {}
            for entry in bond_info:
                grouped.setdefault(entry[0], []).append(entry)

            bond = []
            for bid, entries in grouped.items():
                if len(entries) != 2:
                    raise ValueError(
                        f"Bond identifier {bid!r} has {len(entries)} endpoint(s); "
                        "expected exactly 2 (one on each partner residue).")
                bond.append(entries)

            # Collect the data needed to add to bond information in Sequence
            for bondx in bond:

                if bondx[0][1] > bondx[1][1]:
                    bondx[0], bondx[1] = bondx[1], bondx[0]

                m1, at1, m2, at2 = bondx[0][1], bondx[0][3], bondx[1][1], bondx[1][3]
                slot1 = bondx[0][2]  # 1-based slot for m1's attachment
                slot2 = bondx[1][2]  # 1-based slot for m2's attachment
                mol1 = self.__get_monomer_prop('m_romol', m1)
                mol2 = self.__get_monomer_prop('m_romol', m2)
                _check_bond_chemistry(
                    mol1, at1, mol2, at2,
                    bond_label=f'crosslink bond-id {bondx[0][0]!r} '
                               f'(R{bondx[0][2]} of residue {m1} <-> '
                               f'R{bondx[1][2]} of residue {m2})',
                    warning_sink=self._warning_sink)
                self.__add_bond(m1, at1, m2, at2, slot1=slot1, slot2=slot2)

        # Filter unique bonds
        self.__only_unique_bonds()

    ## end of Sequence.__parse_biln_bonds()

    ########################################################################################
    def __add_bond(self, m_id1, atom1, m_id2, atom2, slot1=None, slot2=None):
        """
        Add an entry in the bond list

        :param m_id1: index of monomer 1
        :param atom1: index of bond atom in monomer 1
        :param m_id2: index of monomer 2
        :param atom2: index of bond atom in monomer 2
        :param slot1: 1-based R-group slot for atom1 (optional; avoids ambiguity
                      when an atom neighbours multiple dummies, e.g. backbone N
                      which has both [1*] and [3*]).
        :param slot2: 1-based R-group slot for atom2 (optional; same rationale).
        """
        slot1 = slot1 if slot1 is not None else _slot_for_attachment(
            self.s_monomers[m_id1]['m_romol'], atom1)
        slot2 = slot2 if slot2 is not None else _slot_for_attachment(
            self.s_monomers[m_id2]['m_romol'], atom2)
        endpoints = ((m_id1, slot1), (m_id2, slot2))
        if endpoints[0] == endpoints[1]:
            raise ValueError("A bond cannot join the same attachment slot to itself.")
        for monomer_id, slot in endpoints:
            if (monomer_id, slot) in self._used_slots:
                raise ValueError(
                    f"Residue {monomer_id} R{slot} is already used by another bond.")
        self._used_slots.update(endpoints)
        self.s_nbonds += 1
        entry = [m_id1, atom1, m_id2, atom2]
        if slot1 is not None:
            entry.extend([slot1, slot2])
        self.s_bonds.append(entry)

    ########################################################################################
    def get_monomer(self, idx):
        """
        Return the dictionary of ith monomer by index

        :param idx: index of monomer
        :type idx: int
        """

        mon = None
        try:
            mon = self.s_monomers[idx]
        except IndexError:
            warnings.warn(f'Cannot access monomer index {idx} in sequence')

        return mon

    ########################################################################################
    def __get_monomer_prop(self, property_val, idx):
        """
        Return a property of the i-th monomer

        :param property_val: which property to fetch
        :type property_val: str
        :param idx: index of monomer
        :type idx: int
        :return: desired property
        """

        m_prop = None
        try:
            m_dict = self.get_monomer(idx)
        except ValueError:
            m_prop = None

        try:
            m_prop = m_dict[property_val]
        except ValueError:
            warnings.warn(
                f'Cannot access property {property_val} in monomer {idx}')

        return m_prop

    ############################################################################
    def length(self):
        """
        Return the length of the sequence
        """
        return len(self.s_monomers)

    ############################################################################
    def __len__(self):
        """
        Return the length of the sequence.
        """
        return self.length()

    ############################################################################
    def __only_unique_bonds(self):
        """
        In __parseBILNBonds it might happen that a bond is added twice
        (once from the amide section, once if it is)
        in addition set explicitly by (1,2) sets.
        Here, we filter these bonds out

        :return: filtered bond list with unique bonds only
        """
        seen = set()
        unique = []
        for b in self.s_bonds:
            key = tuple(b)
            if key not in seen:
                seen.add(key)
                unique.append(b)
        self.s_bonds = unique

    ############################################################################
    def is_valid(self):
        """Flag if the initialized pyPept.Sequence is valid.
        This includes check of R-groups present in the monomer dictionary,
        and R-group connectivity when explicitly given.

        :return: bool
        """
        return self.__is_valid

    ############################################################################
    @classmethod
    def validate(cls, biln, path=SequenceConstants.def_path,
                 monomer_lib=SequenceConstants.def_lib_filename, fmt=None):
        """Dry-run parse and chemistry-check a BILN/CABILN string.

        Constructs the Sequence, captures all warnings and errors, and returns
        a ValidationReport without raising.  Call report.build() to get the
        Sequence object if report.ok is True.

        :param biln: BILN or CABILN sequence string
        :param path: monomer library package path (default: built-in)
        :param monomer_lib: monomer SDF filename (default: monomers.sdf)
        :param fmt: pass ``'biln'`` to accept old BILN crosslink notation.
        :returns: ValidationReport
        """
        caught_warns = []
        errors = []
        seq = None

        try:
            seq = cls(biln, path=path, monomer_lib=monomer_lib, fmt=fmt,
                      warning_sink=caught_warns.append)
        except ValueError as exc:
            errors.append(str(exc))

        bonds = []
        if seq is not None:
            for entry in seq.s_bonds:
                m1, m2 = entry[0], entry[2]
                slot1 = entry[4] if len(entry) > 4 else None
                slot2 = entry[5] if len(entry) > 5 else None
                bonds.append((m1, slot1, m2, slot2))

        ok = seq is not None and not errors
        report = ValidationReport(
            biln=biln, ok=ok, bonds=bonds,
            warnings=caught_warns, errors=errors,
        )
        report._seq = seq
        return report

    ############################################################################
    def __bool__(self):
        """
        Returns if a pyPept.Sequence object has been constructed with a valid
        BILN.

        :seealso: is_valid

        >>> s = Sequence('P-E-P')
        >>> bool(s)
        True
        >>> s = Sequence('P-nonsense-P')
        >>> bool(s)
        False
        """
        return self.is_valid()

    ## End of the Sequence class declaration.


##########################################################################
# Additional functions
##########################################################################

def get_monomer_info(path):
    """Return a detached monomer table from the shared, versioned library."""
    from pyPept.monomer_store import monomer_table

    return monomer_table(path)


########################################################################
def greekify(mol, name_aa):
    """
    Converts the atom names into relative names using the Greek alphabet.

    :param mol: RDKit molecule object
    :type mol: RDKit molecule
    :param name_aa: name of the AA that will be modified
    :type name_aa: str
    """

    # Some local fixes to name correctly the natural amino acids if required
    fix_aa = {'W': {'CE1': 'CE2', 'NE2': 'NE1', 'CZ1': 'CZ3', 'CH': 'CH2',
                    'CD1': 'CD2', 'CD2': 'CD1'},
              'N': {'OD2': 'OD1', 'ND1': 'ND2'},
              'I': {'CD': 'CD1', 'CG1': 'CG2', 'CG2': 'CG1'},
              'T': {'OG2': 'OG1', 'CG1': 'CG2'},
              'P': {'CG2': 'CG', 'CG1': 'CD'}}

    greek = list('ABGDEZHTIKLMNXOPRS')
    greekdex = defaultdict(list)
    ca_atom = get_atom_by_name(mol, 'CA')

    # Recognize if the atom is not part of the backbone
    for atom in mol.GetAtoms():
        is_backbone = (atom.GetPDBResidueInfo() is not None and
                       atom.GetPDBResidueInfo().GetName().strip() in (
                           'LOWER', 'UPPER', 'N', 'CA', 'C', 'H', 'HA', 'O',
                           'OXT'))
        if atom.GetSymbol() != 'H' and atom.GetSymbol()[
            0] != 'R' and not is_backbone:
            n_atom = len(
                Chem.GetShortestPath(mol, ca_atom.GetIdx(), atom.GetIdx())) - 1
            greekdex[n_atom].append(atom)

    # Iterate over the list of atoms and the greek letters
    for k in greekdex:
        if len(greek) <= k:
            pass
        elif len(greekdex[k]) == 0:
            pass
        elif len(greekdex[k]) == 1:
            # Special case
            new_name = f'{greekdex[k][0].GetSymbol()}{greek[k]}'
            if name_aa in fix_aa:
                if new_name in fix_aa[name_aa]:
                    digit = fix_aa[name_aa][new_name][-1]
                    name = f'{greekdex[k][0].GetSymbol(): >2}{greek[k]}{digit}'
                else:
                    name = f'{greekdex[k][0].GetSymbol(): >2}{greek[k]} '
            else:
                name = f'{greekdex[k][0].GetSymbol(): >2}{greek[k]} '

            greekdex[k][0].GetPDBResidueInfo().SetName(name)

        elif len(greekdex[k]) < 36:
            list_atom = list(string.digits + string.ascii_uppercase)[1:]
            for i, atom in enumerate(greekdex[k]):

                new_name = f'{atom.GetSymbol()}{greek[k]}{list_atom[i]}'
                if name_aa in fix_aa:
                    if new_name in fix_aa[name_aa]:
                        if name_aa != 'P':
                            digit = fix_aa[name_aa][new_name][-1]
                            name = f'{atom.GetSymbol(): >2}{greek[k]}{digit}'
                        else:
                            digit = fix_aa[name_aa][new_name][-1]
                            name = f'{atom.GetSymbol(): >2}{digit} '
                    else:
                        name = f'{atom.GetSymbol(): >2}{greek[k]}{list_atom[i]}'
                else:
                    name = f'{atom.GetSymbol(): >2}{greek[k]}{list_atom[i]}'

                greekdex[k][i].GetPDBResidueInfo().SetName(name)

        else:
            pass


############################################################

def split_outside(string, by_element, outside, keep_marker=True):
    """
    Splits a string by delimiter only if outside of a given delimiter

    :param string: string to be split
    :param by_element: delimiter(s) by which to be split
    :param outside: only split if outside of this
    :param keep_marker: if True keep the chunk marker, remove otherwise

    :return splitChains: split string as list
    """
    # by can be more than 1 character
    by_element = list(by_element)

    # if outside is only one character (e.g. ', "), double it for start and end
    if len(outside) == 1:
        outside = outside + outside

    # Special character
    grpsep = chr(29)

    openers = {outside[0]}
    closers = {outside[1]}
    if outside == '[]':
        openers.add('{')
        closers.add('}')

    out = []
    inside = False
    for char_position in range(len(string)):
        i = string[char_position:char_position+1]
        if i in openers:
            if inside:
                if keep_marker:
                    j = i
                else:
                    j = ''
                inside = False
            else:
                inside = True
                if keep_marker:
                    j = i
                else:
                    j = ''
        elif i in closers:
            inside = False
            if keep_marker:
                j = i
            else:
                j = ''
        else:
            if not inside and i in by_element:
                j = grpsep
            else:
                j = i
        out.append(j)

    # Do the final split
    split_chains = _source_join('', out).split(grpsep)
    return split_chains


############################################################

def get_atom_by_name(mol, name):
    """
    Get an atom RDKit object based on a particular name

    :param mol: RDKit object with the molecule
    :type mol: RDKit molecule
    :param name: name of the atom
    :type name: str

    :return: the atom RDKit object
    """

    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info and info.GetName().strip() == name.strip():
            selected_atom = atom

    return selected_atom


########################################################################################

def get_monomer_codes(df_name):
    """
    Get the PDB codes from the monomer dictionary
    :param df_name: monomer dictionary
    :type df_name: dataframe

    :return: dictionary with the monomer PDB codes
    """
    # Get the values from the monomer dictionary
    monomers = {}
    for idx in df_name.index:
        symbol = df_name.at[idx, 'm_abbr']
        pdb_name = df_name.at[idx, 'pdbName']
        monomers[symbol] = pdb_name

    return monomers


############################################################

def correct_pdb_atoms(seq, path=SequenceConstants.def_path,
                      monomer_lib=SequenceConstants.def_lib_filename):
    """
    Pipeline to assign the correct atom names to the pyPept object

    :param seq: pyPept Sequence object
    :type seq: pyPept.sequence
    :return: the modified pyPept Sequence with correct atom names
    """

    # Special case two main N- and C- terminal caps
    names_cap = {'ac': ['CH3', 'C', 'O'], 'am': ['N']}

    # Read the monomer dataframe
    default_monomer_df_filepath = files(SequenceConstants.def_path).joinpath(SequenceConstants.def_lib_filename)
    if (path == SequenceConstants.def_path
            and monomer_lib == SequenceConstants.def_lib_filename):
        from pyPept.monomer_store import library_path
        monomer_df_filepath = library_path()
    else:
        monomer_df_filepath = files(path).joinpath(monomer_lib)

    if monomer_df_filepath.is_file() is False:
        monomer_df_filepath = default_monomer_df_filepath

    new_df = get_monomer_info(str(monomer_df_filepath))
    if 'pdbName' not in new_df.columns:
        raise ValueError(
            f"Monomer library {monomer_df_filepath} lacks the 'pdbName' metadata "
            "required for PDB atom naming. Supply a library with PDB residue "
            "names, or use --noconf for 2D output."
        )

    # Get monomer codes
    monomers = get_monomer_codes(new_df)

    # Iterate over the monomers
    mm_list = seq.s_monomers
    for i, monomer in enumerate(mm_list):
        mol = monomer['m_romol']
        name = monomer['m_abbr']
        backbone_smiles = Chem.MolFromSmiles('NCC(=O)')
        if i == len(mm_list) - 1:
            backbone_smiles = Chem.MolFromSmiles('NCC(=O)O')
        backbone_atoms = mol.GetSubstructMatches(backbone_smiles)
        aa_flag = 0
        if backbone_atoms:
            bb_index = list(backbone_atoms[0])
            type_mon = new_df.loc[new_df['m_abbr'] == name, 'm_type'].item()
            if type_mon == 'aa':
                aa_flag = 1

        # Iterate over the atoms
        counter = 0
        counter_non = 0
        for j, atom in enumerate(mol.GetAtoms()):
            if aa_flag == 1:
                pos = -1
                if atom.GetIdx() in bb_index:
                    pos = bb_index.index(atom.GetIdx())
                if pos == 0:
                    atomname = ' N  '
                elif pos == 1:
                    atomname = ' CA '
                elif pos == 2:
                    atomname = ' C  '
                elif pos == 3:
                    atomname = ' O  '
                elif pos == 4:
                    atomname = ' OXT'
                else:
                    counter += 1
                    atomname = f' {atom.GetSymbol()}{counter} '
            else:
                # Special case for main capping groups
                if name in ('ac', 'am'):
                    if atom.GetSymbol()[0] != 'R':
                        if names_cap[name][j] == 'CH3':
                            atomname = f' {names_cap[name][j]}'
                        else:
                            atomname = f' {names_cap[name][j]}  '
                else:
                    counter_non += 1
                    if counter_non < 10:
                        atomname = f' {atom.GetSymbol()}{counter_non} '
                    else:
                        atomname = f' {atom.GetSymbol()}{counter_non}'

            # Assign the atom object to the peptide molecule
            info = atom.GetPDBResidueInfo()
            if info is None:
                atom.SetMonomerInfo(Chem.AtomPDBResidueInfo(atomName=atomname,
                                                            serialNumber=atom.GetIdx(),
                                                            residueName=f'{monomers[name]}',
                                                            residueNumber=i + 1,
                                                            chainId="A"))

        # Rename atoms using the greek nomenclature
        if aa_flag == 1:
            greekify(mol, name)

    return seq

    # end of definition of correct_pdb_atoms()

############################################################
# End of sequence.py
############################################################