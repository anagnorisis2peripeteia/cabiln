"""Lower CABILN scopes and connections while retaining source occurrence records.

This pass has no monomer-library dependency. It uses the shared notation grammar
and produces the explicit chains consumed by Sequence, with source tracking kept
inside SourceText. It does not infer attachment chemistry or resolve definitions.
"""

import re

from pyPept.source import (
    SourceText, bond_marker as _source_bond_marker, group as _source_group,
    join as _source_join, record as _source_record, sub as _source_sub,
    scope_group as _source_scope_group,
)
from pyPept.notation import (
    BracketArm, bracket_regions, parse_bracket_group, MAX_NOTATION_CHARACTERS,
    INLINE_BOND_RE as _INLINE_BOND_RE, crosslink_declarations,
    validate_crosslink_declarations, split_outside,
)

_INLINE_CAP_RE = re.compile(r'\.([A-Za-z_]\w*)\((\d+),(\d+)\)')


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


def _expand_inline_caps(biln):
    """Pre-process CABILN notation before chain/residue splitting.

    Supported forms:
    - '%' or newline as segment separator (main chain first, branches follow)
    - .boc(3,1)            named cap — auto-assigns bond ID >= 100, appended as pendant
    - .!1(4,2)             crosslink first endpoint — explicit R-groups, defines inverse
    - .!1                  crosslink second endpoint — parens omitted, inverse inferred
    - !1-A-B-C             branch N-terminal marker — attaches via R1 of first residue
    - A-B-C-!1             branch C-terminal marker — attaches via R2 of last residue
    - .[A(r,s).B(t,u)...]  sequential reaction bracket — left-to-right chain of named
                            cap attachments.  For each step after the first, the first
                            R-group number refers to the PRECEDING fragment (not the host
                            residue).  The host receives only the first bond annotation.
                            Auto bond IDs are assigned in left-to-right order.
    - .[A(r,s).[B(t,u)]]   child scope attaches at A and returns to A on close;
                            protected .{...} and legacy [.arm] scopes also nest.

    Returns (expanded_biln, branch_rgroup).  expanded_biln uses '.' as chain
    separator; branch_rgroup maps each !x to its partner rgroup from first occurrence.
    """
    if len(biln) > MAX_NOTATION_CHARACTERS:
        raise ValueError(
            f'Notation exceeds {MAX_NOTATION_CHARACTERS:,} characters; '
            'split it into smaller documents.')
    biln = _preprocess_cabiln(biln)
    original_segments = [s for s in biln.split('%') if s.strip()]
    segments = original_segments

    # Parse each bracket once, then collect declarations in source order from
    # those entries and the surrounding inline text. All endpoint spellings use
    # the same inverse-pair check before terminal/bare markers infer their slots.
    parsed_brackets = []
    declarations = []

    def inline_declarations(text):
        return (
            (m.group(1), m.group(2), m.group(3))
            for m in _INLINE_BOND_RE.finditer(text) if m.group(2) is not None
        )

    for seg in segments:
        groups = []
        previous_end = 0
        for start, end in bracket_regions(seg):
            declarations.extend(inline_declarations(seg[previous_end:start]))
            scope = parse_bracket_group(seg[start:end])
            declarations.extend(crosslink_declarations(scope.entries))
            root_start = seg.rfind('-', 0, start) + 1
            host = seg[root_start:start].split('.', 1)[0]
            if not host:
                raise ValueError('Sequential bracket has no host monomer.')
            groups.append((scope, host))
            previous_end = end
        declarations.extend(inline_declarations(seg[previous_end:]))
        parsed_brackets.append(groups)
    _seen = validate_crosslink_declarations(declarations)

    appended = []
    _ctr = [100]
    _count = {}

    def _marker_slot(tok, hr, cr):
        if _count.get(tok, 0) >= 2:
            raise ValueError(
                f"Bond {tok!r}: 3rd endpoint found (max 2 allowed).")
        if hr is None:
            if tok not in _seen:
                raise ValueError(
                    f"Bond marker {tok!r} used without parens but no prior "
                    f".{tok}(y,z) first-occurrence found.")
            fhr, fcr = _seen[tok]
            hr, cr = fcr, fhr
        _count[tok] = _count.get(tok, 0) + 1
        return str(int(hr))

    def _inline_owner(match):
        root_start = match.string.rfind('-', 0, match.start()) + 1
        host = match.string[root_start:match.start()].split('.', 1)[0]
        return _source_sub(r'\((?:!\w+|\d+),\d+\)', '', host)

    def _sub_bond(m):
        tok = m.group(1)
        hr = _marker_slot(tok, m.group(2), m.group(3))
        _source_bond_marker(_source_group(m), int(hr))
        return f'({tok},{hr})'

    def _sub_cap(m):
        tok, hr, cr = _source_group(m, 1), m.group(2), m.group(3)
        owner = _inline_owner(m)
        _source_record(tok, _source_group(m), 'inline', host=owner)
        _source_scope_group(
            _source_group(m), host=owner, kind='inline', opening='', closing='',
        )
        bid = _ctr[0]; _ctr[0] += 1
        appended.append(tok + f'({bid},{cr})')
        return f'({bid},{hr})'

    def _sub_bracket(scope, host):
        bracket = scope.text
        frag_tokens, frag_bond_parts = [], []
        host_bonds = []
        protected = scope.opening == '{'
        _source_scope_group(
            bracket, host=host, opening=scope.opening, closing=scope.closing,
            protected=protected,
        )
        # Each frame retains its own current monomer. Closing a child simply
        # pops that frame, so following children resume at the parent's cursor.
        stack = [{
            'scope': scope, 'position': 0, 'current': -1, 'parent': None,
            'protected': protected, 'first': None,
        }]
        while stack:
            frame = stack[-1]
            group = frame['scope']
            position = frame['position']
            if position == len(group.entries):
                stack.pop()
                if group.returns_anchor and stack:
                    anchor = frame['first']
                    if anchor is not None:
                        stack[-1]['current'] = anchor
                        if stack[-1]['first'] is None:
                            stack[-1]['first'] = anchor
                continue
            step = group.entries[position]
            frame['position'] += 1
            current = frame['current']
            owner = host if current == -1 else frag_tokens[current]
            bonds = host_bonds if current == -1 else frag_bond_parts[current]
            if isinstance(step, BracketArm):
                protected = frame['protected'] or step.opening == '{'
                _source_scope_group(
                    step.text, host=owner, parent=group.text, kind='arm',
                    opening=step.opening, closing=step.closing,
                    protected=protected, legacy=step.legacy,
                )
                stack.append({
                    'scope': step, 'position': 0, 'current': current,
                    'parent': group.text, 'protected': protected, 'first': None,
                })
                continue
            tok, previous, own = step.token, step.previous_slot, step.own_slot
            arm = group.text if group is not scope else None
            if tok.startswith('!'):
                slot = _marker_slot(tok, previous, own)
                _source_bond_marker(
                    step.text, int(slot), owner=owner, bracket=bracket, arm=arm,
                    scope=group.text,
                )
                bonds.append(f'({tok},{slot})')
                continue
            _source_record(
                tok, step.text, 'bracket', bracket, arm,
                terminal=position == len(group.entries) - 1,
                host=owner, scope=group.text, parent_scope=frame['parent'],
                protected=frame['protected'], legacy_arm=group.legacy,
            )
            bid = _ctr[0]
            _ctr[0] += 1
            bonds.append(f'({bid},{previous})')
            frag_tokens.append(tok)
            frag_bond_parts.append([f'({bid},{own})'])
            frame['current'] = len(frag_tokens) - 1
            if frame['first'] is None:
                frame['first'] = frame['current']

        for tok, bonds in zip(frag_tokens, frag_bond_parts):
            appended.append(tok + ''.join(bonds))
        return ''.join(host_bonds)

    processed_segments = []
    for segment_number, seg in enumerate(segments):
        if isinstance(seg, SourceText):
            seg.tracker.segment = segment_number
            for root in split_outside(original_segments[segment_number], '-', '[]'):
                seg.tracker.root(root, segment_number)
        seg = _handle_terminal_bond_markers(seg, _seen, _count)
        # Terminal markers only alter text outside brackets. Locate the unchanged
        # groups again after that edit, reusing their already parsed entries.
        parts = []
        last_end = 0
        for (start, end), (scope, host) in zip(
                bracket_regions(seg), parsed_brackets[segment_number]):
            parts.append(seg[last_end:start])
            parts.append(_sub_bracket(scope, host))
            last_end = end
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
