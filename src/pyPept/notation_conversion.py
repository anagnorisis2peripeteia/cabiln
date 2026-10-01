"""Library-free compatibility conversions between historical notation layouts.

These text adapters preserve unsupported text. Chemical resolution and canonical
formatting belong to Sequence and Peptide; callers needing those guarantees use
pyPept.inputs instead. Public imports from pyPept.sequence remain supported.
"""

import re

from pyPept.notation import (
    INLINE_BOND_RE as _INLINE_BOND_RE, bracket_chain, parse_chain_entry,
    supports_bracket_token, text_crosslink_declarations,
    validate_crosslink_declarations, legacy_attachment_slot,
)

_OLD_BILN_RE = re.compile(r'(?<![.!\w\[{])([A-Za-z]\w*)\((\d+),(\d+)\)')
_BRACKET_COLOURS = ['\033[93m', '\033[96m', '\033[92m', '\033[95m', '\033[91m']
_ANSI_RESET = '\033[0m'


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
        tok1, tok2 = m1.group(1), m2.group(1)
        # This text adapter historically preserves noncanonical numeric spelling.
        rg1, rg2 = (
            str(legacy_attachment_slot(3)) if value == "3" else value
            for value in (m1.group(3), m2.group(3))
        )
        replacements[m1.start()] = (m1, f'{tok1}.!{bid}({rg1},{rg2})')
        replacements[m2.start()] = (m2, f'{tok2}.!{bid}')

    result = biln
    for pos in sorted(replacements, reverse=True):
        m, replacement = replacements[pos]
        result = result[:m.start()] + replacement + result[m.end():]
    return result


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
    for ch in biln:
        if ch in ('[', '{'):
            colour = _BRACKET_COLOURS[bracket_count % len(_BRACKET_COLOURS)]
            result.append(colour + ch)
        elif ch in (']', '}'):
            result.append(ch + _ANSI_RESET)
            bracket_count += 1
        else:
            result.append(ch)
    return ''.join(result)


def _parse_bracket_items(bracket_content):
    """Parse bracket content into list of (abbr, r_prev, r_this) tuples."""
    parts = bracket_content.split('.')
    items = []
    for p in parts:
        pm = re.fullmatch(r'([^(]+)\((\d+),(\d+)\)', p)
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
    if validate_crosslink_declarations(
            text_crosslink_declarations(cabiln), strict=False) is None:
        return cabiln
    result = cabiln
    branches = []
    search_start = 0
    existing_tags = set(int(x) for x in re.findall(r'!\s*(\d+)', cabiln))
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

    def detach_hub(start, end, anchor, arms):
        """Move either supported hub spelling into one crosslinked segment."""
        nonlocal result
        name, host_slot, own_slot = anchor
        used = {int(value) for value in re.findall(r'!(\d+)', result)}
        number = next(i for i in range(1, len(used) + 2) if i not in used)
        tag = f'!{number}'
        _xlink_ctr[0] = number + 1
        result = result[:start] + f'.{tag}({host_slot},{own_slot})' + result[end:]
        for arm, hub_slot, partner_slot in arms:
            result = re.sub(
                r'\.' + re.escape(arm) + r'(?![\w(])',
                f'.{arm}({partner_slot},{hub_slot})', result, count=1)
        branches.append(name + '.' + '.'.join([tag] + [arm for arm, _, _ in arms]))

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
            hub_m = re.fullmatch(r'([^(]+)\((\d+),(\d+)\)', flat_part.strip())
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
            _xlink_pat = re.compile(r'^\.?(!\d+)\((\d+),(\d+)\)$')
            arm_info = [_xlink_pat.match(a.strip()) for a in sub_arms]
            if (not sub_arms or any(m is None for m in arm_info)
                    or flat_part + ''.join(f'[{arm}]' for arm in sub_arms) != bracket_content):
                search_start = m_end
                continue
            detach_hub(m_start, m_end, hub_m.groups(), [arm.groups() for arm in arm_info])
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
        _xpat = re.compile(r'^!\d+$')
        if all(_xpat.match(abbr) for abbr, rp, rt in items[1:]):
            detach_hub(m_start, m_end, items[0], items[1:])
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

    if validate_crosslink_declarations(
            text_crosslink_declarations(cabiln), strict=False) is None:
        return cabiln

    segments = re.split(r'[%\n]', cabiln)
    segments = [s.strip() for s in segments if s.strip()]
    if len(segments) < 2:
        return cabiln

    main_seg = segments[0]
    branch_segs = segments[1:]
    deferred_branches = [
        branch for branch in branch_segs
        if re.search(r'\.(?:[A-Za-z_]|\[|\{)|<', branch)
    ]
    branch_segs = [branch for branch in branch_segs if branch not in deferred_branches]

    def _normalize_terminal_marker(bs):
        """Convert standalone terminal !n tokens to inline .!n form.

        'G-G-!1'    -> 'G-G.!1'    (C-terminal marker on last monomer)
        '!1-G-G-am' -> 'G.!1-G-am' (N-terminal marker on first monomer)
        """
        parts = re.split(r'(?<!\()[-](?!\))', bs)
        parts = [p.strip() for p in parts]
        if len(parts) >= 2 and re.match(r'^!\d+$', parts[-1]):
            parts[-2] = parts[-2] + '.' + parts[-1]
            parts = parts[:-1]
        elif len(parts) >= 2 and re.match(r'^!\d+$', parts[0]):
            parts[1] = parts[1] + '.' + parts[0]
            parts = parts[1:]
        return '-'.join(parts)

    branch_segs = [_normalize_terminal_marker(bs) for bs in branch_segs]

    # Only rewrite branches whose full token/annotation syntax is understood.
    # Preserve unsupported text byte-for-byte; never strip parentheses to make
    # a malformed symbol or endpoint look valid.
    parsed_branches = {}
    endpoint_counts = {}
    for bs in branch_segs:
        entries = tuple(parse_chain_entry(part.strip()) for part in bs.split('-'))
        if any(
            entry is None or not supports_bracket_token(entry.token)
            for entry in entries
        ):
            return cabiln
        for entry in entries:
            for tag, slots in entry.markers:
                # This legacy adapter folds numeric labels. Named labels remain
                # available to the modern parser/writer without prefix matching.
                if not re.fullmatch(r'!\d+', tag):
                    return cabiln
                endpoint_counts[tag] = endpoint_counts.get(tag, 0) + 1
        parsed_branches[bs] = entries
    for marker in _INLINE_BOND_RE.finditer(main_seg):
        tag = marker.group(1)
        endpoint_counts[tag] = endpoint_counts.get(tag, 0) + 1
    if any(count > 2 for count in endpoint_counts.values()):
        return cabiln

    # Separate branches: crosslink-connected vs positional.
    positional = []
    crosslink = []
    for bs in branch_segs:
        if any(entry.markers for entry in parsed_branches[bs]):
            crosslink.append(bs)
        else:
            positional.append(bs)

    # --- Phase 0: single-monomer chain segments connected by crosslinks ---
    # Pattern: main has .!n(r_host, r_branch); branch segments form a linear
    # chain MONO.!n.!m(a,b) → MONO.!m.!p(c,d) → MONO.!p (terminal).
    # These are emitted by cabiln_to_branch for non-standard (rt≠2) continuations.
    def _parse_chain_seg(bs):
        entries = parsed_branches[bs]
        if len(entries) != 1:
            return None
        entry = entries[0]
        if entry.slots is not None:
            return None
        outgoing = {tag: slots for tag, slots in entry.markers if slots is not None}
        incoming = [tag for tag, slots in entry.markers if slots is None]
        if len(incoming) != 1 or len(outgoing) > 1:
            return None
        return entry.token, incoming, outgoing

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
            if re.search(re.escape(f'.{t}') + r'\((\d+),(\d+)\)', main_seg):
                anchor_tag = t
                break
        if anchor_tag is None:
            new_crosslink.append(bs)
            continue
        host_m = re.search(re.escape(f'.{anchor_tag}') + r'\((\d+),(\d+)\)', main_seg)
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
        entries = parsed_branches[branch_seg]
        all_tags = [tag for entry in entries for tag, _ in entry.markers]
        unique_tags = list(dict.fromkeys(all_tags))

        if len(unique_tags) > 1:
            # Multi-tag hub (e.g., TBMB.!1.!2.!3) — single monomer only
            if len(entries) != 1 or entries[0].slots is not None:
                unconverted_crosslink.append(branch_seg)
                continue
            hub_name = entries[0].token
            tag_info = {}
            for tag in unique_tags:
                host_pat = re.escape(f'.{tag}') + r'\((\d+),(\d+)\)'
                hm = re.search(host_pat, main_seg)
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
                main_seg = re.sub(
                    re.escape(f'.{tag}') + r'\(\d+,\d+\)',
                    f'.{tag}', main_seg, count=1)
            inner = f'{hub_name}({r_host},{r_branch})'
            if not remaining:
                bracket_str = f'.[{inner}]'
            else:
                # Flat dot notation: all additional arms are pure crosslinks,
                # so sub-brackets are unnecessary — [hub.!2(5,4).!3(6,4)]
                flat_arms = ''.join(f'.{r}' for r in remaining)
                bracket_str = f'.[{inner}{flat_arms}]'
            host_pat = re.escape(f'.{anchor_tag}') + r'\(\d+,\d+\)'
            main_seg = re.sub(host_pat, bracket_str, main_seg, count=1)
            continue
        tag = unique_tags[0]
        host_pat = re.escape(f'.{tag}') + r'\((\d+),(\d+)\)'
        host_m = re.search(host_pat, main_seg)
        if not host_m:
            unconverted_crosslink.append(branch_seg)
            continue
        r_host, r_branch = host_m.group(1), host_m.group(2)

        # Find the complete marker on its monomer, never a numeric prefix of
        # another tag. Continuation slots belong to the monomer, not the marker.
        anchor_idx = None
        parsed = []
        for j, entry in enumerate(entries):
            if any(marker == tag for marker, _ in entry.markers):
                anchor_idx = j
            previous, own = entry.slots or (None, None)
            parsed.append((entry.token, previous, own))
        if anchor_idx is None:
            unconverted_crosslink.append(branch_seg)
            continue

        bracket_str = bracket_chain(parsed, anchor_idx, r_host, r_branch).text
        main_seg = (main_seg[:host_m.start()] + bracket_str
                    + main_seg[host_m.end():])

    # --- Phase 2: positional branches (.(r,r) marker) ---
    unconverted_positional = []
    for branch_seg in positional:
        host_m = re.search(r'\.\((\d+),(\d+)\)', main_seg)
        if not host_m:
            unconverted_positional.append(branch_seg)
            continue
        r_host, r_branch = host_m.group(1), host_m.group(2)

        bracket_items = []
        for i, entry in enumerate(parsed_branches[branch_seg]):
            previous, own = (r_host, r_branch) if i == 0 else entry.slots or ('2', '1')
            bracket_items.append(f'{entry.token}({previous},{own})')
        bracket_str = '.[' + '.'.join(bracket_items) + ']'

        main_seg = (main_seg[:host_m.start()] + bracket_str
                    + main_seg[host_m.end():])

    remaining = unconverted_crosslink + unconverted_positional + deferred_branches
    if remaining:
        main_seg += '%' + '%'.join(remaining)
    return main_seg
