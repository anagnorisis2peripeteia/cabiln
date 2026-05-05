#!/usr/bin/env python3
"""Probe candidate test cases for bracket<->branch round-trips."""
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent.parent / 'src'))
from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch

def probe(desc, form, val, expected_out=None):
    try:
        if form == 'bracket':
            out = cabiln_to_branch(val)
            back = cabiln_to_bracket(out)
        else:
            out = cabiln_to_bracket(val)
            back = cabiln_to_branch(out)
        rt = back == val
        mark = 'OK' if rt else 'FAIL'
        print(f'{mark} [{form}] {desc}')
        if not rt:
            print(f'    IN:   {val!r}')
            print(f'    OUT:  {out!r}')
            print(f'    BACK: {back!r}')
        elif expected_out and out != expected_out:
            print(f'    OUT:  {out!r}')
            print(f'    EXP:  {expected_out!r}')
        if expected_out and out != expected_out:
            print(f'    ** expected_out mismatch **')
    except Exception as e:
        print(f'ERROR [{form}] {desc}: {e}')

print('=== Single-attachment ===')
probe('single monomer R4->R2', 'bracket', 'ac-K.[A(4,2)]-G-am')  # passthrough expected
probe('single monomer R4->R1', 'bracket', 'ac-K.[A(4,1)]-G-am')  # passthrough expected

print('\n=== N-terminal branch ===')
probe('N-term branch bracket->branch', 'branch', 'ac-K.!1(4,2)-G-am%!1-G-G-am')
probe('N-term branch branch->bracket', 'bracket', 'ac-K.[?]')  # skip

print('\n=== C-terminal !1 marker (xfail territory) ===')
probe('C-term !1 marker branch->bracket', 'branch', 'ac-K.!1(4,2)-G-am%G-G-!1')
probe('inline .!1 branch->bracket', 'branch', 'ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1')

print('\n=== Two lipid brackets (both (1,2)) ===')
probe('two lipid brackets bracket->branch',
      'bracket',
      'ac-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-am')

print('\n=== PEG/AEEA repeat ===')
probe('two AEEA bracket->branch', 'bracket', 'ac-K.[AEEA(4,2).AEEA(1,2)]-G-am')
probe('three AEEA bracket->branch', 'bracket', 'ac-K.[AEEA(4,2).AEEA(1,2).AEEA(1,2)]-G-am')

print('\n=== Three positional branches ===')
probe('three positional branches branch->bracket', 'branch',
      'ac-K.(4,1)-A-K.(4,1)-A-K.(4,1)-am%G-am%A-am%G-G-am')

print('\n=== Semaglutide-like (simplified 8AA with lipid) ===')
probe('sema-like bracket->branch', 'bracket',
      'His-Aib-Glu-Gly-Thr-Phe-Thr-K.[gGlu(4,4).AEEA(1,2).C18FA(1,2)]-am')
probe('sema-like bracket->branch (no C18FA, use C20FA)', 'bracket',
      'His-Aib-Glu-Gly-Thr-Phe-Thr-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-am')

print('\n=== Retatrutide bracket->branch->bracket (bracket is canonical) ===')
ret_bracket = ('Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-'
               'K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-'
               'Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am')
probe('retatrutide bracket->branch', 'bracket', ret_bracket)

print('\n=== Mixed: crosslink + positional branch ===')
probe('cyclic + positional branch (bracket)', 'bracket',
      '!1-A-K.[G(4,1).am(2,1)]-G-A-!1')
probe('cyclic + positional branch (branch)', 'branch',
      '!1-A-K.(4,1)-G-A-!1%G-am')

print('\n=== No cap variants ===')
probe('no ac/am caps branch->bracket', 'branch', 'K.(4,1)-G-A%G-A-am')
probe('no ac cap bracket->branch', 'bracket', 'K.[G(4,1).A(2,1).am(2,1)]-G-A')

print('\n=== Longer positional branches (3-4 residues) ===')
probe('4-residue branch bracket->branch', 'bracket',
      'ac-K.[G(4,1).A(2,1).G(2,1).am(2,1)]-G-am')
probe('4-residue branch branch->bracket', 'branch',
      'ac-K.(4,1)-G-am%G-A-G-am')

print('\n=== R3 attachment (e.g., isopeptide to sidechain acid) ===')
probe('R3 branch bracket->branch', 'bracket', 'ac-D.[G(4,3)]-A-am')
probe('R3 single-monomer stays as passthrough?', 'bracket', 'ac-E.[G(4,3)]-A-am')

print('\n=== Two crosslink branches (both !n inline) ===')
probe('two crosslink branches branch->bracket', 'branch',
      'ac-K.!1(4,4)-G-K.!2(4,4)-am%C20FA-AEEA-gGlu.!1%C20FA-AEEA-gGlu.!2')
