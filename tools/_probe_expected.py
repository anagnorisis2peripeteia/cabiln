#!/usr/bin/env python3
"""Print exact intermediate forms for new test cases."""
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent.parent / 'src'))
from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch

cases = [
    # (desc, form, val)
    # Two lipid brackets
    ('two_lipid_brackets',   'bracket', 'ac-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-am'),
    # PEG/AEEA
    ('two_AEEA',             'bracket', 'ac-K.[AEEA(4,2).AEEA(1,2)]-G-am'),
    ('three_AEEA',           'bracket', 'ac-K.[AEEA(4,2).AEEA(1,2).AEEA(1,2)]-G-am'),
    # Three positional branches
    ('three_branches',       'branch',  'ac-K.(4,1)-A-K.(4,1)-A-K.(4,1)-am%G-am%A-am%G-G-am'),
    # Semaglutide-like (short)
    ('sema_like',            'bracket', 'His-Aib-Glu-Gly-Thr-Phe-Thr-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-am'),
    # Retatrutide bracket canonical
    ('retatrutide_bracket',  'bracket', 'Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am'),
    # Cyclic + plain positional branch
    ('cyclic_plus_plain',    'bracket', '!1-A-K.[G(4,1).am(2,1)]-G-A-!1'),
    ('cyclic_plus_plain_b',  'branch',  '!1-A-K.(4,1)-G-A-!1%G-am'),
    # No caps
    ('no_caps_branch',       'branch',  'K.(4,1)-G-A%G-A-am'),
    ('no_caps_bracket',      'bracket', 'K.[G(4,1).A(2,1).am(2,1)]-G-A'),
    # 4-residue branch
    ('four_res_bracket',     'bracket', 'ac-K.[G(4,1).A(2,1).G(2,1).am(2,1)]-G-am'),
    ('four_res_branch',      'branch',  'ac-K.(4,1)-G-am%G-A-G-am'),
    # R3 attachment
    ('r3_branch',            'bracket', 'ac-D.[G(4,3)]-A-am'),
    # Two crosslink (1,2) branches
    ('two_crosslink_lipid',  'branch',  'ac-K.!1(4,4)-G-K.!2(4,4)-am%C20FA-AEEA-gGlu.!1%C20FA-AEEA-gGlu.!2'),
]

for desc, form, val in cases:
    try:
        if form == 'bracket':
            out = cabiln_to_branch(val)
            back = cabiln_to_bracket(out)
        else:
            out = cabiln_to_bracket(val)
            back = cabiln_to_branch(out)
        rt = back == val
        status = 'OK' if rt else 'FAIL'
        print(f'--- {desc} [{status}] ---')
        print(f'  INPUT:  {val!r}')
        print(f'  OUT:    {out!r}')
        if not rt:
            print(f'  BACK:   {back!r}')
    except Exception as e:
        print(f'--- {desc} [ERROR] ---')
        print(f'  {e}')
