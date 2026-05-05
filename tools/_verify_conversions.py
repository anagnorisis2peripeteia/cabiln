#!/usr/bin/env python3
"""Verify bracket<->branch conversions for proposed new test cases."""
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).parent.parent / 'src'))

from pyPept.sequence import cabiln_to_bracket, cabiln_to_branch

cases = [
    # (description, form, value)
    # Cyclic + lipid bracket round-trip
    ("cyclic+lipid bracket->branch",
     "bracket",
     "!1-A-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-A-!1"),

    ("cyclic+lipid branch->bracket",
     "branch",
     "!1-A-K.!2(4,4)-G-A-!1%C20FA-AEEA-gGlu.!2"),

    # Two brackets (one plain, one lipid) on same peptide
    ("two brackets: plain+lipid  bracket->branch",
     "bracket",
     "ac-K.[A(4,1)]-G-K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-am"),

    # Retatrutide bracket->branch (K.!1(4,2) means K-R4 to AEEA-R2 or similar)
    ("retatrutide branch->bracket",
     "branch",
     "Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-K.!1(4,2)-A-Q-Aib-A-F-I-E-Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am%C20FA-gGlu-AEEA-!1"),

    # Two-branch plain round-trip: verify intermediate bracket form
    ("two branches branch->bracket",
     "branch",
     "ac-C.(4,4)-A-K.(4,1)-G-am%C-A-am%G-am"),

    # Two brackets: verify intermediate branch form
    ("two brackets bracket->branch",
     "bracket",
     "ac-C.[C(4,4).A(2,1).am(2,1)]-A-K.[G(4,1).am(2,1)]-G-am"),

    # Long chain with two positional branches
    ("long chain two branches branch->bracket",
     "branch",
     "ac-A-G-K.(4,1)-A-A-E.(4,1)-G-am%G-A-am%A-G-am"),
]

for desc, form, val in cases:
    print(f"\n{'='*60}")
    print(f"  {desc}")
    print(f"  INPUT:  {val!r}")
    try:
        if form == "bracket":
            result = cabiln_to_branch(val)
            print(f"  BRANCH: {result!r}")
            back = cabiln_to_bracket(result)
            print(f"  BACK:   {back!r}")
            print(f"  RT OK:  {back == val}")
        else:
            result = cabiln_to_bracket(val)
            print(f"  BRACKET:{result!r}")
            back = cabiln_to_branch(result)
            print(f"  BACK:   {back!r}")
            print(f"  RT OK:  {back == val}")
    except Exception as e:
        print(f"  ERROR:  {e}")
