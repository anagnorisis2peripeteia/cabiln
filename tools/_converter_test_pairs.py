#!/usr/bin/env python3
"""
11 hand-crafted bracket <-> branch conversion pairs.

Rules:
  - Single .mod stays inline (no branch for 1 monomer)
  - {} brackets are immune to conversion
  - All (2,1) continuations: anchor first (N-term), positional .(r,r)
  - All (1,2) continuations: reverse chain, anchor at C-term, crosslink .!n
  - .!n marks the crosslink point wherever it sits in the branch
  - No (r,r) on continuation monomers in branch — dash is implicit R2->R1
  - Multi-tag hubs (TBMB): branch->bracket only (one-way)
"""

# (description, bracket_form, expected_branch_form[, "oneway"])
# If bracket == branch, no conversion expected.
# "oneway" = branch->bracket only (bracket->branch stays unchanged).

PAIRS = [
    # --- Should NOT convert ---

    # 1. Single mod — only 1 monomer in bracket
    (
        "single_mod_no_convert",
        "ac-K.[G(4,1)]-G-am",
        "ac-K.[G(4,1)]-G-am",
    ),

    # 2. Protected {} bracket — immune
    (
        "protected_bracket_no_convert",
        "ac-K.{G(4,1).A(2,1).am(2,1)}-G-am",
        "ac-K.{G(4,1).A(2,1).am(2,1)}-G-am",
    ),

    # --- (2,1) continuations: anchor at N-term, positional ---

    # 3. Two monomers, host R4->anchor R1
    (
        "two_mono_21_positional",
        "ac-K.[G(4,1).A(2,1)]-G-am",
        "ac-K.(4,1)-G-am%G-A",
    ),

    # 4. Three monomers, standard peptide branch
    (
        "three_mono_21_positional",
        "ac-K.[G(4,1).A(2,1).am(2,1)]-G-am",
        "ac-K.(4,1)-G-am%G-A-am",
    ),

    # 5. Four monomers off D side chain
    (
        "four_mono_21_positional",
        "A-D.[G(4,1).A(2,1).L(2,1).am(2,1)]-A-am",
        "A-D.(4,1)-A-am%G-A-L-am",
    ),

    # 6. R4->R4 crosslink anchor with (2,1) continuations
    #    gGlu joins host via R4, chain extends right via R2->R1
    #    Anchor at N-term, positional .(4,4)
    (
        "R4R4_anchor_21_positional",
        "ac-K.[gGlu(4,4).A(2,1).am(2,1)]-G-am",
        "ac-K.(4,4)-G-am%gGlu-A-am",
    ),

    # 7. Two hosts both with convertible (2,1) brackets
    (
        "two_hosts_both_convert",
        "ac-K.[G(4,1).A(2,1)]-K.[G(4,1).L(2,1)]-am",
        "ac-K.(4,1)-K.(4,1)-am%G-A%G-L",
    ),

    # --- (1,2) continuations: reverse chain, anchor at C-term, crosslink ---

    # 8. Lipid linker — all (1,2), reversed chain
    #    gGlu(4,4).AEEA(1,2).C20FA(1,2) reversed = C20FA-AEEA-gGlu
    #    gGlu.!1 at C-term, crosslink back to host
    (
        "lipid_linker_12_reversed",
        "K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-am",
        "K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1",
    ),

    # --- Branch -> Bracket (crosslink cases) ---

    # 9. Crosslink midpoint anchor — gGlu in middle of branch
    #    A-gGlu.!1-L: A is N-term side, L is C-term side
    #    Bracket: anchor gGlu first, then C-term (2,1), then N-term reversed (1,2)
    (
        "crosslink_midpoint_anchor",
        "K.[gGlu(4,4).L(2,1).A(1,2)]-G-am",
        "K.!1(4,4)-G-am%A-gGlu.!1-L",
    ),

    # 10. Crosslink anchor at C-term of branch
    #     C20FA-AEEA-gGlu.!1: anchor gGlu at end
    #     Bracket: gGlu first, AEEA and C20FA reversed (1,2)
    (
        "crosslink_anchor_Cterm",
        "K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-am",
        "K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1",
    ),

    # --- Multi-tag hub (one-way: branch -> bracket only) ---

    # 11. TBMB 3-way hub — single monomer with 3 crosslinks
    #     Branch: TBMB.!1.!2.!3 — three tags, one per C
    #     Bracket: first tag becomes anchor, rest become .!n(r,r) entries
    #     cabiln_to_branch does NOT convert this back (crosslink entries can't be dash-chained)
    (
        "tbmb_multitag_hub_oneway",
        "ac-C.[TBMB(4,4).!2(5,4).!3(6,4)]-A-A-C.!2-A-A-C.!3-am",
        "ac-C.!1(4,4)-A-A-C.!2(4,5)-A-A-C.!3(4,6)-am%TBMB.!1.!2.!3",
        "oneway",
    ),
]


if __name__ == "__main__":
    import sys
    sys.path.insert(0, "src")
    import warnings
    warnings.filterwarnings("ignore")
    from pyPept.sequence import cabiln_to_branch, cabiln_to_bracket

    all_pass = True

    print("=" * 70)
    print("BRACKET -> BRANCH")
    print("=" * 70)
    for pair in PAIRS:
        desc, bracket, expected_branch = pair[0], pair[1], pair[2]
        oneway = len(pair) > 3 and pair[3] == "oneway"
        if oneway:
            # One-way (branch->bracket only): bracket form should stay unchanged
            got = cabiln_to_branch(bracket)
            ok = got == bracket
            expected_display = bracket + " (unchanged)"
        else:
            got = cabiln_to_branch(bracket)
            ok = got == expected_branch
            expected_display = expected_branch
        if not ok:
            all_pass = False
        status = "PASS" if ok else "FAIL"
        print(f"\n[{status}] {desc}")
        print(f"  input:    {bracket}")
        print(f"  expected: {expected_display}")
        if not ok:
            print(f"  got:      {got}")

    print("\n" + "=" * 70)
    print("BRANCH -> BRACKET (reverse, skip no-convert pairs)")
    print("=" * 70)
    for pair in PAIRS:
        desc, bracket, branch = pair[0], pair[1], pair[2]
        if bracket == branch:
            continue
        got = cabiln_to_bracket(branch)
        ok = got == bracket
        if not ok:
            all_pass = False
        status = "PASS" if ok else "FAIL"
        print(f"\n[{status}] {desc}")
        print(f"  input:    {branch}")
        print(f"  expected: {bracket}")
        if not ok:
            print(f"  got:      {got}")

    print("\n" + "=" * 70)
    print("ALL PASSED" if all_pass else "SOME FAILED")
    print("=" * 70)
