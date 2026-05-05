import sys
sys.path.insert(0, "src")
from pyPept.sequence import cabiln_to_branch, cabiln_to_bracket

# (1,2) bracket should stay as bracket
r1 = cabiln_to_branch("K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-G-am")
print("(1,2) bracket ->", r1)
assert "[" in r1, "Should stay as bracket!"

# (2,1) bracket should convert to N->C branch
r2 = cabiln_to_branch("ac-A-K.[G(4,1).A(2,1).am(2,1)]-G-am")
print("(2,1) bracket ->", r2)
assert "%" in r2, "Should convert to branch!"

# roundtrip for (2,1)
r3 = cabiln_to_bracket(r2)
print("roundtrip ->", r3)
assert r3 == "ac-A-K.[G(4,1).A(2,1).am(2,1)]-G-am", f"roundtrip failed: {r3}"

# crosslink branch -> bracket still works
r4 = cabiln_to_bracket("ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1")
print("xlink branch -> bracket:", r4)

# retatrutide bracket stays as bracket
reta = ("Y-Aib-Q-G-T-F-T-S-D-Y-S-I-aMeLeu-L-D-"
        "K.[gGlu(4,4).AEEA(1,2).C20FA(1,2)]-A-Q-Aib-A-F-I-E-"
        "Y-L-L-E-G-G-P-S-S-G-A-P-P-P-S-am")
r5 = cabiln_to_branch(reta)
print("retatrutide ->", r5[:80])
assert "[" in r5, "Retatrutide should keep bracket!"
assert "%" not in r5, "Retatrutide should NOT have %!"

print("\nALL OK")
