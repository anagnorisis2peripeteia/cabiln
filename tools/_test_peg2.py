import sys
sys.path.insert(0, "src")
from rdkit import Chem
suppl = Chem.SDMolSupplier("src/pyPept/data/monomers.sdf", removeHs=False)
for m in suppl:
    if m and m.GetPropsAsDict().get("m_abbr") == "PEG":
        p = m.GetPropsAsDict()
        print("Found PEG")
        print("  Rgroups:", p.get("m_Rgroups"))
        print("  chem_types:", p.get("m_chem_types"))
        break
else:
    print("PEG NOT FOUND")
