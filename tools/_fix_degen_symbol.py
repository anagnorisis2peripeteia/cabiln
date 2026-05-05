#!/usr/bin/env python3
"""Fix degenerate caps: remove entries without symbol, re-add with symbol set."""
import sys
sys.path.insert(0, "src")
from rdkit import Chem
from rdkit.Chem import AllChem

SDF_PATH = "src/pyPept/data/monomers.sdf"
DEGEN_ABBRS = {"Bn", "Me", "Et", "NHBn"}

DEGEN_CAPS = [
    ("Bn",   "[1*]C([2*])c1ccccc1",    "Benzyl (degenerate)"),
    ("Me",   "[1*][C]([2*])([H])[H]",  "Methyl (degenerate)"),
    ("Et",   "[1*]C([2*])C",            "Ethyl (degenerate)"),
    ("NHBn", "[1*]N([2*])Cc1ccccc1",    "NH-Benzyl (degenerate)"),
]

print("Loading SDF...")
suppl = Chem.SDMolSupplier(SDF_PATH, removeHs=False)
existing = []
for m in suppl:
    if m is None:
        continue
    abbr = m.GetPropsAsDict().get("m_abbr", "")
    if abbr in DEGEN_ABBRS:
        print(f"  Removing old degenerate entry: {abbr}")
        continue
    existing.append(m)

print(f"  Kept {len(existing)} monomers")

print("\nCreating new degenerate entries...")
new_mols = []
for abbr, smiles, name in DEGEN_CAPS:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  {abbr}: FAILED to parse {smiles}")
        continue
    mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)

    mol.SetProp("symbol", abbr)
    mol.SetProp("_Name", abbr)
    mol.SetProp("m_abbr", abbr)
    mol.SetProp("m_name", name)
    mol.SetProp("m_type", "cap")
    mol.SetProp("m_subtype", "cap")
    mol.SetProp("m_Rgroups", "[H],[H],None,None,None,None")
    mol.SetProp("m_degenerate", "true")

    # Verify isotopes
    isos = [a.GetIsotope() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    print(f"  {abbr}: OK, isotopes={sorted(isos)}, atoms={mol.GetNumAtoms()}")
    new_mols.append(mol)

print(f"\nWriting SDF with {len(existing) + len(new_mols)} total monomers...")
writer = Chem.SDWriter(SDF_PATH)
writer.SetForceV3000(False)
for m in existing:
    writer.write(m)
for m in new_mols:
    writer.write(m)
writer.close()

# Verify
suppl2 = Chem.SDMolSupplier(SDF_PATH, removeHs=False)
count = sum(1 for m in suppl2 if m)
print(f"  Done: {count} monomers in SDF")

# Check symbol
from rdkit.Chem import PandasTools
df = PandasTools.LoadSDF(SDF_PATH)
for abbr in DEGEN_ABBRS:
    row = df[df['m_abbr'] == abbr]
    if not row.empty:
        print(f"  {abbr}: symbol={row.iloc[0].get('symbol', 'MISSING')}")
