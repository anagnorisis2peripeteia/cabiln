#!/usr/bin/env python3
"""Check if README bracket and branch forms produce identical molecules."""

import sys
sys.path.insert(0, "src")

from pyPept.sequence import Sequence, cabiln_to_branch, cabiln_to_bracket
from pyPept.molecule import Molecule
from rdkit.Chem import MolToSmiles
import warnings
warnings.filterwarnings("ignore")

def info(label, cabiln):
    try:
        seq = Sequence(cabiln)
        mol = Molecule(seq)
        romol = mol.get_molecule(fmt="ROMol")
        smi = MolToSmiles(romol) if romol else "NO MOL"
        bonds = [(b[0], b[4] if len(b)>4 else '?', b[2], b[5] if len(b)>4 else '?')
                 for b in seq.s_bonds if len(b) > 4]
        print(f"  {label}: {cabiln}")
        print(f"    SMILES: {smi}")
        print(f"    Crosslink bonds: {bonds}")
    except Exception as e:
        print(f"  {label}: {cabiln}")
        print(f"    ERROR: {e}")

# Test the lipid linker forms
print("=== Lipid Linker ===")
info("bracket (2,1)", "ac-K.[gGlu(4,4).AEEA(2,1).C20FA(2,1)]-G-am")
info("branch positional", "ac-K.(4,4)-G-am%gGlu-AEEA(2,1)-C20FA(2,1)")
info("branch !1 reversed", "ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1")

# Converter outputs
print("\n=== Converter ===")
b = "ac-K.[gGlu(4,4).AEEA(2,1).C20FA(2,1)]-G-am"
print(f"  bracket -> branch: {cabiln_to_branch(b)}")
r = "ac-K.!1(4,4)-G-am%C20FA-AEEA-gGlu.!1"
print(f"  branch !1 -> bracket: {cabiln_to_bracket(r)}")
p = "ac-K.(4,4)-G-am%gGlu-AEEA(2,1)-C20FA(2,1)"
print(f"  branch pos -> bracket: {cabiln_to_bracket(p)}")

# Check if C20FA has R1
print("\n=== C20FA R-groups ===")
from rdkit import Chem
sdf_path = "src/pyPept/data/monomers.sdf"
suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
for m in suppl:
    if m and m.GetPropsAsDict().get("m_abbr") == "C20FA":
        props = m.GetPropsAsDict()
        print(f"  m_Rgroups: {props.get('m_Rgroups', '')}")
        print(f"  m_chem_types: {props.get('m_chem_types', '')}")
        # Show dummy atoms
        for atom in m.GetAtoms():
            if atom.GetAtomicNum() == 0:
                print(f"  Dummy R{atom.GetIsotope()} -> {[n.GetSymbol() for n in atom.GetNeighbors()]}")
        break
