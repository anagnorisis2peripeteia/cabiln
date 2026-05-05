#!/usr/bin/env python3
"""Test degenerate caps at both N-term and C-term positions."""
import sys
sys.path.insert(0, "src")
import warnings
warnings.filterwarnings("ignore")

from pyPept.sequence import Sequence
from pyPept.molecule import Molecule
from rdkit.Chem import MolToSmiles


def smiles_for(label, cabiln):
    try:
        seq = Sequence(cabiln)
        mol = Molecule(seq)
        romol = mol.get_molecule(fmt="ROMol")
        smi = MolToSmiles(romol) if romol else "NO MOL"
        print(f"  {label}: {smi[:120]}")
        return smi
    except Exception as e:
        print(f"  {label}: FAIL - {e}")
        return None


print("=== Bn degenerate vs paired ===")
s1 = smiles_for("Bn  at N-term", "Bn-A-G-am")
s2 = smiles_for("Bn_ at N-term", "Bn_-A-G-am")
if s1 and s2:
    print(f"  N-term MATCH: {s1 == s2}")

s3 = smiles_for("Bn  at C-term", "ac-A-G-Bn")
s4 = smiles_for("_Bn at C-term", "ac-A-G-_Bn")
if s3 and s4:
    print(f"  C-term MATCH: {s3 == s4}")

print("\n=== Me degenerate vs paired ===")
s5 = smiles_for("Me  at N-term", "Me-A-G-am")
s6 = smiles_for("Me_ at N-term", "Me_-A-G-am")
if s5 and s6:
    print(f"  N-term MATCH: {s5 == s6}")

s7 = smiles_for("Me  at C-term", "ac-A-G-Me")
s8 = smiles_for("_Me at C-term", "ac-A-G-_Me")
if s7 and s8:
    print(f"  C-term MATCH: {s7 == s8}")

print("\n=== Et degenerate vs paired ===")
s9 = smiles_for("Et  at N-term", "Et-A-G-am")
s10 = smiles_for("Et_ at N-term", "Et_-A-G-am")
if s9 and s10:
    print(f"  N-term MATCH: {s9 == s10}")

s11 = smiles_for("Et  at C-term", "ac-A-G-Et")
s12 = smiles_for("_Et at C-term", "ac-A-G-_Et")
if s11 and s12:
    print(f"  C-term MATCH: {s11 == s12}")

print("\n=== NHBn degenerate vs paired ===")
s13 = smiles_for("NHBn  at N-term", "NHBn-A-G-am")
s14 = smiles_for("NHBn_ at N-term", "NHBn_-A-G-am")
if s13 and s14:
    print(f"  N-term MATCH: {s13 == s14}")

s15 = smiles_for("NHBn  at C-term", "ac-A-G-NHBn")
s16 = smiles_for("_NHBn at C-term", "ac-A-G-_NHBn")
if s15 and s16:
    print(f"  C-term MATCH: {s15 == s16}")

print("\n=== Both-ends test ===")
smiles_for("Bn-A-G-K-Bn", "Bn-A-G-K-Bn")
smiles_for("Me-A-G-Me",   "Me-A-G-Me")
