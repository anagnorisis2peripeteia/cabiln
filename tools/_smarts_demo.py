#!/usr/bin/env python3
"""Generate SMARTS matching demo output for README."""
import sys, os, warnings, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
warnings.filterwarnings('ignore')

from pyPept.interfaces.monomer_pipeline import pre_activate

EXAMPLES = [
    # (label, smiles, description)
    # --- Standard amino acids ---
    ("Alanine", "N[C@@H](C)C(=O)O", "backbone only, no sidechain reactive group"),
    ("Cysteine", "N[C@@H](CS)C(=O)O", "thiol sidechain"),
    ("Selenocysteine", "N[C@@H](C[SeH])C(=O)O", "selenol sidechain"),
    ("Lysine", "NCCCC[C@@H](N)C(=O)O", "primary amine sidechain (e-NH2)"),
    ("Aspartate", "N[C@@H](CC(=O)O)C(=O)O", "sidechain carboxyl vs backbone carboxyl"),
    ("Serine", "N[C@@H](CO)C(=O)O", "hydroxyl sidechain"),
    ("Tyrosine", "N[C@@H](Cc1ccc(O)cc1)C(=O)O", "phenolic hydroxyl (label_only)"),
    ("Histidine", "N[C@@H](Cc1c[nH]cn1)C(=O)O", "aromatic NH (label_only)"),
    ("Arginine", "N[C@@H](CCCNC(=N)N)C(=O)O", "guanidinium multi-atom match + amine"),
    ("Proline", "O=C(O)[C@@H]1CCCN1", "ring N - no R3 (only 1 H on N)"),

    # --- Bioconjugation handles ---
    ("Propargylglycine", "N[C@@H](CC#C)C(=O)O", "terminal alkyne (alkyne_c)"),
    ("Azidoalanine", "N[C@@H](CN=[N+]=[N-])C(=O)O", "azide (azide_alpha_c)"),
    ("Chloroacetylalanine", "N[C@@H](CCl)C(=O)O", "alkyl halide (alkyl_halide_c)"),
    ("NHS-alanine", "N[C@@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O", "NHS ester"),
    ("Aldehyde-alanine", "N[C@@H](CC=O)C(=O)O", "aldehyde sidechain"),
    ("Aminooxy-alanine", "N[C@@H](CON)C(=O)O", "aminooxy"),
    ("Hydrazide-alanine", "N[C@@H](CC(=O)NN)C(=O)O", "hydrazide"),
    ("Allylglycine", "N[C@@H](CC=C)C(=O)O", "terminal alkene (RCM stapling)"),

    # --- Edge cases ---
    ("Glycolic acid", "OCC(=O)O", "backbone_o instead of backbone_n (depsipeptide)"),
    ("Acetyl cap", "CC(=O)O", "cap - no N, COOH only -> R2=backbone_c"),
    ("Amide cap", "N", "cap - amine only, no COOH -> R1=backbone_n"),
]

for label, smi, desc in EXAMPLES:
    try:
        result = pre_activate(smi)
        chuckles = result.chuckles
        slots = []
        for slot, ct in sorted(result.chem_types.items()):
            lg = result.leaving.get(slot, '?')
            slots.append(f"R{slot}={ct} (LG={lg})")
        print(f"### {label}")
        print(f"SMILES:   {smi}")
        print(f"CHUCKLES: {chuckles}")
        print(f"Slots:    {', '.join(slots)}")
        print(f"Note:     {desc}")
        print()
    except Exception as e:
        print(f"### {label}")
        print(f"SMILES:   {smi}")
        print(f"ERROR:    {str(e)[:200]}")
        print(f"Note:     {desc}")
        print()
