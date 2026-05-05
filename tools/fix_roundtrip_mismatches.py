#!/usr/bin/env python3
"""
Fix all round-trip mismatches by:
1. Targeted LG/slot fixes for caps with wrong reagent-form LGs or slot direction
2. Accept pipeline output for monomers where pre_activate finds valid extra sites
3. Report remaining ambiguous cases
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem, SDWriter, rdDepictor
from pyPept.interfaces.monomer_pipeline import pre_activate

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')

df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')
print(f"Loaded {len(df)} monomers\n")

# ── Phase 1: Targeted SDF data fixes ────────────────────────────────────────
# These fix the stored data so that restore_full_smiles produces unambiguous
# input for pre_activate.

TARGETED_FIXES = {
    # Alkylation caps: store reagent-form LG so restored SMILES has halide
    'Me_':  {'lg': {2: '[I]'}},       # MeI
    'Et_':  {'lg': {2: '[Br]'}},      # EtBr
    'Bn_':  {'lg': {2: '[Br]'}},      # BnBr
    # Protecting groups: fix slot direction (electrophilic → R2)
    'boc':  {'rename_slot': {1: 2}},   # COOH cap → R2
    'pbf':  {'rename_slot': {1: 2}},   # sulfonyl cap → R2
    # acm: alkylation cap, fix both slot and LG
    'acm':  {'rename_slot': {1: 2}, 'lg': {2: '[Cl]'}},
    # Bioconjugation handles: rename R4→R3 to match pipeline sequential numbering
    'Mal':  {'rename_slot': {4: 3}},
    'TCO':  {'rename_slot': {4: 3}},
    'Tz':   {'rename_slot': {4: 3}},
}

def parse_rgroups(rg):
    if isinstance(rg, list): return rg
    return [None if v.strip() == 'None' else v.strip() for v in rg.split(',')]

def rgroups_to_str(rg_list):
    return ', '.join(v if v else 'None' for v in rg_list)

def apply_targeted_fix(sym, mol, rg_list, fix):
    """Apply a targeted fix to a monomer's SDF entry."""
    emol = Chem.RWMol(mol)
    modified = False

    if 'rename_slot' in fix:
        for old_slot, new_slot in fix['rename_slot'].items():
            for atom in emol.GetAtoms():
                if atom.GetAtomicNum() == 0 and atom.GetIsotope() == old_slot:
                    atom.SetIsotope(new_slot)
                    modified = True
            # Remap LG: move old slot's LG to new slot position
            if old_slot <= len(rg_list) and rg_list[old_slot - 1] is not None:
                lg_val = rg_list[old_slot - 1]
                # Extend list if needed
                while len(rg_list) < new_slot:
                    rg_list.append(None)
                rg_list[new_slot - 1] = lg_val
                rg_list[old_slot - 1] = None
                modified = True

    if 'lg' in fix:
        for slot, new_lg in fix['lg'].items():
            while len(rg_list) < slot:
                rg_list.append(None)
            rg_list[slot - 1] = new_lg
            modified = True

    # Trim trailing None entries
    while rg_list and rg_list[-1] is None:
        rg_list.pop()

    if modified:
        try:
            Chem.SanitizeMol(emol)
        except:
            pass
    return emol.GetMol(), rg_list, modified


# ── Phase 2: Restore + pre_activate + compare ───────────────────────────────

def restore_full_smiles(mol, rgroups):
    emol = Chem.RWMol(Chem.RWMol(mol))
    try:
        Chem.Kekulize(emol, clearAromaticFlags=True)
    except:
        pass
    dummies = [(a.GetIdx(), a.GetIsotope()) for a in emol.GetAtoms() if a.GetAtomicNum() == 0]
    rm = []
    for d, s in dummies:
        if s < 1 or s > len(rgroups): continue
        lg = rgroups[s - 1]
        if lg is None: continue
        if lg == '[H]':
            for nb in emol.GetAtomWithIdx(d).GetNeighbors():
                nb.SetNoImplicit(False)
            rm.append(d)
        elif lg == '[OH]':
            emol.GetAtomWithIdx(d).SetAtomicNum(8)
            emol.GetAtomWithIdx(d).SetIsotope(0)
            emol.GetAtomWithIdx(d).SetNumExplicitHs(1)
        elif lg in ('[Cl]', '[Br]', '[I]'):
            emol.GetAtomWithIdx(d).SetAtomicNum({'[Cl]': 17, '[Br]': 35, '[I]': 53}[lg])
            emol.GetAtomWithIdx(d).SetIsotope(0)
        else:
            rm.append(d)
    for i in sorted(rm, reverse=True):
        emol.RemoveAtom(i)
    try:
        Chem.SanitizeMol(emol)
        return Chem.MolToSmiles(emol)
    except:
        return None


def canonical(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m else None


def strip_dummies(smi):
    """Remove all dummy atoms from a SMILES — used for safety check."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return None
    emol = Chem.RWMol(mol)
    for idx in sorted([a.GetIdx() for a in emol.GetAtoms() if a.GetAtomicNum() == 0], reverse=True):
        for nb in emol.GetAtomWithIdx(idx).GetNeighbors():
            nb.SetNoImplicit(False)
        emol.RemoveAtom(idx)
    try:
        Chem.SanitizeMol(emol)
        return Chem.MolToSmiles(emol)
    except:
        return None


# ── Main loop ────────────────────────────────────────────────────────────────

phase1_fixed = []
phase2_fixed = []
still_failing = []
total_mismatch = 0

for sym in sorted(df.index):
    mol = df.loc[sym, 'ROMol']
    if mol is None: continue
    rg = parse_rgroups(df.loc[sym, 'm_Rgroups'])
    stored_chuckles = canonical(Chem.MolToSmiles(mol))

    # Phase 1: apply targeted fix if applicable
    if sym in TARGETED_FIXES:
        mol, rg, _ = apply_targeted_fix(sym, mol, rg, TARGETED_FIXES[sym])
        df.at[sym, 'ROMol'] = mol
        df.at[sym, 'm_Rgroups'] = rgroups_to_str(rg)
        stored_chuckles = canonical(Chem.MolToSmiles(mol))

    # Restore full SMILES
    full_smi = restore_full_smiles(mol, rg)
    if full_smi is None: continue

    # Run pre_activate
    try:
        result = pre_activate(full_smi)
    except:
        continue

    result_canon = canonical(result.chuckles)
    stored_lgs = {i + 1: l for i, l in enumerate(rg) if l is not None}

    if stored_chuckles == result_canon and result.leaving == stored_lgs:
        if sym in TARGETED_FIXES:
            phase1_fixed.append(sym)
        continue

    total_mismatch += 1

    # Phase 2: safety check — restoring LGs on pipeline output gives same molecule?
    result_mol_tmp = Chem.MolFromSmiles(result.chuckles)
    result_rg_tmp = [None] * max(result.leaving.keys())
    for s, lg in result.leaving.items():
        result_rg_tmp[s - 1] = lg
    result_restored = restore_full_smiles(result_mol_tmp, result_rg_tmp) if result_mol_tmp else None
    original_restored = canonical(full_smi)

    if result_restored and canonical(result_restored) == original_restored:
        # Same underlying molecule — safe to accept pipeline output
        new_mol = Chem.MolFromSmiles(result.chuckles)
        if new_mol is not None:
            rdDepictor.Compute2DCoords(new_mol)

            # Set m_Rgroups from pipeline leaving dict
            max_slot = max(result.leaving.keys())
            new_rg = [None] * max_slot
            for s, lg in result.leaving.items():
                new_rg[s - 1] = lg

            df.at[sym, 'ROMol'] = new_mol
            df.at[sym, 'm_Rgroups'] = rgroups_to_str(new_rg)
            phase2_fixed.append(sym)
            continue

    still_failing.append((sym, stored_chuckles, result_canon,
                          str(stored_lgs), str(result.leaving)))


print("=" * 80)
print("RESULTS")
print("=" * 80)
print(f"  Phase 1 (targeted fix):     {len(phase1_fixed)} — {phase1_fixed}")
print(f"  Phase 2 (accept pipeline):  {len(phase2_fixed)} — {phase2_fixed}")
print(f"  Still failing:              {len(still_failing)}")
print(f"  Total mismatches found:     {total_mismatch}")

if still_failing:
    print(f"\n--- Still failing ({len(still_failing)}) ---")
    for sym, st, res, slg, rlg in still_failing:
        print(f"  {sym:25s}")
        print(f"    stored:  {st[:70]}")
        print(f"    result:  {res[:70]}")
        print(f"    s_lgs:   {slg}")
        print(f"    r_lgs:   {rlg}")

# ── Write updated SDF ────────────────────────────────────────────────────────
if phase1_fixed or phase2_fixed:
    print(f"\nWriting updated SDF with {len(phase1_fixed) + len(phase2_fixed)} fixes...")
    writer = SDWriter(SDF_PATH)
    for sym in df.index:
        mol = df.loc[sym, 'ROMol']
        if mol is None: continue
        mol.SetProp('symbol', sym)
        mol.SetProp('m_Rgroups', df.loc[sym, 'm_Rgroups'])
        for col in df.columns:
            if col not in ('ROMol', 'symbol', 'm_Rgroups', 'ID'):
                val = df.loc[sym, col]
                if val is not None and str(val) != 'nan':
                    mol.SetProp(col, str(val))
        writer.write(mol)
    writer.close()
    print("Done.")
