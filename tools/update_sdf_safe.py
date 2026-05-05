#!/usr/bin/env python3
"""
Safe SDF update: only update monomers where the pipeline adds genuinely new
R-group slots WITHOUT changing the existing slot structure.

Skips:
- Hydrocarbon caps that can't be auto-detected
- COOH placement changes (C([n*])=O → C(=O)O[n*])
- Protecting group caps with slot reassignment
- Crosslinkers with slot renumbering
"""
import sys, os, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem
from pyPept.interfaces.monomer_pipeline import pre_activate

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'src', 'pyPept', 'data')
SDF_PATH = os.path.join(DATA_DIR, 'monomers.sdf')

print("Loading SDF...")
df = PandasTools.LoadSDF(SDF_PATH)
df = df.set_index('symbol')
print(f"  {len(df)} monomers\n")

SKIP_CAPS = {'Dmb', 'Ph', '_Bn', '_Et', '_Me', 'otbu'}

LG_ATOMS = {
    '[H]': (1, 0), '[OH]': (8, 1),
    '[Cl]': (17, 0), '[Br]': (35, 0), '[I]': (53, 0),
}

def restore_full_smiles(mol, rgroups):
    emol = Chem.RWMol(Chem.RWMol(mol))
    try:
        Chem.Kekulize(emol, clearAromaticFlags=True)
    except Exception:
        pass
    dummies = []
    for atom in emol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            dummies.append((atom.GetIdx(), atom.GetIsotope()))
    replacements = []
    for dummy_idx, slot in dummies:
        if slot < 1 or slot > len(rgroups):
            continue
        lg = rgroups[slot - 1]
        if lg is None:
            continue
        if lg == '[H]':
            replacements.append((dummy_idx, 'remove', None))
        elif lg == '[OH]':
            replacements.append((dummy_idx, 'replace_oh', None))
        elif lg in ('[Cl]', '[Br]', '[I]'):
            replacements.append((dummy_idx, 'replace_halide', LG_ATOMS[lg][0]))
        else:
            replacements.append((dummy_idx, 'unknown', lg))
    atoms_to_remove = []
    for dummy_idx, action, param in replacements:
        if action == 'remove':
            for nb in emol.GetAtomWithIdx(dummy_idx).GetNeighbors():
                nb.SetNoImplicit(False)
            atoms_to_remove.append(dummy_idx)
        elif action == 'replace_oh':
            emol.GetAtomWithIdx(dummy_idx).SetAtomicNum(8)
            emol.GetAtomWithIdx(dummy_idx).SetIsotope(0)
            emol.GetAtomWithIdx(dummy_idx).SetNumExplicitHs(1)
        elif action == 'replace_halide':
            emol.GetAtomWithIdx(dummy_idx).SetAtomicNum(param)
            emol.GetAtomWithIdx(dummy_idx).SetIsotope(0)
        elif action == 'unknown':
            atoms_to_remove.append(dummy_idx)
    for idx in sorted(atoms_to_remove, reverse=True):
        emol.RemoveAtom(idx)
    try:
        Chem.SanitizeMol(emol)
        return Chem.MolToSmiles(emol)
    except Exception:
        return None

def parse_rgroups(rg):
    if isinstance(rg, list):
        return rg
    return [None if v.strip() == 'None' else v.strip() for v in rg.split(',')]

def canonical(smi):
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol) if mol else None

def count_slots(smiles):
    return len(re.findall(r'\[\d+\*\]', smiles))

def get_slot_numbers(smiles):
    return sorted(set(int(m) for m in re.findall(r'\[(\d+)\*\]', smiles)))

def is_safe_update(stored_canon, result_canon, stored_lgs, result_leaving):
    """Check if this update only ADDS new slots without changing existing ones."""
    stored_slots = get_slot_numbers(stored_canon)
    result_slots = get_slot_numbers(result_canon)

    # All existing slots must still be present
    for s in stored_slots:
        if s not in result_slots:
            return False, f"slot {s} removed"

    # Existing LGs must match
    for slot, lg in stored_lgs.items():
        if slot in result_leaving and result_leaving[slot] != lg:
            return False, f"slot {slot} LG changed: {lg} -> {result_leaving[slot]}"

    # Strip all dummies from both CHUCKLES and compare the backbone
    stored_bare = re.sub(r'\[\d+\*\]', '[*]', stored_canon)
    result_bare = re.sub(r'\[\d+\*\]', '[*]', result_canon)

    # Check for COOH placement change: C([*])=O vs C(=O)O[*]
    # This is unsafe — changes where bonds form
    stored_cooh_on_c = len(re.findall(r'C\(\[\*\]\)=O', stored_bare))
    result_cooh_on_o = len(re.findall(r'C\(=O\)O\[\*\]', result_bare))
    if stored_cooh_on_c > 0 and result_cooh_on_o > 0:
        return False, "COOH placement change (C→O)"

    # If same slot count but different structure, check if it's just renumbering
    if len(stored_slots) == len(result_slots) and stored_slots != result_slots:
        return False, "slot renumbering"

    if stored_canon == result_canon:
        return False, "no change needed"

    return True, "new slots added"

updated = []
skipped = []
already_match = []
unsafe = []
failed = []

for sym in sorted(df.index):
    if sym in SKIP_CAPS:
        skipped.append((sym, 'hydrocarbon cap'))
        continue

    mol = df.loc[sym, 'ROMol']
    if mol is None:
        continue

    rgroups = parse_rgroups(df.loc[sym, 'm_Rgroups'])
    stored_chuckles = Chem.MolToSmiles(mol)

    full_smi = restore_full_smiles(mol, rgroups)
    if full_smi is None:
        failed.append((sym, 'restore failed'))
        continue

    try:
        result = pre_activate(full_smi)
    except Exception as e:
        failed.append((sym, str(e)[:40]))
        continue

    stored_canon = canonical(stored_chuckles)
    result_canon = canonical(result.chuckles)

    stored_lgs = {}
    for i, lg in enumerate(rgroups, 1):
        if lg is not None:
            stored_lgs[i] = lg

    if stored_canon == result_canon and result.leaving == stored_lgs:
        already_match.append(sym)
        continue

    safe, reason = is_safe_update(stored_canon, result_canon, stored_lgs, result.leaving)

    if not safe:
        unsafe.append((sym, reason))
        continue

    # Safe to update — create new mol from result CHUCKLES
    new_mol = Chem.MolFromSmiles(result.chuckles)
    if new_mol is None:
        failed.append((sym, f'cant parse: {result.chuckles[:40]}'))
        continue

    AllChem.Compute2DCoords(new_mol)

    max_slot = max(result.leaving.keys()) if result.leaving else 0
    max_slot = max(max_slot, 6)
    rg_parts = []
    for i in range(1, max_slot + 1):
        if i in result.leaving:
            rg_parts.append(result.leaving[i])
        else:
            rg_parts.append('None')
    new_rgroups = ','.join(rg_parts)

    df.at[sym, 'ROMol'] = new_mol
    df.at[sym, 'm_Rgroups'] = new_rgroups
    updated.append(sym)

print(f"Already match:    {len(already_match)}")
print(f"Safe updates:     {len(updated)}")
print(f"Unsafe (skipped): {len(unsafe)}")
print(f"Skipped caps:     {len(skipped)}")
print(f"Failed:           {len(failed)}")
print(f"Total:            {len(df)}")

if updated:
    print(f"\n--- Safely updated ({len(updated)}) ---")
    for sym in updated:
        rg = df.loc[sym, 'm_Rgroups']
        smi = Chem.MolToSmiles(df.loc[sym, 'ROMol'])
        print(f"  {sym:20s} {smi[:50]}  rg={rg[:40]}")

if unsafe:
    print(f"\n--- Unsafe (not updated) ({len(unsafe)}) ---")
    by_reason = {}
    for sym, reason in unsafe:
        by_reason.setdefault(reason, []).append(sym)
    for reason, syms in sorted(by_reason.items()):
        print(f"  {reason}: {len(syms)} monomers")
        for s in syms[:3]:
            print(f"    {s}")
        if len(syms) > 3:
            print(f"    ... and {len(syms)-3} more")

if failed:
    print(f"\n--- Failed ({len(failed)}) ---")
    for sym, err in failed:
        print(f"  {sym}: {err}")

# Save
print(f"\nSaving SDF...")
df_save = df.copy().reset_index()
df_save = df_save.rename(columns={'index': 'symbol'})
PandasTools.WriteSDF(df_save, SDF_PATH, molColName='ROMol',
                     properties=list(df_save.columns.drop(['ROMol'])))
print("  Done")
