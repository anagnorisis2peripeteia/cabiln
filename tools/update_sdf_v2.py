#!/usr/bin/env python3
"""
Safe SDF update v2: only update monomers where the pipeline's result
contains extra slots AND removing those extra slots yields exactly
the stored CHUCKLES. This guarantees no structural changes.
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

def get_slot_numbers(smiles):
    return sorted(set(int(m) for m in re.findall(r'\[(\d+)\*\]', smiles)))

def strip_extra_slots(result_chuckles, stored_slots):
    """Remove extra dummy atoms from result CHUCKLES, keeping only stored slots.
    Returns canonical SMILES or None on failure."""
    mol = Chem.MolFromSmiles(result_chuckles)
    if mol is None:
        return None
    emol = Chem.RWMol(mol)

    to_remove = []
    for atom in emol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            iso = atom.GetIsotope()
            if iso not in stored_slots:
                for nb in atom.GetNeighbors():
                    nb.SetNoImplicit(False)
                to_remove.append(atom.GetIdx())

    for idx in sorted(to_remove, reverse=True):
        emol.RemoveAtom(idx)

    try:
        Chem.SanitizeMol(emol)
        return Chem.MolToSmiles(emol)
    except Exception:
        return None

updated = []
skipped_caps = []
already_match = []
not_safe = []
failed = []

for sym in sorted(df.index):
    if sym in SKIP_CAPS:
        skipped_caps.append(sym)
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
        failed.append((sym, str(e)[:60]))
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

    # Safety check: result must have MORE slots
    stored_slots = set(get_slot_numbers(stored_canon))
    result_slots = set(get_slot_numbers(result_canon))
    extra_slots = result_slots - stored_slots

    if not extra_slots:
        not_safe.append((sym, 'no new slots'))
        continue

    # All stored slots must exist in result
    if not stored_slots.issubset(result_slots):
        not_safe.append((sym, 'lost stored slots'))
        continue

    # Existing LGs must match
    lg_mismatch = False
    for slot, lg in stored_lgs.items():
        if slot in result.leaving and result.leaving[slot] != lg:
            not_safe.append((sym, f'LG changed for slot {slot}'))
            lg_mismatch = True
            break
    if lg_mismatch:
        continue

    # KEY CHECK: strip extra dummies from result → must match stored CHUCKLES
    stripped = strip_extra_slots(result.chuckles, stored_slots)
    if stripped is None:
        not_safe.append((sym, 'strip failed'))
        continue

    if stripped != stored_canon:
        not_safe.append((sym, 'structure changed'))
        continue

    # All checks passed — safe to update
    new_mol = Chem.MolFromSmiles(result.chuckles)
    if new_mol is None:
        failed.append((sym, f'parse failed: {result.chuckles[:40]}'))
        continue

    AllChem.Compute2DCoords(new_mol)

    max_slot = max(result.leaving.keys()) if result.leaving else 0
    max_slot = max(max_slot, 6)
    rg_parts = []
    for i in range(1, max_slot + 1):
        rg_parts.append(result.leaving.get(i, 'None'))
    new_rgroups = ','.join(rg_parts)

    df.at[sym, 'ROMol'] = new_mol
    df.at[sym, 'm_Rgroups'] = new_rgroups
    updated.append(sym)

print(f"Already match:     {len(already_match)}")
print(f"Safely updated:    {len(updated)}")
print(f"Not safe (kept):   {len(not_safe)}")
print(f"Skipped caps:      {len(skipped_caps)}")
print(f"Failed:            {len(failed)}")
print(f"Total:             {len(df)}")

if not_safe:
    print(f"\n--- Not safe ({len(not_safe)}) ---")
    by_reason = {}
    for sym, reason in not_safe:
        by_reason.setdefault(reason, []).append(sym)
    for reason in sorted(by_reason):
        syms = by_reason[reason]
        print(f"  {reason}: {len(syms)} monomers")
        for s in syms[:5]:
            print(f"    {s}")
        if len(syms) > 5:
            print(f"    ... and {len(syms)-5} more")

if failed:
    print(f"\n--- Failed ({len(failed)}) ---")
    for sym, err in failed:
        print(f"  {sym}: {err}")

# Save
if updated:
    print(f"\nSaving SDF with {len(updated)} updates...")
    df_save = df.copy().reset_index()
    df_save = df_save.rename(columns={'index': 'symbol'})
    PandasTools.WriteSDF(df_save, SDF_PATH, molColName='ROMol',
                         properties=list(df_save.columns.drop(['ROMol'])))
    print("  Done")
else:
    print("\nNo safe updates to apply.")
