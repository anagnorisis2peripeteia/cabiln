#!/usr/bin/env python3
"""
Populate m_chem_types (and fix m_Rgroups) for all monomers missing the field.
Flags selection-rule conflicts to stdout.
"""

import sys, re, pathlib, collections
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / 'src'))

from rdkit import Chem
from rdkit.Chem import SDWriter, rdDepictor
from pyPept.interfaces.monomer_pipeline import pre_activate, ActivationError

SDF  = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
rdDepictor.SetPreferCoordGen(True)


def restore_plain(chuckles, rgroups_str):
    rgroups = [r.strip() for r in rgroups_str.split(',')]
    s = chuckles
    for i, lg in enumerate(rgroups):
        slot = i + 1
        if lg.lower() not in ('none', ''):
            s = re.sub(r'\[' + str(slot) + r'\*\]', lg, s)
    return s


def chem_types_str(ct_dict):
    return ','.join(f'{k}:{v}' for k, v in sorted(ct_dict.items(), key=lambda x: int(x[0])))


def rgroups_str(leaving_dict):
    lst = ['None'] * 6
    for slot, lg in leaving_dict.items():
        idx = int(slot) - 1
        if 0 <= idx < 6:
            lst[idx] = lg
    return ','.join(lst)


# ── conflict rules ────────────────────────────────────────────────────────────
def check_conflicts(abbr, m_type, detected_ct):
    issues = []
    types_set = set(detected_ct.values())

    # backbone_o is the C-terminal analog in amino alcohols (_ol series)
    has_bb_c = 'backbone_c' in types_set or 'backbone_o' in types_set or 'lactone_c' in types_set
    # cyclooctyne N is ring-embedded — backbone_n SMARTS won't fire but N is present
    has_bb_n = 'backbone_n' in types_set or 'cyclooctyne_c' in types_set
    has_backbone = has_bb_n and has_bb_c

    if m_type == 'aa':
        # Full backbone required; _ol AAs count backbone_o as their C-terminus
        if not has_backbone:
            issues.append('m_type=aa but no complete backbone (N+C/O) detected')
        # Partial backbone on an aa is always a problem
        if has_bb_n and not has_bb_c:
            issues.append('backbone_n detected without backbone_c or backbone_o')
        if has_bb_c and not has_bb_n:
            issues.append('backbone_c/o detected without backbone_n')

    if m_type == 'cap':
        # Caps have exactly ONE terminus by design; having BOTH is suspicious
        if has_backbone:
            issues.append('m_type=cap but full backbone (N+C) detected — may be miscategorized')
        # Single-terminus caps (N-only or C-only) are correct — no partial-backbone check

    # linker: omega-amino acids legitimately have full backbone — no check

    # Slot count sanity
    if len(detected_ct) > 7:
        issues.append(f'{len(detected_ct)} R-groups detected — unusually many, check structure')

    return issues


# ── main ──────────────────────────────────────────────────────────────────────
suppl   = Chem.SDMolSupplier(str(SDF), removeHs=False)
all_mols = list(suppl)

updated   = 0
skipped   = 0
failed    = []
conflicts = []

print(f'Total entries in SDF: {len(all_mols)}')
print()

for i, mol in enumerate(all_mols):
    if mol is None:
        continue
    p    = mol.GetPropsAsDict()
    abbr = p.get('m_abbr', '?')
    ct   = p.get('m_chem_types', '').strip()
    rg   = p.get('m_Rgroups', '').strip()
    mtype = p.get('m_type', '')

    if ct:  # already has chem_types
        skipped += 1
        continue

    chuckles = Chem.MolToSmiles(mol)
    plain    = restore_plain(chuckles, rg)

    try:
        result = pre_activate(plain)

        # Conflict checks
        issues = check_conflicts(abbr, mtype, result.chem_types)
        if issues:
            conflicts.append((abbr, mtype, result.chem_types, issues))

        # Build new mol from CHUCKLES output
        new_mol = Chem.MolFromSmiles(result.chuckles)
        if new_mol is None:
            failed.append((abbr, 'invalid CHUCKLES returned by pre_activate'))
            continue

        rdDepictor.Compute2DCoords(new_mol)

        # Copy all existing properties
        for k, v in p.items():
            new_mol.SetProp(k, str(v))

        # Overwrite / add the pipeline fields
        new_mol.SetProp('m_chem_types', chem_types_str(result.chem_types))
        new_mol.SetProp('m_Rgroups',    rgroups_str(result.leaving))

        all_mols[i] = new_mol
        updated += 1

    except ActivationError as e:
        failed.append((abbr, str(e)))
    except Exception as e:
        failed.append((abbr, f'{type(e).__name__}: {e}'))

# ── write SDF ─────────────────────────────────────────────────────────────────
with open(SDF, 'w') as fh:
    writer = SDWriter(fh)
    writer.SetKekulize(False)
    for mol in all_mols:
        if mol is not None:
            writer.write(mol)
    writer.close()

# ── report ────────────────────────────────────────────────────────────────────
print(f'Updated:  {updated}')
print(f'Skipped (already had chem_types): {skipped}')
print(f'Failed:   {len(failed)}')
print(f'Conflicts flagged: {len(conflicts)}')
print()

if failed:
    print('--- FAILURES (pre_activate could not process) ---')
    fc = collections.Counter(e.split(':')[0] for _, e in failed)
    for msg, n in fc.most_common():
        print(f'  [{n}x] {msg}')
    print('  First 20:')
    for abbr, e in failed[:20]:
        print(f'    {abbr}: {e}')
    print()

if conflicts:
    print('--- CONFLICTS (selection rule issues) ---')
    for abbr, mtype, ct, issues in conflicts:
        print(f'  {abbr} ({mtype}):')
        for iss in issues:
            print(f'    ! {iss}')
        print(f'    detected: {ct}')
