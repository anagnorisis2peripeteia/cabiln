"""Audit all 1001 monomers for chem_type counts using pre_smarts patterns."""
import sys
import io
from collections import Counter
from rdkit import Chem
from rdkit.Chem import RWMol
from pyPept.interfaces.reaction_library import _CHEM_TYPE_REGISTRY

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

suppl = Chem.SDMolSupplier("src/pyPept/data/monomers.sdf", removeHs=False)

type_counts = Counter()
type_examples = {}


def restore_raw(mol):
    rg_str = mol.GetProp("m_Rgroups")
    lgs = [x.strip() for x in rg_str.split(",")]

    rw = RWMol(mol)
    dummies = []
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() == 0:
            dummies.append((atom.GetIdx(), atom.GetIsotope()))
    dummies.sort(key=lambda x: x[0], reverse=True)

    for idx, iso in dummies:
        slot = iso - 1
        if slot < 0 or slot >= len(lgs):
            continue
        lg = lgs[slot]
        if lg == "None" or lg == "":
            continue

        atom = rw.GetAtomWithIdx(idx)
        nbrs = [n for n in atom.GetNeighbors()]
        if not nbrs:
            rw.RemoveAtom(idx)
            continue
        nbr = nbrs[0]
        nbr_idx = nbr.GetIdx()

        rw.RemoveAtom(idx)
        # After removal, adjust neighbor index
        adj_nbr = nbr_idx if nbr_idx < idx else nbr_idx

        if lg == "[H]":
            pass  # implicit H
        elif lg == "[OH]":
            new_idx = rw.AddAtom(Chem.Atom(8))
            rw.AddBond(adj_nbr, new_idx, Chem.BondType.SINGLE)
        elif lg == "[Cl]":
            new_idx = rw.AddAtom(Chem.Atom(17))
            rw.AddBond(adj_nbr, new_idx, Chem.BondType.SINGLE)

    try:
        Chem.SanitizeMol(rw)
        return rw.GetMol()
    except Exception:
        return None


# Build patterns from registry
patterns = []
for entry in _CHEM_TYPE_REGISTRY:
    ctype, pre_sma, lg, infer_sma, label_only = entry
    if pre_sma:
        pat = Chem.MolFromSmarts(pre_sma)
        if pat:
            patterns.append((ctype, pat, label_only))

count = 0
failed = 0
fail_names = []
for mol in suppl:
    if mol is None:
        continue
    name = mol.GetProp("symbol")
    raw = restore_raw(mol)
    if raw is None:
        failed += 1
        fail_names.append(name)
        continue

    count += 1
    raw_h = Chem.AddHs(raw)

    matched_types = set()
    for ctype, pat, label_only in patterns:
        if raw_h.HasSubstructMatch(pat):
            if ctype not in matched_types:
                matched_types.add(ctype)
                type_counts[ctype] += 1
                if ctype not in type_examples:
                    type_examples[ctype] = []
                if len(type_examples[ctype]) < 10:
                    type_examples[ctype].append(name)

print(f"Processed: {count}, failed restore: {failed}")
if fail_names:
    print(f"Failed: {', '.join(fail_names[:20])}")
print()
print("=== Pre-activate chem_type counts (from restored raw SMILES) ===")
for typ in sorted(type_counts, key=type_counts.get, reverse=True):
    ex = ", ".join(type_examples.get(typ, [])[:6])
    lo = ""
    for entry in _CHEM_TYPE_REGISTRY:
        if entry[0] == typ and entry[4]:
            lo = " [LABEL_ONLY]"
    print(f"{typ}{lo}: {type_counts[typ]}  ex: {ex}")
