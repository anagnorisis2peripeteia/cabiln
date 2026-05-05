#!/usr/bin/env python3
"""
Remove degenerate cap entries (m_degenerate=true) from monomers.sdf.

These entries have two dummy atoms and m_degenerate=true. They are being replaced
by automatic degeneracy detection that groups paired caps (e.g. Bn_/_Bn) by their
core structure and resolves the base name (Bn) to the correct variant at parse time.
"""

import pathlib
from rdkit import Chem


def main():
    sdf_path = pathlib.Path(__file__).parent.parent / 'src' / 'pyPept' / 'data' / 'monomers.sdf'
    print(f"Loading SDF: {sdf_path}")

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    kept = []
    removed = []

    for mol in suppl:
        if mol is None:
            continue
        props = mol.GetPropsAsDict()
        if props.get('m_degenerate') == 'true':
            removed.append(props.get('m_abbr', props.get('symbol', '???')))
        else:
            kept.append(mol)

    print(f"Removed {len(removed)} degenerate entries: {removed}")
    print(f"Keeping {len(kept)} entries")

    # Write back
    writer = Chem.SDWriter(str(sdf_path))
    for mol in kept:
        writer.write(mol)
    writer.close()

    print(f"SDF rewritten: {sdf_path}")


if __name__ == '__main__':
    main()
