"""PDB residue and atom naming for pyPept Sequence objects.

This module retains pyPept's naming conventions independently of sequence parsing.
"""

__credits__ = ["Rodrigo Ochoa", "J.B. Brown", "Thomas Fox"]
__license__ = "MIT"

import string
from collections import defaultdict

from rdkit import Chem

from pyPept.monomer_store import library_resource, monomer_table


def greekify(mol, name_aa):
    """
    Converts the atom names into relative names using the Greek alphabet.

    :param mol: RDKit molecule object
    :type mol: RDKit molecule
    :param name_aa: name of the AA that will be modified
    :type name_aa: str
    """

    # Some local fixes to name correctly the natural amino acids if required
    fix_aa = {'W': {'CE1': 'CE2', 'NE2': 'NE1', 'CZ1': 'CZ3', 'CH': 'CH2',
                    'CD1': 'CD2', 'CD2': 'CD1'},
              'N': {'OD2': 'OD1', 'ND1': 'ND2'},
              'I': {'CD': 'CD1', 'CG1': 'CG2', 'CG2': 'CG1'},
              'T': {'OG2': 'OG1', 'CG1': 'CG2'},
              'P': {'CG2': 'CG', 'CG1': 'CD'}}

    greek = list('ABGDEZHTIKLMNXOPRS')
    greekdex = defaultdict(list)
    ca_atom = get_atom_by_name(mol, 'CA')

    # Recognize if the atom is not part of the backbone
    for atom in mol.GetAtoms():
        is_backbone = (atom.GetPDBResidueInfo() is not None and
                       atom.GetPDBResidueInfo().GetName().strip() in (
                           'LOWER', 'UPPER', 'N', 'CA', 'C', 'H', 'HA', 'O',
                           'OXT'))
        if atom.GetSymbol() != 'H' and atom.GetSymbol()[
            0] != 'R' and not is_backbone:
            n_atom = len(
                Chem.GetShortestPath(mol, ca_atom.GetIdx(), atom.GetIdx())) - 1
            greekdex[n_atom].append(atom)

    # Iterate over the list of atoms and the greek letters
    for k in greekdex:
        if len(greek) <= k:
            pass
        elif len(greekdex[k]) == 0:
            pass
        elif len(greekdex[k]) == 1:
            # Special case
            new_name = f'{greekdex[k][0].GetSymbol()}{greek[k]}'
            if name_aa in fix_aa:
                if new_name in fix_aa[name_aa]:
                    digit = fix_aa[name_aa][new_name][-1]
                    name = f'{greekdex[k][0].GetSymbol(): >2}{greek[k]}{digit}'
                else:
                    name = f'{greekdex[k][0].GetSymbol(): >2}{greek[k]} '
            else:
                name = f'{greekdex[k][0].GetSymbol(): >2}{greek[k]} '

            greekdex[k][0].GetPDBResidueInfo().SetName(name)

        elif len(greekdex[k]) < 36:
            list_atom = list(string.digits + string.ascii_uppercase)[1:]
            for i, atom in enumerate(greekdex[k]):

                new_name = f'{atom.GetSymbol()}{greek[k]}{list_atom[i]}'
                if name_aa in fix_aa:
                    if new_name in fix_aa[name_aa]:
                        if name_aa != 'P':
                            digit = fix_aa[name_aa][new_name][-1]
                            name = f'{atom.GetSymbol(): >2}{greek[k]}{digit}'
                        else:
                            digit = fix_aa[name_aa][new_name][-1]
                            name = f'{atom.GetSymbol(): >2}{digit} '
                    else:
                        name = f'{atom.GetSymbol(): >2}{greek[k]}{list_atom[i]}'
                else:
                    name = f'{atom.GetSymbol(): >2}{greek[k]}{list_atom[i]}'

                greekdex[k][i].GetPDBResidueInfo().SetName(name)

        else:
            pass


def get_atom_by_name(mol, name):
    """
    Get an atom RDKit object based on a particular name

    :param mol: RDKit object with the molecule
    :type mol: RDKit molecule
    :param name: name of the atom
    :type name: str

    :return: the atom RDKit object
    """

    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info and info.GetName().strip() == name.strip():
            selected_atom = atom

    return selected_atom


def get_monomer_codes(df_name):
    """
    Get the PDB codes from the monomer dictionary
    :param df_name: monomer dictionary
    :type df_name: dataframe

    :return: dictionary with the monomer PDB codes
    """
    # Get the values from the monomer dictionary
    monomers = {}
    for idx in df_name.index:
        symbol = df_name.at[idx, 'm_abbr']
        pdb_name = df_name.at[idx, 'pdbName']
        monomers[symbol] = pdb_name

    return monomers


def correct_pdb_atoms(seq, path="pyPept.data",
                      monomer_lib="monomers.sdf"):
    """
    Pipeline to assign the correct atom names to the pyPept object

    :param seq: pyPept Sequence object
    :type seq: pyPept.sequence
    :return: the modified pyPept Sequence with correct atom names
    """

    # Special case two main N- and C- terminal caps
    names_cap = {'ac': ['CH3', 'C', 'O'], 'am': ['N']}

    # Read the monomer dataframe
    monomer_df_filepath = library_resource(path, monomer_lib)
    new_df = monomer_table(str(monomer_df_filepath))
    if 'pdbName' not in new_df.columns:
        raise ValueError(
            f"Monomer library {monomer_df_filepath} lacks the 'pdbName' metadata "
            "required for PDB atom naming. Supply a library with PDB residue "
            "names, or use --noconf for 2D output."
        )

    # Get monomer codes
    monomers = get_monomer_codes(new_df)

    # Iterate over the monomers
    mm_list = seq.s_monomers
    for i, monomer in enumerate(mm_list):
        mol = monomer['m_romol']
        name = monomer['m_abbr']
        backbone_smiles = Chem.MolFromSmiles('NCC(=O)')
        if i == len(mm_list) - 1:
            backbone_smiles = Chem.MolFromSmiles('NCC(=O)O')
        backbone_atoms = mol.GetSubstructMatches(backbone_smiles)
        aa_flag = 0
        if backbone_atoms:
            bb_index = list(backbone_atoms[0])
            type_mon = new_df.loc[new_df['m_abbr'] == name, 'm_type'].item()
            if type_mon == 'aa':
                aa_flag = 1

        # Iterate over the atoms
        counter = 0
        counter_non = 0
        for j, atom in enumerate(mol.GetAtoms()):
            if aa_flag == 1:
                pos = -1
                if atom.GetIdx() in bb_index:
                    pos = bb_index.index(atom.GetIdx())
                if pos == 0:
                    atomname = ' N  '
                elif pos == 1:
                    atomname = ' CA '
                elif pos == 2:
                    atomname = ' C  '
                elif pos == 3:
                    atomname = ' O  '
                elif pos == 4:
                    atomname = ' OXT'
                else:
                    counter += 1
                    atomname = f' {atom.GetSymbol()}{counter} '
            else:
                # Special case for main capping groups
                if name in ('ac', 'am'):
                    if atom.GetSymbol()[0] != 'R':
                        if names_cap[name][j] == 'CH3':
                            atomname = f' {names_cap[name][j]}'
                        else:
                            atomname = f' {names_cap[name][j]}  '
                else:
                    counter_non += 1
                    if counter_non < 10:
                        atomname = f' {atom.GetSymbol()}{counter_non} '
                    else:
                        atomname = f' {atom.GetSymbol()}{counter_non}'

            # Assign the atom object to the peptide molecule
            info = atom.GetPDBResidueInfo()
            if info is None:
                atom.SetMonomerInfo(Chem.AtomPDBResidueInfo(atomName=atomname,
                                                            serialNumber=atom.GetIdx(),
                                                            residueName=f'{monomers[name]}',
                                                            residueNumber=i + 1,
                                                            chainId="A"))

        # Rename atoms using the greek nomenclature
        if aa_flag == 1:
            greekify(mol, name)

    return seq
