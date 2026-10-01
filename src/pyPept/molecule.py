"""
Class to create a rdkit molecule from BILN sequence information.

From publication: pyPept: a python library to generate atomistic 2D and 3D representations of peptides
Journal of Cheminformatics, 2023
"""


__credits__ = ["Rodrigo Ochoa", "J.B. Brown", "Thomas Fox"]
__license__ = "MIT"


# pylint: disable=E1101

import warnings

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdDepictor

from pyPept.attachments import reaction_for_types
from pyPept.leaving_groups import restore_leaving_groups
from pyPept.peptide import Endpoint, Peptide
from pyPept.structure import require_supported_stereo

# SanitizeMol fires "not removing hydrogen atom without neighbors" when it
# reconciles H-counts on atoms adjacent to removed dummy atoms.  The output
# molecule is chemically correct; suppress the noise.
RDLogger.DisableLog('rdApp.warning')


class Molecule:
    """
    Wrapper class around a rdkit ROMol object, with customization for
    peptides.
    """

    def __init__(self, sequence=None, depiction='local'):
        """
        Initialize from a Sequence or an already resolved Peptide.
        :param sequence: parsed source or resolved peptide to assemble
        :type sequence: pyPept.Sequence or pyPept.peptide.Peptide
        :param depiction: method to generate a 2D image
                          The local method is used by default.
                          Use 'rdkit' for RDKit depiction, or None to assemble
                          chemistry without generating coordinates.
        :type depiction: str or None
        """

        self.mol = []
        self.offset = []
        self.bondlist = []
        self.monomers = []
        if depiction not in ('rdkit', 'local', None):
            raise ValueError(
                f"Depiction was {depiction}, expected 'rdkit', 'local', or None.")
        self.depiction = depiction

        self.__from_sequence(sequence)

        if not isinstance(self.mol, Chem.rdchem.Mol):
            raise RuntimeError('pyPept.Molecule initialization failure: ' +
                'problem initializing rdkit.ROMol')

    def __assemble(self, peptide):
        """Join resolved endpoints and restore their unconsumed leaving groups.

        Reaction isotopes are temporary labels allocated per endpoint. They do
        not encode occurrence IDs or slots, which can be large or nonconsecutive.
        """
        from pyPept.interfaces.reaction_library import run_bond_smirks

        if not peptide.occurrences:
            raise ValueError("Cannot assemble an empty peptide")
        labels, chemistry, leaving_groups, pool = {}, {}, {}, {}
        for node in peptide.occurrences:
            if node.definition is None:
                raise ValueError("Assembly requires resolved monomer definitions")
            molecule = node.definition.copy_template()
            require_supported_stereo(molecule)
            metadata = {
                slot: (chem_type, leaving)
                for slot, chem_type, leaving, _ in node.definition.attachment_metadata
            }
            for atom in molecule.GetAtoms():
                atom.SetIntProp('_residue_idx', node.id)
                if atom.GetAtomicNum() != 0:
                    continue
                endpoint = Endpoint(node.id, atom.GetIsotope())
                label = len(labels) + 1
                labels[endpoint] = label
                chemistry[endpoint], leaving_groups[label] = metadata[endpoint.slot]
                atom.SetIsotope(label)
                if atom.GetIsotope() != label:
                    raise ValueError(
                        "Too many attachment sites for RDKit isotope labels"
                    )
            pool[node.id] = molecule

        parent = {identity: identity for identity in pool}

        def find_root(identity):
            while parent[identity] != identity:
                parent[identity] = parent[parent[identity]]
                identity = parent[identity]
            return identity

        for edge in peptide.connections:
            left, right = edge.endpoints
            root1 = find_root(left.occurrence_id)
            root2 = find_root(right.occurrence_id)
            type1, type2 = chemistry[left], chemistry[right]
            reaction = reaction_for_types(type1, type2)
            if reaction is None:
                raise ValueError(
                    f"No reaction defined for chem_type pair ({type1!r}, {type2!r}) "
                    f"between monomer {left.occurrence_id} (slot {left.slot}) "
                    f"and monomer {right.occurrence_id} (slot {right.slot}). "
                    "This attachment chemistry is not supported by the "
                    "current reaction library."
                )
            intramolecular = root1 == root2
            try:
                product = run_bond_smirks(
                    pool[root1], pool[root2], labels[left], labels[right],
                    reaction, intramolecular,
                )
            except ValueError as exc:
                raise ValueError(
                    f"Bond formation failed between monomer {left.occurrence_id} "
                    f"(slot {left.slot}, {type1}) and monomer {right.occurrence_id} "
                    f"(slot {right.slot}, {type2}): {exc}"
                ) from exc
            pool[root1] = product
            if not intramolecular:
                del pool[root2]
                parent[root2] = root1

        fragments = iter(pool.values())
        combined = next(fragments)
        for fragment in fragments:
            combined = Chem.CombineMols(combined, fragment)
        return restore_leaving_groups(combined, leaving_groups, sanitize=False)

    def __fixDihedrals(self):
        """
        Fast fix to get a reasonable 2D representation of the Molecule.

        :return: None
        """

        psi_mol = Chem.MolFromSmiles('NCC(=O)N')
        phi_mol = Chem.MolFromSmiles('C(=O)NCC(=O)')
        amid_mol = Chem.MolFromSmiles('CC(=O)NC')
        sc_mol = Chem.MolFromSmarts('[CH][CX4][CH2][#6]')

        psi_matches = self.mol.GetSubstructMatches(psi_mol)
        phi_matches = self.mol.GetSubstructMatches(phi_mol)
        amid_matches = list(self.mol.GetSubstructMatches(amid_mol))
        sidechain_matches = list(self.mol.GetSubstructMatches(sc_mol))
        sidechain_matches1 = list(self.mol.GetSubstructMatches(
            Chem.MolFromSmarts('N[CX4][CX4][#6]')))
        [sidechain_matches.append(sm1) for sm1 in sidechain_matches1]

        # Fix for Proline
        proline = Chem.MolFromSmiles('CC(=O)N1C(C(=O))CCC1')
        proline_matches = list(self.mol.GetSubstructMatches(proline))
        [amid_matches.append(sm1) for sm1 in proline_matches]

        # Make sure all open-chain angles are 120 deg
        any_angle = Chem.MolFromSmarts('[*]~[*x0]~[*]')
        all_angles = list(self.mol.GetSubstructMatches(any_angle))

        # This is for exocyclic bonds (e.g. phenol OH)
        any_angle1 = Chem.MolFromSmarts('[*]~[*]~[*x0]')
        all_angle1 = list(self.mol.GetSubstructMatches(any_angle1))

        [all_angles.append(sm1) for sm1 in all_angle1]

        bond4 = Chem.MolFromSmarts('[*D4]')
        bond4 = list(self.mol.GetSubstructMatches(bond4))
        new_angles = []
        for b4 in bond4:
            b4 = list(b4)[0]

            for idx in range(len(all_angles) - 1, -1, -1):
                angle = all_angles[idx]

                if angle[1] == b4:
                    all_angles.remove(angle)
                    new_angles.append(angle)

        new_angles = set(new_angles)

        confs = self.mol.GetConformers()
        if len(confs) != 1:
            warnings.warn(f"{len(confs)} conformers for molecule, expected one.")

        count_error = 0
        for conf in confs:
            for pm in psi_matches:
                psi = [pm[0], pm[1], pm[2], pm[4]]
                try:
                    Chem.rdMolTransforms.SetDihedralDeg(conf, *psi, 180.)
                except:
                    count_error+=1

            for pm in phi_matches:
                phi = [pm[0], pm[2], pm[3], pm[4]]
                try:
                    Chem.rdMolTransforms.SetDihedralDeg(conf, *phi, 180.)
                except:
                    count_error+=1

            for pm in amid_matches:
                phi = [pm[0], pm[1], pm[3], pm[4]]
                try:
                    Chem.rdMolTransforms.SetDihedralDeg(conf, *phi, 180.)
                except:
                    count_error+=1

            for pm in sidechain_matches:
                phi = [pm[0], pm[1], pm[2], pm[3]]
                try:
                    Chem.rdMolTransforms.SetDihedralDeg(conf, *phi, 180.)
                except:
                    count_error+=1

            for am in all_angles:
                try:
                    Chem.rdMolTransforms.SetAngleDeg(conf, *am, 120.)
                except:
                    count_error+=1

    def __from_sequence(self, sequence):
        """
        Convert parsed source or a resolved peptide into a rdkit mol object.

        :param sequence: Input sequence
        :type sequence: pyPept.Sequence, pyPept.peptide.Peptide, or None

        :return: rdkit.Chem.ROMol object.
        """

        if sequence is None:
            return None
        if isinstance(sequence, Peptide):
            peptide = sequence
        else:
            if not sequence.is_valid():
                warnings.warn('Invalid sequence object given! Returning None.')
                return None
            # Preserve the legacy inspection attributes for Sequence callers.
            self.monomers = sequence.s_monomers
            self.bondlist = sequence.s_bonds
            peptide = Peptide.from_sequence(sequence)
        self.mol = self.__assemble(peptide)

        # Sanitize the completed product before generating coordinates.
        try:
            Chem.SanitizeMol(self.mol)
        except Exception as exc:
            raise ValueError(
                "Molecule sanitization failed after assembly — the assembled "
                "structure has invalid valence. This usually means incompatible "
                "R-group attachment points were joined. Check R-group assignments "
                f"in the BILN sequence. RDKit detail: {exc}"
            ) from exc

        if self.depiction == 'rdkit':
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(self.mol, clearConfs=True)

        if self.depiction=='local':
            AllChem.Compute2DCoords(self.mol, clearConfs=True)
            Chem.AssignAtomChiralTagsFromStructure(self.mol)
            self.__fixDihedrals()

        return self.mol

    def write_molecule(self, fmt=None, out_file=None):
        """
        Write a molecule to a file or stream in a supported format.

        :param fmt: molecule output format. Allowed formats: SDF, PDB, Smiles, ROMol
                          - SDF: string in the format of a sd-file which can be written to file
                          - PDB: string in the format of a pdb-file
                          - Smiles: smiles representation of sequence
                          - ROMol: rdkit.Chem.ROMol object
        :type fmt: string
        :param out_file: output file to be written
                            if outfile = '-': result is written to stdout.
                            if outfile = None: result is returned to calling context.
        :type out_file: string
        """
        mol = self.get_molecule(fmt=fmt)

        if out_file is None:
            return mol
        if out_file == '-':
            print(mol)
        else:
            with open(out_file, 'w', encoding='utf-8') as out_f:
                out_f.write(mol)

    def get_molecule(self, fmt=None):
        """
        Get a molecule representation of the initialized Sequence object.

        :param fmt: molecule output format. Allowed formats: SDF, PDB, Smiles, ROMol
                          - SDF: string in the format of a sd-file which can be written to file
                          - PDB: string in the format of a pdb-file
                          - Smiles: smiles representation of sequence
                          - ROMol: rdkit.Chem.ROMol object
        :type fmt: string
        """

        if fmt == 'SDF':
            mol = Chem.MolToMolBlock(self.mol, includeStereo=True)
        elif fmt == 'PDB':
            mol = Chem.MolToPDBBlock(self.mol)
        elif fmt == 'Smiles':
            mol = Chem.MolToSmiles(self.mol)
        elif fmt == 'ROMol':
            mol = self.mol
        else:
            warnings.warn('fmt must be one of ' +
                '"SDF", "PDB", "Smiles", "ROMol" but was ' +
                f'{fmt}. Return value will be None.')
            mol = None

        return mol

    def get_residue_atom_map(self):
        """
        Return {occurrence_id: [atom_indices]} for the assembled molecule.

        Sequence inputs use their original zero-based residue indices. Resolved
        Peptide inputs retain their occurrence IDs even when reordered.

        Each atom inherits its original monomer's property at every reaction
        step. Leaving-group replacement preserves it as well. Unmapped atoms
        are omitted rather than guessed from a neighboring occurrence.
        """
        if self.mol is None:
            return {}
        mapping = {}
        for atom in self.mol.GetAtoms():
            if atom.HasProp('_residue_idx'):
                mapping.setdefault(atom.GetIntProp('_residue_idx'), []).append(atom.GetIdx())
        return mapping
