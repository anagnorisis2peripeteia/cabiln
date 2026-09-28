"""
Class to create a rdkit molecule from BILN sequence information.

From publication: pyPept: a python library to generate atomistic 2D and 3D representations of peptides
Journal of Cheminformatics, 2023
"""

########################################################################################
# Authorship
########################################################################################

__credits__ = ["Rodrigo Ochoa", "J.B. Brown", "Thomas Fox"]
__license__ = "MIT"


########################################################################################
# Modules
########################################################################################
# pylint: disable=E1101

# System libraries
import warnings

# Third-party libraries
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdDepictor

from pyPept.leaving_groups import restore_leaving_groups
from pyPept.structure import require_supported_stereo

# SanitizeMol fires "not removing hydrogen atom without neighbors" when it
# reconciles H-counts on atoms adjacent to removed dummy atoms.  The output
# molecule is chemically correct; suppress the noise.
RDLogger.DisableLog('rdApp.warning')

##########################################################################
# Functions and classes
##########################################################################

class Molecule:
    """
    Wrapper class around a rdkit ROMol object, with customization for
    peptides.
    """

    ############################################################################
    def __init__(self, sequence=None, depiction='local'):
        """
        Initialize a Molecule object, optionally with a Sequence object.
        :param sequence:  input sequence to be converted to a molecule
        :type sequence: pyPept.Sequence object.
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

        # Main function, will modify data members defined above.
        self.__from_sequence(sequence)

        if not isinstance(self.mol, Chem.rdchem.Mol):
            raise RuntimeError('pyPept.Molecule initialization failure: ' +
                'problem initializing rdkit.ROMol')

    ############################################################################
    def __add_bonds_to_mol(self, sequence):
        """
        Form inter-monomer bonds using SMIRKS reactions from reactions.yaml.

        Each dummy atom in every monomer is relabeled to a globally unique
        isotope  (m_idx + 1) * 100 + original_isotope  before any reactions
        run.  This means that even when multiple monomers of the same type
        (e.g. two Cys) are present, each dummy is uniquely addressable and
        the targeted SMIRKS can match exactly the right pair.

        Supported connections are defined by the reaction library.
        """
        from pyPept.attachments import resolve_connection
        from pyPept.interfaces.reaction_library import run_bond_smirks
        from pyPept.sequence import _slot_for_attachment

        monomers_orig = [mon['m_romol'] for mon in sequence.s_monomers]
        for monomer in monomers_orig:
            require_supported_stereo(monomer)

        # ── Step 1: relabel all dummies to globally unique isotopes ──────────
        # unique_iso(m_idx, slot) = (m_idx + 1) * 100 + slot  (slot is 1-based)
        tagged = []
        for m_idx, mol in enumerate(monomers_orig):
            rw = Chem.RWMol(mol)
            for atom in rw.GetAtoms():
                atom.SetIntProp('_residue_idx', m_idx)
                if atom.GetAtomicNum() == 0:
                    orig = atom.GetIsotope()      # 1-based slot index
                    atom.SetIsotope((m_idx + 1) * 100 + orig)
            tagged.append(rw.GetMol())

        # ── Step 2: fragment pool with union-find ─────────────────────────────
        pool = {i: tagged[i] for i in range(len(tagged))}
        parent = list(range(len(tagged)))

        def find_root(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        # ── Step 3: process each bond ─────────────────────────────────────────
        for bond in self.bondlist:
            m1, at1, m2, at2 = bond[:4]
            root1, root2 = find_root(m1), find_root(m2)
            frag1, frag2 = pool[root1], pool[root2]

            # Globally unique isotopes for this bond's two attachment dummies.
            # Use the stored slot when present (avoids _slot_for_attachment
            # ambiguity when an atom neighbours multiple dummies, e.g. backbone N).
            slot1 = bond[4] if len(bond) > 4 else _slot_for_attachment(monomers_orig[m1], at1)
            slot2 = bond[5] if len(bond) > 4 else _slot_for_attachment(monomers_orig[m2], at2)
            iso1 = (m1 + 1) * 100 + slot1
            iso2 = (m2 + 1) * 100 + slot2

            # Look up SMIRKS reaction
            ct1, ct2, entry = resolve_connection(
                monomers_orig[m1], slot1, monomers_orig[m2], slot2,
                leaving_groups1=sequence.s_monomers[m1].get('m_Rgroups'),
                leaving_groups2=sequence.s_monomers[m2].get('m_Rgroups'),
            )

            if entry is None:
                raise ValueError(
                    f"No reaction defined for chem_type pair ({ct1!r}, {ct2!r}) "
                    f"between monomer {m1} (slot {slot1}) and monomer {m2} "
                    f"(slot {slot2}). Check the slot indices in your CABILN — "
                    f"this is most likely a wrong R-group index."
                )
            intramol = (root1 == root2)

            try:
                product = run_bond_smirks(frag1, frag2, iso1, iso2, entry, intramol)
            except ValueError as exc:
                raise ValueError(
                    f"Bond formation failed between monomer {m1} (slot {slot1}, "
                    f"{ct1}) and monomer {m2} (slot {slot2}, {ct2}): {exc}"
                ) from exc

            pool[root1] = product
            if not intramol:
                del pool[root2]
                parent[root2] = root1

        # ── Step 4: combine any remaining disconnected fragments ──────────────
        frags = list(pool.values())
        combined = frags[0]
        for frag in frags[1:]:
            combined = Chem.CombineMols(combined, frag)

        self.mol = Chem.RWMol(combined)

    ############################################################################
    def __restore_and_remove_rgroups(self, sequence):
        """
        Restore leaving-group atoms on unbound R-group slots, then remove all
        remaining dummy atoms.

        After SMIRKS-based assembly, bonded dummies are already consumed by the
        reactions.  Only unbound dummies remain in self.mol.  Each is identified
        by its globally unique isotope  (m_idx + 1) * 100 + slot  (slot 1-based),
        which maps back to the monomer's m_Rgroups list to determine whether
        to remove the dummy ([H]) or replace it with a leaving-group atom ([OH]).
        """
        # Build leaving-group lookup: globally-unique isotope → leaving SMILES
        leaving_for_iso: dict = {}
        for m_idx, mon in enumerate(sequence.s_monomers):
            for slot_idx, lg in enumerate(mon['m_Rgroups']):
                iso = (m_idx + 1) * 100 + (slot_idx + 1)
                leaving_for_iso[iso] = lg

        return restore_leaving_groups(self.mol, leaving_for_iso, sanitize=False)

    ########################################################################################
    def __fixDihedrals(self):
        """
        Fast fix to get a reasonable 2D representation of the Molecule.

        :return: None
        """

        # Define a number of substructures for phi, psi, amide group and sidechains
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

        # Search for atoms that have 4 connection points, then remove these
        #   from the all_angles list
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

        # Get conformers
        confs = self.mol.GetConformers()
        if len(confs) != 1:
            warnings.warn(f"{len(confs)} conformers for molecule, expected one.")

        # Update the dihedrals
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
    # end of Molecule.__fixDihedrals()

    ########################################################################################
    def __from_sequence(self, sequence):
        """
        Function to convert a pyPept.Sequence object into a rdkit mol object.

        :param sequence: Input sequence
        :type sequence: pyPept.Sequence or None

        :return: rdkit.Chem.ROMol object.
        """

        if sequence is None:
            return None
        if not sequence.is_valid():
            warnings.warn('Invalid sequence object given! Returning None.')
            return None

        self.monomers = sequence.s_monomers
        self.bondlist = sequence.s_bonds

        # Step 1: form all inter-monomer bonds via SMIRKS, building self.mol
        # as the combined product.  Uses globally-unique dummy isotopes so that
        # each reactive site is unambiguously targeted even when multiple
        # monomers of the same type are present.
        self.__add_bonds_to_mol(sequence)

        # Step 2: restore unbound R-group slots (remove [H] dummies, replace
        # others with their leaving-group atom).
        self.mol = self.__restore_and_remove_rgroups(sequence)

        # Step 4: sanitize and generate 2D coords
        try:
            Chem.SanitizeMol(self.mol)
        except Exception as exc:
            raise ValueError(
                "Molecule sanitization failed after assembly — the assembled "
                "structure has invalid valence. This usually means incompatible "
                "R-group attachment points were joined. Check R-group assignments "
                f"in the BILN sequence. RDKit detail: {exc}"
            ) from exc

        # Compute 2D coordinates
        if self.depiction == 'rdkit':
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(self.mol, clearConfs=True)

        if self.depiction=='local':
            AllChem.Compute2DCoords(self.mol, clearConfs=True)
            Chem.AssignAtomChiralTagsFromStructure(self.mol)
            self.__fixDihedrals()

        return self.mol

    ########################################################################################
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

    ########################################################################################
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

    ########################################################################################
    def get_residue_atom_map(self):
        """
        Return {residue_idx: [atom_indices]} for the assembled molecule.

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

    # End of the Molecule class declaration.

############################################################
# End of molecule.py
############################################################
