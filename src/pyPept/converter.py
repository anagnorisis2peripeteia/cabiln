"""
Class to to manipulate BILN and convert to HELM (and vice versa)

From publication: pyPept: a python library to generate atomistic 2D and 3D representations of peptides
Journal of Cheminformatics, 2023
"""


__credits__ = ["Rodrigo Ochoa", "J.B. Brown", "Thomas Fox"]
__license__ = "MIT"
__version__ = "1.0"

# pylint: disable=E1101

import copy
import re
import warnings

from pyPept.sequence import SequenceConstants
from pyPept.notation import split_outside
from pyPept.notation import legacy_attachment_slot

class Converter:
    """
    Class to convert between BILN and HELM (and viceversa)
    """
    def __init__(self, biln=None, helm=None, chuckles=None):
        """Construct a BILN, HELM, or CHUCKLES molecule.

        :param biln: BILN sequence string
        :type biln: str
        :param helm: HELM sequence string
        :type helm: str
        :param chuckles: CHUCKLES sequence string (dot-separated monomers,
            pipe-separated chains, same cyclisation syntax as BILN)
        :type chuckles: str
        """

        self.polymerinfo = {"chains": [], "bonds": []}

        given = sum(x is not None for x in (biln, helm, chuckles))
        if given > 1:
            raise ValueError(
                "Converter accepts only one of biln, helm, or chuckles.")
        elif biln is not None:
            self.eval_biln(biln=biln)
        elif helm is not None:
            self.eval_helm(helm=helm)
        elif chuckles is not None:
            self.eval_chuckles(chuckles=chuckles)

    @staticmethod
    def __remove_brackets(sequence):
        """Remove brackets around individual residues

        :param sequence: the peptide sequence
        :type sequence: list

        :return newseq: new list without the brackets
        :type newseq: list
        """
        newseq = []
        for residue in sequence:
            pat = re.sub(r'\[(.*)\]', r'\1', residue)
            newseq.append(pat)

        return newseq

    @staticmethod
    def __split_helm(helm):
        """split HELM string into parts

        :param helm: the helm peptide sequence
        :type sequence: string

        :return helm_parts: a list containing relevant parts from helm
        :type helm_parts: list
        """

        helm_parts = []
        while len(helm):
            match = re.search(r'\$[^,;]', helm)  # this excludes any $ signs in the cxsmiles
            if match is None:
                match = re.search(r'\$', helm)

            if match is None:
                helm_parts.append(helm)
                break

            helm_parts.append(helm[:match.span()[0]])
            helm = helm[match.span()[0] + 1:]

        if len(helm_parts) == 4:
            helm_parts.append('')

        list_of_simple_polymers = helm_parts[0]
        if SequenceConstants.helm_polymer in list_of_simple_polymers:
            list_of_simple_polymers = list_of_simple_polymers.split(
                SequenceConstants.helm_polymer)
        else:
            list_of_simple_polymers = [list_of_simple_polymers]
        helm_parts[0] = list_of_simple_polymers

        list_of_connections = helm_parts[1]
        if list_of_connections != '':
            if SequenceConstants.helm_polymer in list_of_connections:
                list_of_connections = list_of_connections.split(
                    SequenceConstants.helm_polymer)
            else:
                list_of_connections = [list_of_connections]
        else:
            list_of_connections = []

        helm_parts[1] = list_of_connections

        return helm_parts

    def __to_biln(self):
        """Generate BILN from polymerInfo

        :return biln: the generated biln from polymer dictionary
        :type biln: string
        """

        chains = copy.deepcopy(self.polymerinfo["chains"])
        bonds = copy.deepcopy(self.polymerinfo["bonds"])

        # move bond info into monomers using CABILN .!n(y,z) inline notation
        for ibond, bond in enumerate(bonds):
            c1_value, res1, rgroup1, c2_value, res2, rgroup2 = bond
            bid = ibond + 1
            slot1 = legacy_attachment_slot(rgroup1)
            slot2 = legacy_attachment_slot(rgroup2)
            chains[c1_value][res1] = (
                f"{chains[c1_value][res1]}.!{bid}({slot1},{slot2})")
            chains[c2_value][res2] = (
                f"{chains[c2_value][res2]}.!{bid}")

        list_ofsimple_polymers = []
        for chain in chains:
            poly = SequenceConstants.monomer_join.join(chain)
            list_ofsimple_polymers.append(poly)

        biln = SequenceConstants.chain_separator.join(list_ofsimple_polymers)

        if len(biln) == 0:
            biln = None

        return biln

    def __to_helm(self):
        """Generate a HELM from an internal PolymerInfo dictionary.

        :return helm: the generated helm v2.0 from polymer dictionary
        :type helm: string
        """

        chains = copy.deepcopy(self.polymerinfo["chains"])
        bonds = copy.deepcopy(self.polymerinfo["bonds"])

        for count_chain, chain in enumerate(chains):
            for count_res, res in enumerate(chain):
                if len(res) > 1:
                    chains[count_chain][count_res] = f"[{res}]"

        list_of_simple_polymers = []
        for count_chain, chain in enumerate(chains):
            poly = ".".join(chain)
            list_of_simple_polymers.append(f'PEPTIDE{count_chain + 1}{{{poly}}}')

        list_of_connections = []
        for bond in bonds:
            c1_val, r1_val, g1_val, c2_val, r2_val, g2_val = bond
            bond_info = f'PEPTIDE{c1_val + 1},PEPTIDE{c2_val + 1},{r1_val + 1}:R{g1_val}-{r2_val + 1}:R{g2_val}'
            list_of_connections.append(bond_info)

        list_of_connections = SequenceConstants.helm_polymer.join(list_of_connections)
        list_of_simple_polymers = SequenceConstants.helm_polymer.join(list_of_simple_polymers)

        if not list_of_simple_polymers:
            helm = None
        else:
            helm = f"{list_of_simple_polymers}${list_of_connections}$$$V2.0"

        return helm

    def eval_helm(self, helm):
        """Generate internal PolymerInfo dictionary from a HELM string.

        This method fills the polymerinfo dictionary:
        polymerinfo = {"chains": [chain1, chain2, ...],
                       "bonds": [b1, b2, b3]}
        where chains contains the monomer abbreviation for each chain
        and bonds contains the "extra" bonds explicitly given in HELM or BILN as a list

        :param helm: the peptide helm string
        :type helm: string
        """

        if not isinstance(helm, str) or not helm.strip():
            raise ValueError("HELM must be a non-empty string.")
        try:
            list_of_simple_polymers, list_of_connections, groups, annotations,\
            version = self.__split_helm(helm)

        except (ValueError, IndexError) as exc:
            raise ValueError("Invalid HELM: expected five sections separated by '$'.") from exc

        if not list_of_simple_polymers:
            raise ValueError("HELM contains no simple polymers.")


        pattern = re.compile(r'{.*}')

        id_val = []
        polymer = []

        for idx, chain in enumerate(list_of_simple_polymers):

            chain = chain.strip()
            match = pattern.search(chain)
            if match is None:
                raise ValueError(f"HELM polymer contains no sequence: {chain!r}")

            seq = match.span()
            id_chain = chain[:seq[0]]
            if not re.fullmatch(r'PEPTIDE[1-9]\d*', id_chain):
                raise ValueError(f"Unsupported HELM polymer identifier: {id_chain!r}")
            id_chain = int(re.sub('PEPTIDE', '', id_chain))

            poly = chain[seq[0] + 1:seq[1] - 1]

            if not poly:
                raise ValueError(f"HELM polymer PEPTIDE{id_chain} contains no residues.")

            poly = split_outside(poly, SequenceConstants.chain_separator, '[]')

            poly = self.__remove_brackets(poly)

            if id_chain in id_val:
                raise ValueError(f"Duplicate HELM polymer identifier: PEPTIDE{id_chain}")
            id_val.append(id_chain)
            polymer.append(poly)


        # bond information is of the type: 'PEPTIDE1,PEPTIDE2,1:R1-4:R3'
        # split into individual parts 'id1, id2, res1:rgroup1-res2:rgroup2'
        bonds = []
        if len(list_of_connections) > 0:
            for idx, conn in enumerate(list_of_connections):
                id1, id2, bond = conn.split(',')

                res1, rgroup1, res2, rgroup2 = re.split(r'[-:]', bond)

                id1 = int(id1.replace('PEPTIDE', ''))
                id2 = int(id2.replace('PEPTIDE', ''))

                # translate back to current order of chains in polymerinfo["chains"]
                id1 = id_val.index(id1)
                id2 = id_val.index(id2)

                # start counting of residues at 0
                res1 = int(res1) - 1
                res2 = int(res2) - 1

                if not (0 <= res1 < len(polymer[id1]) and
                        0 <= res2 < len(polymer[id2])):
                    raise ValueError("HELM connection refers to an out-of-range residue.")

                # get the number of the Rgroup (keep numbering 1-3)
                rgroup1 = int(rgroup1.replace('R', ''))
                rgroup2 = int(rgroup2.replace('R', ''))

                bonds.append([id1, res1, rgroup1, id2, res2, rgroup2])

        self.polymerinfo = {"chains": polymer, "bonds": bonds}

    def eval_chuckles(self, chuckles):
        """Generate internal PolymerInfo dictionary from a CHUCKLES sequence.

        CHUCKLES format (Siani et al. JCICS 1994):
          - monomers separated by '.' within a chain
          - chains separated by '|'
          - cyclisation bonds via same ``(bondID,RgroupNum)`` syntax as BILN

        Example::

            A.G.K                      # linear tripeptide
            A.G.K(1,3).E(1,3)          # sidechain-cyclised
            A.G.K|D.E                  # two chains

        :param chuckles: CHUCKLES sequence string
        :type chuckles: str
        """
        self._read_legacy(chuckles, "|", ".", "CHUCKLES")

    def __to_chuckles(self):
        """Generate a CHUCKLES sequence string from polymerinfo.

        :return: CHUCKLES string, or None if empty
        :rtype: str or None
        """
        chains = copy.deepcopy(self.polymerinfo["chains"])
        bonds = copy.deepcopy(self.polymerinfo["bonds"])

        for ibond, bond in enumerate(bonds):
            c1_value, res1, rgroup1, c2_value, res2, rgroup2 = bond
            chains[c1_value][res1] = (
                f"{chains[c1_value][res1]}({ibond + 1},{rgroup1})")
            chains[c2_value][res2] = (
                f"{chains[c2_value][res2]}({ibond + 1},{rgroup2})")

        chain_strings = ['.'.join(chain) for chain in chains]
        chuckles = '|'.join(chain_strings)
        return chuckles if chuckles else None

    def eval_biln(self, biln):
        """Generate internal PolymerInfo dictionary from a BILN string.

        This method fills the PolymerInfo dictionary:
        polymerinfo = {"chains": [chain1, chain2, ...],
                       "bonds": [b1, b2, b3]}
        where chains contains the monomer abbreviation for each chain
        and bonds contains the "extra" bonds explicitly given in HELM or 
        BILN as a list.

        :param biln: the peptide BILN string
        :type biln: string

        :return: None
        """

        if not isinstance(biln, str) or not biln.strip():
            raise ValueError("BILN must be a non-empty string.")
        if '%' in biln or re.search(r'\.[!\[{]', biln):
            raise ValueError(
                "Converter(biln=...) accepts legacy BILN, not CABILN attachments. "
                "Use Sequence for CABILN parsing.")
        self._read_legacy(biln, ".", SequenceConstants.monomer_join, "BILN")

    def _read_legacy(self, text, chain_separator, monomer_separator, kind):
        """Read paired legacy endpoints once for BILN and CHUCKLES."""
        chains, endpoints = [], {}
        for chain_index, chain in enumerate(text.split(chain_separator)):
            residues = []
            for residue_index, residue in enumerate(chain.split(monomer_separator)):
                for bond_id, slot in re.findall(r"\((\d+),(\d+)\)", residue):
                    endpoints.setdefault(bond_id, []).extend(
                        (chain_index, residue_index, int(slot))
                    )
                residues.append(re.sub(r"\(.*\)", "", residue))
            chains.append(residues)
        self.polymerinfo["chains"] = chains
        for bond_id, bond in endpoints.items():
            if len(bond) != 6:
                raise ValueError(f"{kind} bond {bond_id} must have exactly two endpoints.")
        self.polymerinfo["bonds"] = list(endpoints.values())

    def get_helm(self):
        """Return a HELM string from the PolymerInfo generated at
        instantiation.
        
        :return: str"""
        helm = self.__to_helm()
        return helm

    def get_biln(self):
        """Return CABILN text under the retained historical method name.

        Legacy input and HELM attachment R3 become modern CABILN R4. Pass this
        result to Sequence, not Converter(biln=...), which reads legacy BILN.

        :return: str or None"""
        biln = self.__to_biln()
        return biln

    def get_chuckles(self):
        """Return a CHUCKLES sequence string from the PolymerInfo generated
        at instantiation.

        :return: str or None
        """
        return self.__to_chuckles()
