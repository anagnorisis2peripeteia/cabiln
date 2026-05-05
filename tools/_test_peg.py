import sys
sys.path.insert(0, "src")
from pyPept.sequence import Sequence
from pyPept.molecule import Molecule
from rdkit.Chem import MolToSmiles
import warnings
warnings.filterwarnings("ignore")

try:
    seq = Sequence("ac-K.[PEG(4,1).am(2,1)]-G-am")
    mol = Molecule(seq)
    romol = mol.get_molecule(fmt="ROMol")
    print("bracket OK:", MolToSmiles(romol) if romol else "NO MOL")
except Exception as e:
    print("bracket FAIL:", e)

try:
    seq2 = Sequence("ac-K.!1(4,1)-G-am%!1-PEG-am")
    mol2 = Molecule(seq2)
    romol2 = mol2.get_molecule(fmt="ROMol")
    print("branch OK:", MolToSmiles(romol2) if romol2 else "NO MOL")
except Exception as e:
    print("branch FAIL:", e)
