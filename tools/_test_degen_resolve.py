import sys
sys.path.insert(0, "src")
import warnings
warnings.filterwarnings("ignore")
from pyPept.sequence import Sequence
from pyPept.molecule import Molecule
from rdkit.Chem import MolToSmiles

pairs = [
    ("Bn-A-G-Bn", "Bn_-A-G-_Bn"),
    ("Me-A-G-Me", "Me_-A-G-_Me"),
    ("Et-A-G-Et", "Et_-A-G-_Et"),
    ("OMe-A-G-OMe", "OMe_-A-G-_OMe"),
    ("OBn-A-G-OBn", "OBn_-A-G-_OBn"),
    ("NHBn-A-G-NHBn", "NHBn_-A-G-_NHBn"),
]

for base, explicit in pairs:
    s1 = MolToSmiles(Molecule(Sequence(base)).get_molecule(fmt="ROMol"))
    s2 = MolToSmiles(Molecule(Sequence(explicit)).get_molecule(fmt="ROMol"))
    print(f"{base:20s} == {explicit:20s} : {s1 == s2}")
