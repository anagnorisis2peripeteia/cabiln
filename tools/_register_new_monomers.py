"""Design, validate, and register new monomers for README walkthrough coverage.

Each type with <5 library monomers gets enough new ones to reach 5.
Run with --dry-run to validate without registering.
"""
import sys
import io
import argparse

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

from pyPept.interfaces.monomer_pipeline import pre_activate
from pyPept.interfaces.cli_monomer import register_monomer

NEW_MONOMERS = [
    # selenol (have Sec, seC → need 3)
    ("Se_Hcy", "N[C@@H](CC[SeH])C(=O)O"),
    ("D_Se_Hcy", "N[C@H](CC[SeH])C(=O)O"),
    ("Se_Pen", "NC(C)(C[SeH])C(=O)O"),
    # cyclooctyne_c (have CyoAla, CyoLys → need 3)
    ("CyoOrn", "N[C@@H](CCCC1CCCCC#CC1)C(=O)O"),
    ("CyoDab", "N[C@@H](CCC1CCCCC#CC1)C(=O)O"),
    ("Cyo", "OC(=O)CC1CCCCC#CC1"),
    # maleimide_c (have MalAla, MalLys, Mal → need 2)
    ("MalOrn", "N[C@@H](CCCN1C(=O)C=CC1=O)C(=O)O"),
    ("MalDab", "N[C@@H](CCN1C(=O)C=CC1=O)C(=O)O"),
    # tetrazine_c (have TzAla, TzLys, Tz → need 2)
    ("TzOrn", "N[C@@H](CCCCc1nncnn1)C(=O)O"),
    ("TzDab", "N[C@@H](CCCc1nncnn1)C(=O)O"),
    # tco_c (have TcoAla, TcoLys, TCO → need 2)
    ("TcoOrn", "N[C@@H](CCCC1=CCCCCCC1)C(=O)O"),
    ("TcoDab", "N[C@@H](CCC1=CCCCCCC1)C(=O)O"),
    # aminooxy (have Aoa → need 4)
    ("AoaLys", "N[C@@H](CCCCON)C(=O)O"),
    ("AoaOrn", "N[C@@H](CCCON)C(=O)O"),
    ("AoaDab", "N[C@@H](CCON)C(=O)O"),
    ("D_Aoa", "N[C@H](CON)C(=O)O"),
    # formamide_c: REMOVED — amide_nh protects the C atom via _SC_PATTERNS,
    # so formamide_c is never the winning pre_activate type. It only appears
    # at assembly time via infer_chem_type.
    # nhs_ester (have NHSAla → need 4)
    ("NHSLys", "N[C@@H](CCCCC(=O)ON1C(=O)CCC1=O)C(=O)O"),
    ("NHSOrn", "N[C@@H](CCCC(=O)ON1C(=O)CCC1=O)C(=O)O"),
    ("NHSDab", "N[C@@H](CCC(=O)ON1C(=O)CCC1=O)C(=O)O"),
    ("D_NHSAla", "N[C@H](CC(=O)ON1C(=O)CCC1=O)C(=O)O"),
    # hydrazide (have HzAla → need 4)
    ("HzLys", "N[C@@H](CCCCC(=O)NN)C(=O)O"),
    ("HzOrn", "N[C@@H](CCCC(=O)NN)C(=O)O"),
    ("HzDab", "N[C@@H](CCC(=O)NN)C(=O)O"),
    ("D_HzAla", "N[C@H](CC(=O)NN)C(=O)O"),
]


def validate_all():
    ok = 0
    fail = 0
    for symbol, smiles in NEW_MONOMERS:
        try:
            result = pre_activate(smiles)
            types = [ct for ct in result.chem_types.values() if ct]
            sc_types = {
                k: v for k, v in result.chem_types.items() if k >= 4 and v
            }
            print(f"OK  {symbol:15s} {smiles}")
            print(f"    CHUCKLES: {result.chuckles}")
            print(f"    sidechain types: {sc_types}")
            ok += 1
        except Exception as e:
            print(f"FAIL {symbol:15s} {smiles}")
            print(f"     {e}")
            fail += 1
    print(f"\n{ok} OK, {fail} FAIL out of {len(NEW_MONOMERS)}")
    return fail == 0


def register_all():
    for symbol, smiles in NEW_MONOMERS:
        try:
            register_monomer(smiles, symbol)
            print(f"Registered {symbol}")
        except Exception as e:
            print(f"FAILED {symbol}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate only, don't register",
    )
    args = parser.parse_args()

    if args.dry_run:
        validate_all()
    else:
        if validate_all():
            print("\nAll validated. Registering...")
            register_all()
        else:
            print("\nFix failures before registering.")
            sys.exit(1)
