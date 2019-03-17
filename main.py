from rdkit import Chem
import argparse
import sys
import os
from molecule_set import MoleculeSet
import pprint

parser = argparse.ArgumentParser(prog='ATClass',
                                 description="Atom type assigning based on chemical properties of atoms."
                                             " Used for parametrization of empirical methods for partial atomic"
                                             " charges calculation."
                                             "\nEXIT STATUS:",
                                 epilog='End of help block. Now try it yourself. Good luck!'
                                 )

parser.add_argument("input",
                    type=str,
                    help="SDF file to process.")
parser.add_argument("classificator",
                    type=str,
                    help="Atomic property that atom type assigned to atom is derived from."
                         "\nSelect one of these:"
                         "\n\t'hybrid' "
                         "\n\t'hbo' "
                         "\n\t'group' "
                         "\n\t'partners' "
                         "\n\t'global' "
                    )
args = parser.parse_args()


def check_arguments():
    if not os.path.exists(args.input):
        print(
            "ERROR: Wrong \"input\" argument. Check if the file and path directories exist."
            "\nEnd of program.")
        sys.exit(1)

    if not os.path.isfile(args.input):
        print(
            "ERROR: Wrong \"input\" argument. Regular file not found in given path."
            "\nEnd of program.")
        sys.exit(1)

    if args.classificator not in ['hybrid', 'hbo', 'group', 'partners', 'global']:
        print(
            "ERROR: Wrong \"classificator\" argument. Use one of: 'hybrid', 'hbo', 'group', 'partners', 'global'."
            "\nEnd of program.")
        sys.exit(1)

    return True


if __name__ == "__main__":
    if check_arguments():
        # solve in different way, create Classification class?
        sdf_set = MoleculeSet(args.input)
        atom_types = sdf_set.classify_atoms(args.classificator)
        print_final = pprint.PrettyPrinter(indent=2)
        print_final.pprint(atom_types)

        # just my curiosity:
        molecs = [set(mol) for mol in atom_types]
        ret = molecs[0]
        for i in range(1, len(molecs)):
            ret = ret.union(molecs[i])
        print(sorted(list(ret)))


