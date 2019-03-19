import argparse
import sys
import os
from classification import Classification
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
                         "\n\t'substruct' "
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

    if args.classificator not in ['hybrid', 'hbo', 'substruct', 'partners', 'global']:
        print(
            "ERROR: Wrong \"classificator\" argument. Use one of: 'hybrid', 'hbo', 'substruct', 'partners', 'global'."
            "\nEnd of program.")
        sys.exit(1)

    return True


if __name__ == "__main__":
    if check_arguments():
        atom_types = Classification(args.input, args.classificator).get_atom_types()
        print_final = pprint.PrettyPrinter(indent=2)
        # print_final.pprint(atom_types)

        # just my curiosity:
        # molecs = [set(mol) for mol in atom_types]
        # ret = molecs[0]
        # for i in range(1, len(molecs)):
        #     ret = ret.union(molecs[i])
        # print(sorted(list(ret)))


