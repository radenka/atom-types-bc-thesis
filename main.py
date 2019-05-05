#!/usr/bin/env python3
import argparse
import attyc
from attyc.io import get_available_classifiers, create_output_file_for_parametrization, count_atom_types
import os, sys
from attyc.exceptions import *

from rdkit import Chem
from collections import Counter
import pprint


def load_arguments():
    parser = argparse.ArgumentParser(prog='ATTYC: ATom TYpe Classification',
                                     description='Atom type assigning based on chemical properties of atoms.'
                                                 ' Used for parametrization of empirical methods for partial atomic'
                                                 ' charges calculation.',
                                     epilog='End of help block. Try it yourself. Good luck!'
                                     )

    parser.add_argument('input_sdf',
                        type=str,
                        help='SDF file to process.')
    parser.add_argument('classifier',
                        type=str,
                        help="Atomic property that atom types assigned to atoms are derived from. "
                             "Select one of these:" + get_available_classifiers()
                        )
    parser.add_argument('--file_output',
                        action='store_true',
                        help='Atom types will be written into text file.')
    parser.add_argument('--screen_output',
                        action='store_true',
                        help='Atom types will be printed on screen.')
    args = parser.parse_args()
    return args


def pdb_trial():
    c_attyp = Counter()
    residues = Counter()
    file = os.path.join(os.path.dirname(__file__), 'pdbs', '2dhc.pdb')
    # removeHs=False doesn't change anything. Depends on pdbfile structure? (no Hs in 2dhc.pdb)
    pdb = Chem.MolFromPDBFile(file)
    pdb2 = Chem.AddHs(pdb, addCoords=True)  # explicit Hs addition works fine, for loop over all atoms doesn't  detect them
    # Chem.MolToPDBFile(pdb2, '2dhc_with_Hs.pdb')  # creates file without Hs
    # print(pdb2.GetNumAtoms())   # returns atoms including Hs
    # res = Chem.SplitMolByPDBResidues(pdb)   # dictionary, doesn't return residues in pdb order but in alphabetical one
    for i in range(pdb2.GetNumAtoms()):
        atom = pdb2.GetAtomWithIdx(i)
        if atom.GetPDBResidueInfo() is None:
            # added H doesn't have any PDBResidueInfo, so that it doesn't have any residue name
            continue
        name = atom.GetPDBResidueInfo().GetName().strip()
        c_attyp[name] += 1
        residues[atom.GetPDBResidueInfo().GetResidueName()] += 1
        if name == 'O2':
            print(atom.GetPDBResidueInfo().GetResidueName())

        # print(atom.GetPDBResidueInfo().GetResidueName())    # works fine
        # name = atom.GetPDBResidueInfo().GetName().strip()
        # print(f'\'{name}\'')
    print_final = pprint.PrettyPrinter(indent=2)
    # print_final.pprint(c_attyp)
    # print(c_attyp)
    print(residues)


if __name__ == "__main__":
    args = load_arguments()
    if args:
        # classifier = None
        try:
            attyc.classify_atoms(args.input_sdf, args.classifier, args.file_output, args.screen_output)
        except (InputSDFileError, FileScreenOutputError, ClassifierNameError, ClassifierClassError, InputSMARTSError) as ex:
            print(ex.message)
            sys.exit(1)