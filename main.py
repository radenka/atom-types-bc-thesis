#!/usr/bin/env python3
import argparse
import attyc
import time
from attyc.io import get_available_classifiers
import os, sys
from attyc.exceptions import *
# imports for pdb files


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


if __name__ == "__main__":
    start = time.time()
    args = load_arguments()
    if args:
        try:
            attyc.classify_atoms(args.input_sdf, args.classifier, args.file_output, args.screen_output)
        # what to do with caught exceptions?
        except InputSDFileError as ex:
            print('Problems with input sdfile. End of prog.')
            print(ex.message)
            sys.exit(1)
    end = time.time()
    print('\nExecution time:', end-start)