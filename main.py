#!/usr/bin/env python3
import argparse
import attyc
import sys
from attyc.io import get_available_classifiers
from attyc.exceptions import *


def load_arguments():
    """
    Parses command line arguments while ATTYC being called through this module.
    :return: parsed arguments in arg parser
    """
    parser = argparse.ArgumentParser(prog='ATTYC: ATom TYpe Classification',
                                     description='Atom type assigning based on chemical properties of atoms.'
                                                 ' Used for parametrization of empirical methods for partial atomic'
                                                 ' charges calculation.',
                                     epilog='End of help block. Try it yourself. Good luck!'
                                     )

    parser.add_argument('input_file',
                        type=str,
                        help='File to process.')
    parser.add_argument('classifier',
                        type=str,
                        help="Atomic property that atom types assigned to atoms are derived from. "
                             "Select one of these:" + get_available_classifiers() + "."
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
    """
    Runs classification.
    """
    args = load_arguments()
    if args:
        try:
            attyc.classify_atoms(args.input_file, args.classifier, args.file_output, args.screen_output)
        except (InputFileError, FileScreenOutputError, ClassifierNameError, ClassifierClassError, ExternalTypesInputFileError) as ex:
            print(ex.message)
            sys.exit(1)
