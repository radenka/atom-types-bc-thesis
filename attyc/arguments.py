import argparse
import os
import sys


def load_arguments():
    parser = argparse.ArgumentParser(prog='ATTYC: ATom TYpe Classification',
                                     description='Atom type assigning based on chemical properties of atoms.'
                                                 ' Used for parametrization of empirical methods for partial atomic'
                                                 ' charges calculation.'
                                                 '\nEXIT STATUS:',
                                     epilog='End of help block. Now try it yourself. Good luck!'
                                     )

    parser.add_argument('input_sdf',
                        type=str,
                        help='SDF file to process.')
    parser.add_argument('classifier',
                        type=str,
                        help="Atomic property that atom types assigned to atoms are derived from."
                             "\nSelect one of these:"
                             "\n\t'hybrid' "
                             "\n\t'hbo' "
                             "\n\t'substruct' "
                             "\n\t'partners' "
                             "\n\t'global' "
                        )
    parser.add_argument('--file_output',
                        action='store_true',
                        help='Atom types output will be written to text file.')
    parser.add_argument('--verbose',
                        action='store_true',
                        help='Atom types output will be printed on screen.')
    args = parser.parse_args()
    return args


def check_arguments(input_sdf, classifier, file_output, screen_output):
    if not os.path.exists(input_sdf):
        print(
            "ERROR: Wrong \"input\" argument. Check if the file and path directories exist."
            "\nEnd of program.")
        sys.exit(1)

    if not os.path.isfile(input_sdf):
        print(
            "ERROR: Wrong \"input\" argument. Regular file not found in given path."
            "\nEnd of program.")
        sys.exit(1)

    if type(file_output) != bool or type(screen_output) != bool:
        print(
            "ERROR: Wrong \"output\" argument. Use boolean only."
            "\nEnd of program.")
        sys.exit(1)

    if classifier not in ['hybrid', 'hbo', 'substruct', 'partners', 'global']:
        print(
            "ERROR: Wrong \"classifier\" argument. Use one of these: 'hybrid', 'hbo', 'substruct', 'partners', 'global'."
            "\nEnd of program.")
        sys.exit(1)

    return True