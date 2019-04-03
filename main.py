import argparse
import attyc
import time


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
                        help="Atomic property that atom types assigned to atoms are derived from."
                             "\nSelect one of these:"
                             "\n\t'hybrid' "
                             "\n\t'hbo' "
                             "\n\t'substruct' "
                             "\n\t'partners' "
                        )
    parser.add_argument('--file_output',
                        action='store_true',
                        help='Atom types output will be written to text file.')
    parser.add_argument('--screen_output',
                        action='store_true',
                        help='Atom types output will be printed on screen.')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    start = time.time()
    args = load_arguments()
    if args:
        attyc.classify_atoms(args.input_sdf, args.classifier, args.file_output, args.screen_output)
    end = time.time()
    print('\nExecution time:', end-start)