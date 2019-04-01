from attyc.arguments import load_arguments, check_arguments
from attyc.classifiers.hybrid import HybridClassifier
from attyc.classifiers.hbo import HBOClassifier
from attyc.classifiers.partners import PartnersClassifier
from attyc.classifiers.substruct import SubstructClassifier
import attyc.io as IO
import pprint
import time


def perform_classification(input_sdf, classifier, file_output, screen_output):
    # TODO:
    # cutting name of sdfile from given path (/home/radenka/Documents/set01.sdf)
    # get_statistics() function - prints out detected atom types and their count

    if check_arguments(input_sdf, classifier, file_output, screen_output):
        SMARTS_and_atom_types = IO.load_atom_types()
        # very ugly! change it ASAP!
        cl = None
        if classifier == 'hybrid':   # sign '#'
            cl = HybridClassifier(input_sdf, SMARTS_and_atom_types)
        elif classifier == 'hbo':      # sign '~'
            cl = HBOClassifier(input_sdf, SMARTS_and_atom_types)
        elif classifier == 'partners':  # sign ':'
            cl = PartnersClassifier(input_sdf, SMARTS_and_atom_types)
        elif classifier == 'substruct':   # sign '*'
            cl = SubstructClassifier(input_sdf, SMARTS_and_atom_types)
        cl.classify_atoms()

        if file_output:
            print('Output will be written into a text file.')
            cl.create_atom_types_file()
        if screen_output:
            print('Output will be printed on screen.')
            print_final = pprint.PrettyPrinter(indent=2)
            print_final.pprint(cl.get_atom_types())


        # molecs = [set(mol) for mol in cl.all_atom_types]
        # ret = molecs[0]
        # for i in range(1, len(molecs)):
        #     ret = ret.union(molecs[i])
        # print(sorted(list(ret)))

        return cl.get_atom_types()


if __name__ == "__main__":
    start = time.time()
    args = load_arguments()
    if args:
        perform_classification(args.input_sdf, args.classifier, args.file_output, args.verbose)
    end = time.time()
    print(end-start)