from modules.arguments import load_arguments, check_arguments
from modules.classifiers.hybrid import HybridClassifier
from modules.classifiers.hbo import HBOClassifier
from modules.classifiers.partners import PartnersClassifier
from modules.classifiers.substruct import SubstructClassifier
import pprint


def perform_classification(input_sdf, classifier, output_file, screen_output):
    # TODO:
    # add file_output check (if ATTYC is being called from another script and user misspells "True" or "False")
    # cutting name of sdfile from given path (/home/radenka/Documents/set01.sdf)
    # get_statistics() function - prints out detected atom types and their count

    if check_arguments(input_sdf, classifier):
        cl = None
        if classifier == 'hybrid':
            cl = HybridClassifier(input_sdf)
        if classifier == 'hbo':
            cl = HBOClassifier(input_sdf)
        if classifier == 'partners':
            cl = PartnersClassifier(input_sdf)
        if classifier == 'substruct':
            cl = SubstructClassifier(input_sdf)
        cl.classify_atoms()

        if output_file:
            print('Output will be written into a text file.')
            cl.create_atom_types_file()
        if screen_output:
            print('Output will be printed on screen.')
            print_final = pprint.PrettyPrinter(indent=2)
            print_final.pprint(cl.get_atom_types())


        molecs = [set(mol) for mol in cl.all_atom_types]
        ret = molecs[0]
        for i in range(1, len(molecs)):
            ret = ret.union(molecs[i])
        print(sorted(list(ret)))

        return cl.get_atom_types()


if __name__ == "__main__":
    args = load_arguments()
    if args:
        perform_classification(args.input_sdf, args.classifier, args.output_file, args.verbose)