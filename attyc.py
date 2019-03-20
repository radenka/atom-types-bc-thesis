from modules.arguments import load_arguments, check_arguments
from modules.classifiers.hybrid import HybridClassifier
from modules.classifiers.hbo import HBOClassifier
from modules.classifiers.partners import PartnersClassifier
import pprint


def perform_classification(input_sdf, classifier, output_file):
    # TODO:
    # add file_output check (if ATTYC is being called from another script and user misspells "True" or "False")
    if check_arguments(input_sdf, classifier):
        cl = None
        if classifier == 'hybrid':
            cl = HybridClassifier(input_sdf)
        if classifier == 'hbo':
            cl = HBOClassifier(input_sdf)
        if classifier == 'partners':
            cl = PartnersClassifier(input_sdf)
        # if classifier == 'substruct':
        #     cl = SubstructClassifier(input_sdf)
        # if classifier == 'global':
        #     cl = GlobalClassifier(input_sdf)
        cl.classify_atoms()

        if output_file:
            cl.create_atom_types_file()
        else:
            print_final = pprint.PrettyPrinter(indent=2)
            print_final.pprint(cl.get_atom_types())


if __name__ == "__main__":
    args = load_arguments()
    if args:
        perform_classification(args.input_sdf, args.classifier, args.output_file)