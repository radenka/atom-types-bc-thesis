from .classifiers.hbo import HBOClassifier as hbo
from .classifiers.partners import PartnersClassifier as partners
from .classifiers.substruct import SubstructClassifier as substruct
from .classifiers.hybrid import HybridClassifier as hybrid
import pprint
import sys
import os


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

    if type(classifier) == str:
        classifier_name = classifier.split('.')
        if len(classifier_name) != 2 or classifier_name[0] != 'attyc':
            print(
                "ERROR: Wrong \"classifier\" argument. Use 'attyc.name_of_classifier' only').\n"
                "\nEnd of program.")
            sys.exit(1)

        try:
            classifier = getattr(sys.modules[__name__], classifier_name[-1])
        except (AttributeError, IndexError):
            print(
                "ERROR: Wrong \"classifier\" argument. Possible cause: misspelled classifier name (use "
                "'attyc.name_of_classifier only').\n"
                "Check if requested .py module located in 'classifiers' directory is imported into '__init__.py'"
                " of 'attyc' directory and that you use its alias defined in '__init__.py'."
                "\nEnd of program.")
            sys.exit(1)

    return classifier


def classify_atoms(input_sdf, classifier, file_output=False, screen_output=False):
    classifier = check_arguments(input_sdf, classifier, file_output, screen_output)
    supplier = io.load_molecules(input_sdf)

    cl = classifier()
    cl.classify_atoms(supplier)

    # outputs
    if file_output:
        print('Output will be written into a text file.')
        io.create_output_file_for_parametrization(input_sdf, cl.get_atom_types(), cl.get_name())
    if screen_output:
        print('Output will be printed on screen.')
        print_final = pprint.PrettyPrinter(indent=2)
        print('Hello!')

    return cl.get_atom_types()