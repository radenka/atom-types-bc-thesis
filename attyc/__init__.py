from .classifiers.hbo import HBOClassifier as hbo
from .classifiers.partners import PartnersClassifier as partners
from .classifiers.substruct import SubstructClassifier as substruct
from .classifiers.hybrid import HybridClassifier as hybrid
from .exceptions import *

import pprint
import sys
import os


# move to io.py?
def check_arguments(input_sdf, classifier, file_output, screen_output):
    if not os.path.exists(input_sdf):
        raise InputSDFileError('input_sdf',
                               f'Check if {input_sdf} file and all path directories exist.')

    if not os.path.isfile(input_sdf):
        raise InputSDFileError('input_sdf',
                               f'Regular file {input_sdf} not found in given path.')

    if type(classifier) == str:
        classifier_name = classifier.split('.')
        if len(classifier_name) != 2 or classifier_name[0] != 'attyc':
            raise ClassifierNameError('classifier', 'Use "attyc.name_of_classifier" only!')

        try:
            # print(sys.modules[__name__]), possible to call sys.modules.keys()
            # <module 'attyc' from '/home/radenka/Dokumenty/SKOLA/6_semestr_Bc/BcThesis/atom-types-bc-thesis/attyc/__init__.py'>
            classifier = getattr(sys.modules[__name__], classifier_name[-1])
        except (AttributeError, IndexError):
            # dotaz ####
            raise ClassifierNameError('classifier', 'Possible cause: misspelled classifier name (use '
                                      '"attyc.name_of_classifier" only). Check if requested .py module located in '
                                      '"classifiers" directory is imported into "__init__.py" in "attyc" directory '
                                      'and that you use its alias defined in "__init__.py".')
    else:
        raise AttributeError('Wrong argument: Classifier name must be string!')

    if type(file_output) != bool or type(screen_output) != bool:
        raise FileScreenOutputError('output', 'Use boolean only.')

    return classifier


def classify_atoms(input_sdf, classifier, file_output=False, screen_output=False):
    classifier = check_arguments(input_sdf, classifier, file_output, screen_output)
    supplier = io.load_molecules(input_sdf)

    cl = classifier()
    cl.classify_atoms(supplier)

    if file_output:
        print('Output will be written into a text file.')
        io.create_output_file_for_parametrization(input_sdf, cl.get_assigned_atom_types(), cl.get_name())
    if screen_output:
        print('Output will be printed on screen.')
        print_final = pprint.PrettyPrinter(indent=2)
        print_final.pprint(io.count_atom_types(cl.get_assigned_atom_types()))

    return cl.get_assigned_atom_types()