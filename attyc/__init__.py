from .exceptions import *
from .io import create_output_file_for_parametrization, count_atom_types
import pprint
import sys
import os
import importlib
import inspect


# move to io.py?
# protein incorporation in progress
def check_arguments(input_sdf, classifier, file_output, screen_output):
    is_pdb = input_sdf.endswith('.pdb')
    # if is_pdb and classifier != 'peptide':
    #     raise InputSDFileError('input', 'Wrong usage of classifiers! PDB format supports only "peptide" classifier.')
    #
    # if not is_pdb and classifier == 'peptide':
    #     raise InputSDFileError('input', 'Wrong usage of classifiers! Do not se "peptide" classifier for SDFiles.')

    if not os.path.exists(input_sdf):
        raise InputSDFileError('input_sdf',
                               f'Check if {input_sdf} file and all path directories exist.')

    if not os.path.isfile(input_sdf):
        raise InputSDFileError('input_sdf',
                               f'Regular file {input_sdf} not found in given path.')

    if type(classifier) == str:
        current_dir = os.path.basename(os.path.dirname(__file__))

        try:
            module = importlib.import_module('.classifiers.' + classifier, current_dir)
            classifier = [m[1] for m in inspect.getmembers(module, inspect.isclass) if m[0].lower().startswith(classifier)][0]

        except ModuleNotFoundError:
            raise ClassifierNameError('classifier', f'Check if module named {classifier} is located in "classifiers" '
                                                    'directory.')
        except IndexError:
            raise ClassifierClassError("classifier", f"It seems that {classifier}.py doesn't contain subclass of class "
                                                     f"Classifier. Another possibility: Class name in {classifier}.py "
                                                     f"doesn't start with '{classifier}'.")

    else:
        raise ClassifierNameError('classifier', 'Wrong argument: Classifier name must be string!')

    if type(file_output) != bool or type(screen_output) != bool:
        raise FileScreenOutputError('output', 'Use boolean only when specifying screen/file output.')

    return classifier, is_pdb


def classify_atoms(input_sdf, classifier, file_output, screen_output):
    classifier, is_pdb = check_arguments(input_sdf, classifier, file_output, screen_output)
    supplier = io.load_molecules(input_sdf, is_pdb)
    cl = classifier()
    cl.classify_atoms(supplier, is_pdb)

    if file_output:
        print('Output will be written into a text file.')
        create_output_file_for_parametrization(input_sdf, is_pdb, cl.get_assigned_atom_types(), cl.get_name())
    if screen_output:
        print('Output (atom type statistics) will be printed on screen.')
        print_final = pprint.PrettyPrinter(indent=2)
        print_final.pprint(count_atom_types(cl.get_assigned_atom_types()))
    return cl
