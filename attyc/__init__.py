import pprint
from .io import check_arguments, create_output_file_for_parametrization, count_atom_types


def classify_atoms(input_file, classifier, file_output, screen_output):
    """
    Runs atom classification. Function is called using attyc.classify_atoms(...) command.
    :param input_file: molecule file with atoms to classify
    :param classifier: chosen classifier
    :param file_output, screen_output: boolean to specify type of classification output
    :return: None
    """
    classifier, is_pdb = check_arguments(input_file, classifier, file_output, screen_output)
    supplier = io.load_molecules(input_file, is_pdb)
    cl = classifier()
    cl.classify_atoms(supplier, is_pdb)

    if file_output:
        print('Output will be written into a text file.')
        create_output_file_for_parametrization(input_file, is_pdb, cl.get_assigned_atom_types(), cl.get_name())
    if screen_output:
        print('Output (atom type statistics) will be printed on screen.')
        print_final = pprint.PrettyPrinter(indent=2)
        print_final.pprint(count_atom_types(cl.get_assigned_atom_types()))

