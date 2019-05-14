import importlib
import inspect
import os
import re
from pathlib import Path
from .exceptions import *
from rdkit import Chem
from collections import Counter
from .exceptions import ExternalTypesInputFileError, InputFileError


def check_arguments(input_file, classifier, file_output, screen_output):
    """
    Validates input arguments.
    :param input_file: path to input file
    :param classifier: name of required classifier to use
    :param file_output, screen_output: boolean to specify type of classification output
    :return: tuple containing appropriate Classifier subclass and boolean value is_pdb
    """
    if not os.path.exists(input_file):
        raise InputFileError('input_file',
                             f'Check if {input_file} file and all path directories exist.')

    if not os.path.isfile(input_file):
        raise InputFileError('input_file',
                             f'Regular file {input_file} not found in given path.')

    if not input_file.endswith('.sdf') and not input_file.endswith('.pdb'):
        raise InputFileError('input_file', f'Name of file "{input_file}" doesn\'t end with ".sdf" nor ".pdb". Rename '
                                           f'file appropriately to run classification.')
    is_pdb = input_file.endswith('.pdb')

    if type(classifier) == str:
        current_dir = os.path.basename(os.path.dirname(__file__))

        if not is_pdb and classifier == 'peptide':
            raise ClassifierNameError('classifier', '"Peptide" classifier can\'t be used to classify SDF files.')

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


def get_available_classifiers():
    """
    :return: string containing names of .py files located in 'classifiers' directory
    """
    classifiers = ""
    for file in os.listdir(os.path.join(os.path.dirname(__file__), 'classifiers')):
        # ignores automatically generated python file __init__ and __pycache__
        if not file.startswith('__') and file.endswith('.py'):
            classifiers += f'\'{file[:-3]}\' '
    return classifiers


def load_molecules(input_file, is_pdb):
    """
    :param input_file: file contains molecules to classify
    :param is_pdb: boolean to distinguish SDF and PDB file
    :return: instance of class Mol for PDB files, SDMolSupplier for SDF files.
    """
    supplier = None
    if is_pdb:
        supplier = Chem.MolFromPDBFile(input_file, removeHs=False)
    else:
        supplier = Chem.SDMolSupplier(input_file, removeHs=False)
        if len(supplier) < 1:
            raise InputFileError("input_file", f"Input SDFile {input_file} doesn't contain any molecules. "
                                               f"Check it and try again.")
    return supplier


def load_SMARTS(file):
    """
    Loades SMARTS and corresponding atom types from external text file
    :param file: file to load SMARTS and atom types from
    :return: list of tuples that contain SMARTS and corresponding defined atom types
    """
    SMARTS_and_atom_types = []
    for line in file:
        try:
            SMARTS, atom_types, _ = re.split(r'\s+', line, 2)
            atom_types = atom_types.split(',')
            SMARTS_and_atom_types.append((SMARTS, atom_types))

        except ValueError:
            raise ExternalTypesInputFileError('smarts', 'Something went wrong while loading SMARTS and atom types from text file.'
                                   ' Check the input file (each line should contain at least 2 columns). ')
    return SMARTS_and_atom_types


def load_PDB(file):
    """
    Loads IUPAC amino acid atom types and corresponding designed atom types according to structural motif membership
    :param file: file to extract data from
    :return: dictionary
    """
    PDB_atom_types = {}
    for line in file:
        PDB_atom_type, residues_with_assigned_atom_types = re.split(r'\s+', line, 1)
        residues_with_assigned_atom_types = residues_with_assigned_atom_types.strip()
        values = []
        for residue_and_atom_type in residues_with_assigned_atom_types.split(','):
            residue, defined_atom_type = residue_and_atom_type.split(':')
            values.append((residue, defined_atom_type))
        PDB_atom_types[PDB_atom_type] = values
    return PDB_atom_types


def load_external_atom_types(notation):
    """
    Loads data from external files to provide "substruct" and "peptide" classification
    :param notation: string that specifies which external file extract data from
    :return: list (SMARTS) or dictionary (PDB) containing data necessary for "substruct" or "peptide" classification
    """
    current_dir = os.path.dirname(__file__)
    if not os.path.exists(os.path.join(current_dir, notation + '_atom_types.txt')):
        raise ExternalTypesInputFileError('external_file',
                               f'{notation} input file not found. Check if file "{notation}_atom_types.txt" '
                               f'is located in the same directory as "io.py".')

    with open(os.path.join(current_dir, notation + '_atom_types.txt')) as file:
        if notation == 'SMARTS':
            return load_SMARTS(file)
        return load_PDB(file)


def create_output_file_for_parametrization(input_file, is_pdb, molset_atom_types, classifier_name):
    """
    Creates file with assigned atom types with respect to classifier used. Writes one molecule per line.
    :param input_file: path of source molfile which atoms are to classify
    :param is_pdb: boolean to distinguish SDF and PDB file
    :param molset_atom_types: list of all assigned atom types in molecule set
    :param classifier_name: classifier used
    :return: None
    """
    input_filename = os.path.basename(input_file)
    if input_filename.endswith('.sdf') or input_filename.endswith('.pdb'):
        input_filename = input_filename[:-4]
    # returns parent directory of directory where io.py is saved
    parent_dir = Path(__file__).resolve().parents[1]
    output_dirname = 'attyc_outputs'
    if not os.path.isdir(os.path.join(parent_dir, output_dirname)):
        print(f'Creating directory {output_dirname}...')
        os.mkdir(os.path.join(parent_dir, output_dirname))

    if is_pdb:
        file_extension = 'PDB'
    else:
        file_extension = 'SDF'
    output_filename = f'{input_filename}{file_extension}_{classifier_name}.txt'
    print(f'Output filename: {output_filename},\n'
          f'path: {os.path.join(parent_dir, output_dirname, output_filename)}')
    with open(os.path.join(parent_dir, output_dirname, output_filename), 'w') as file:
        file.writelines(','.join(mol_atom_types) + os.linesep for mol_atom_types in molset_atom_types)
    print('Finished successfully.')


def count_atom_types(molset_atom_types):
    """
    Counts assigned atom types of given molecule set.
    :param molset_atom_types: list of lists that contain assigned atom types of molecules in moleculeset
    :return: counter containing statistics of atom types of moleculeset
    """
    counter = Counter()
    for mol_atom_types in molset_atom_types:
        for atom_type in mol_atom_types:
            counter[atom_type] += 1
    return counter
