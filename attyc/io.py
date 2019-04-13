from rdkit import Chem
import os
import re
from pathlib import Path
from collections import Counter
from .exceptions import InputSMARTSError

#### just for my use ###########################

def create_output_paths(sdf_name):
    if not os.path.isdir('outputs'):
        os.mkdir('outputs')
    if not os.path.isfile(f'outputs/{sdf_name}_patterns.txt'):
        os.mknod(f'outputs/{sdf_name}_patterns.txt')


def create_substruct_outputs(sdf_name, mol_idx, molecule):
    # with open(f'outputs/{sdf_name}_patterns.txt', 'a') as f:
    #     f.write(f'Molecule n. {mol_idx}\n')
    #     for pattern, atoms in all_pattern_atoms.items():
    #         f.write(f'{pattern}: {atoms}\n')
    #     f.write('\n')
    make_sdf(mol_idx, molecule)


def make_sdf(mol_idx, molecule):
    writer = Chem.SDWriter(f'outputs/mol{mol_idx}.sdf')
    writer.write(molecule)

################################################


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


def load_molecules(input_sdf):
    """
    :param input_sdf: sdfile containing molecules to classify
    :return: SDMolSupplier containing Mol instances of molecules in given sdfile
    """
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    if len(supplier) < 1:
        print("Input SDFile doesn't contain any molecules. Check it and try again.\n"
              "End of program.")
        exit(1)
    return supplier


def load_SMARTS_and_atom_types():
    current_dir = os.path.dirname(__file__)
    # print('Current directory of io.py:', current_dir) -> attyc
    if not os.path.exists(os.path.join(current_dir, 'SMARTS_atom_types.txt')):
        raise InputSMARTSError('SMARTS input file not found. Check if file "SMARTS_atom_types.txt" is located '
                               'in the same directory as "io.py".')

    SMARTS_and_atom_types = []
    with open(os.path.join(current_dir, 'SMARTS_atom_types.txt')) as file:
        for line in file:
            try:
                SMARTS, atom_types, _ = re.split(r'\s+', line, 2)
                atom_types = atom_types.split(',')
                SMARTS_and_atom_types.append((SMARTS, atom_types))

            except ValueError:
                raise InputSMARTSError('Something went wrong while loading SMARTS and atom types from text file.'
                                       ' Check the input file (each line should contain at least 2 columns). ')
        return SMARTS_and_atom_types


def create_output_file_for_parametrization(input_sdf, set_atom_types, classifier_name):
    input_filename = os.path.basename(input_sdf)
    if input_filename.endswith('.sdf'):
        input_filename = input_filename[:-4]
    # returns parent directory of directory where io.py is saved
    parent_dir = Path(__file__).resolve().parents[1]
    output_dirname = 'ATTYC_outputs'
    if not os.path.isdir(os.path.join(parent_dir, output_dirname)):
        print(f'Creating directory {output_dirname}...')
        os.mkdir(os.path.join(parent_dir, output_dirname))
    output_filename = f'{input_filename}SDF_{classifier_name}.txt'
    print(f'Output filename: {output_filename},\n'
          f'path: {os.path.join(parent_dir, output_dirname, output_filename)}')
    with open(os.path.join(parent_dir, output_dirname, output_filename), 'w') as file:
        file.writelines(','.join(mol_atom_types) + os.linesep for mol_atom_types in set_atom_types)
    print('Finished successfully.')


def count_atom_types(set_atom_types):
    counter = Counter()
    for mol_atom_types in set_atom_types:
        for atom_type in mol_atom_types:
            counter[atom_type] += 1
    return counter
