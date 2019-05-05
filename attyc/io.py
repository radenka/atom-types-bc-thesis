from rdkit import Chem
import os, sys
import re
from pathlib import Path
from collections import Counter
from .exceptions import InputSMARTSError, InputSDFileError

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


def load_molecules(input_sdf, is_pdb):
    """
    :param input_sdf: sdfile contains molecules of atoms to classify
    :return: SDMolSupplier contains Mol instances of molecules of given sdfile
    """
    supplier = None
    if is_pdb:
        supplier = Chem.MolFromPDBFile(input_sdf, removeHs=False)
    else:
        supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
        if len(supplier) < 1:
            raise InputSDFileError("input_file", f"Input SDFile {input_sdf} doesn't contain any molecules. "
                                                 f"Check it and try again.")
    return supplier


def load_SMARTS(file):
    SMARTS_and_atom_types = []
    for line in file:
        try:
            SMARTS, atom_types, _ = re.split(r'\s+', line, 2)
            atom_types = atom_types.split(',')
            SMARTS_and_atom_types.append((SMARTS, atom_types))

        except ValueError:
            raise InputSMARTSError('smarts', 'Something went wrong while loading SMARTS and atom types from text file.'
                                   ' Check the input file (each line should contain at least 2 columns). ')
    return SMARTS_and_atom_types


def load_PDB(file):
    PDB_atom_types_supplier = {}
    for line in file:
        PDB_atom_type, residues_with_assigned_atom_types = re.split(r'\s+', line, 1)
        residues_with_assigned_atom_types = residues_with_assigned_atom_types.strip()
        values = []
        for residue_and_atom_type in residues_with_assigned_atom_types.split(','):
            residue, defined_atom_type = residue_and_atom_type.split(':')
            values.append((residue, defined_atom_type))
            # print(PDB_atom_type, values)
        PDB_atom_types_supplier[PDB_atom_type] = values
    return PDB_atom_types_supplier


def load_external_atom_types(notation):
    current_dir = os.path.dirname(__file__)
    if not os.path.exists(os.path.join(current_dir, notation + '_atom_types.txt')):
        # rewrite to ExternalTypesInputFileError
        raise InputSMARTSError('input', f'{notation} input file not found. Check if file "{notation}_atom_types.txt" '
                                        f'is located in the same directory as "io.py".')

    with open(os.path.join(current_dir, notation + '_atom_types.txt')) as file:
        if notation == 'SMARTS':
            return load_SMARTS(file)
        return load_PDB(file)


def load_SMARTS_and_atom_types():
    current_dir = os.path.dirname(__file__)
    # print('Current directory of io.py:', current_dir) -> attyc
    if not os.path.exists(os.path.join(current_dir, 'SMARTS_atom_types.txt')):
        raise InputSMARTSError('smarts','SMARTS input file not found. Check if file "SMARTS_atom_types.txt" is located '
                               'in the same directory as "io.py".')

    SMARTS_and_atom_types = []
    with open(os.path.join(current_dir, 'SMARTS_atom_types.txt')) as file:
        for line in file:
            try:
                SMARTS, atom_types, _ = re.split(r'\s+', line, 2)
                atom_types = atom_types.split(',')
                SMARTS_and_atom_types.append((SMARTS, atom_types))

            except ValueError:
                raise InputSMARTSError('smarts', 'Something went wrong while loading SMARTS and atom types from text file.'
                                       ' Check the input file (each line should contain at least 2 columns). ')
        return SMARTS_and_atom_types


def create_output_file_for_parametrization(input_sdf, is_pdb, set_atom_types, classifier_name):
    input_filename = os.path.basename(input_sdf)
    if input_filename.endswith('.sdf') or input_filename.endswith('.pdb'):
        input_filename = input_filename[:-4]
    # returns parent directory of directory where io.py is saved
    parent_dir = Path(__file__).resolve().parents[1]
    output_dirname = 'ATTYC_outputs'
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
        file.writelines(','.join(mol_atom_types) + os.linesep for mol_atom_types in set_atom_types)
    print('Finished successfully.')


def count_atom_types(set_atom_types):
    counter = Counter()
    for mol_atom_types in set_atom_types:
        for atom_type in mol_atom_types:
            counter[atom_type] += 1
    print(len(counter.items()))
    return counter
