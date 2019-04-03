from rdkit import Chem
import os
import sys
import re

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


def load_molecules(input_sdf):
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)
    if len(supplier) < 1:
        print("ERROR: Input SDFile doesn't contain any molecules. Check it and try again.\n"
              "End of program.")
        exit(1)
    return supplier


def load_SMARTS_and_atom_types():
    # print(os.getcwd())  # returns absolute path to atom-types-bc-thesis, why?
    # print(os.path.exists('SMARTS_atom_types.txt')) # False if the file is saved in the same directory as this script

    # print(os.path.abspath(os.path.realpath(__file__)))
    current_dir = os.path.basename(os.path.dirname(__file__))
    if not os.path.exists(os.path.join(current_dir, 'SMARTS_atom_types.txt')):
        print('ERROR: SMARTS input file not found. Check if file "SMARTS_atom_types.txt" is located in the same directory'
              'as "io.py". \n'
              'End of program.')
        sys.exit(1)

    SMARTS_and_atom_types = []
    with open(os.path.join(current_dir, 'SMARTS_atom_types.txt')) as file:
        for line in file:
            try:
                SMARTS, atom_types, functional_group = re.split(r'\s+', line, 2)
                atom_types = atom_types.split(',')
                SMARTS_and_atom_types.append((SMARTS, atom_types))

            except ValueError as err:
                print('ERROR: Something went wrong while loading SMARTS and atom types from text file.')
                if hasattr(err, 'message'):
                    print('Error caused by:')
                    print('Interpret message:', err.message)
                else:
                    print(err)
                print('Check the input file.')
                sys.exit(1)

        return SMARTS_and_atom_types


def create_output_file_for_parametrization(input_sdf, moleculeset_atom_types, classifier):
    input_filename = os.path.basename(input_sdf)
    if input_filename.endswith('.sdf'):
        input_filename = input_filename[:-4]
    output_dir = 'ATTYC_outputs'
    if not os.path.isdir(output_dir):
        print(f'Creating directory {output_dir}...')
        os.mkdir(output_dir)
    output_filename = f'{input_filename}SDF_{classifier}.txt'
    print(f'Output filename: {output_filename},\n'
          f'path: {os.path.abspath(os.path.join(output_dir, output_filename))}')

    with open(os.path.join(output_dir, output_filename), 'w') as file:
        file.writelines(','.join(mol_atom_types) + os.linesep for mol_atom_types in moleculeset_atom_types)
    print('Finished successfully.')
