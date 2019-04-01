from rdkit import Chem
import os
import sys
import re


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


def load_atom_types():
    if not os.path.exists('SMARTS_atom_types.txt'):
        print('ERROR: SMARTS input file not found. Check if file "SMARTS_atom_types.txt" is located in the same directory'
              'as "io.py". \n'
              'End of program.')
        sys.exit(1)

    SMARTS_and_atom_types = []
    with open('SMARTS_atom_types.txt') as file:
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
                print('Check the input file carefully.')
                sys.exit(1)

        return SMARTS_and_atom_types

