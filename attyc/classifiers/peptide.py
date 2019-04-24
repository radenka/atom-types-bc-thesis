from ..classifier import Classifier
from ..io import load_external_atom_types
from collections import Counter


class PeptideClassifier(Classifier):

    def __init__(self):
        super().__init__(__file__)
        self.PDB_atom_types = load_external_atom_types('PDB')
        self.standard_residues = ['GLY', 'ALA', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR', 'VAL', 'CYS',
                                  'SEC', 'PRO', 'SER', 'THR', 'ASN', 'GLN', 'ARG', 'HIS', 'LYS', 'ASP',
                                  'GLU', 'PYL']

    def complete_classification(self, mol_atom_types, mol):
        # delete later
        counter = Counter()
        for atm_idx, atom_type in enumerate(mol_atom_types):
            # some of the atom types still can be None
            # these atoms have same PDB name as standard amino acid PDB atoms, but belong to non standard residues
            if atom_type is None:
                atom_type = 'plain'

            element_symbol = mol.GetAtomWithIdx(atm_idx).GetSymbol()
            # sulphur and selen belong to one atom type, they are not distinguished by element symbol
            if element_symbol != 'S' and element_symbol != 'SE':
                atom_type = f'{element_symbol}_{atom_type}'

            mol_atom_types[atm_idx] = atom_type
            # print(atom_type)
            counter[atom_type] += 1
        # print(counter)

    def get_atom_types(self, mol):
        mol_atom_types = [None] * mol.GetNumAtoms()  #
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            try:
                # print(atom.GetPDBResidueInfo().GetName().strip())
                values = self.PDB_atom_types[atom.GetPDBResidueInfo().GetName().strip()]
                if atom.GetPDBResidueInfo().GetName().strip() == 'SE':
                    print(atom.GetSymbol())
                # print(values)
                for residue, assigned_atom_type in values:
                    if mol_atom_types[i] is None and (residue == atom.GetPDBResidueInfo().GetResidueName()
                            or residue == '*'):
                        mol_atom_types[i] = assigned_atom_type
            except KeyError:
                mol_atom_types[i] = 'plain'

        self.complete_classification(mol_atom_types, mol)
        # print(mol_atom_types)
        return mol_atom_types

