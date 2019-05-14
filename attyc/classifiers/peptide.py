from ..classifier import Classifier
from ..io import load_external_atom_types


class PeptideClassifier(Classifier):

    def __init__(self):
        super().__init__(__file__)
        self.PDB_atom_types = load_external_atom_types('PDB')
        self.standard_residues = ['GLY', 'ALA', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'TYR', 'VAL', 'CYS',
                                  'SEC', 'PRO', 'SER', 'THR', 'ASN', 'GLN', 'ARG', 'HIS', 'LYS', 'ASP',
                                  'GLU', 'PYL']

    def complete_classification(self, mol_atom_types, mol):
        """
        Classifies all atoms that have no assigned atom type yet.
        :param mol_atom_types: list of assigned atom types
        :param mol: molecule which atoms are being classified
        :return: None
        """
        for atm_idx, atom_type in enumerate(mol_atom_types):
            # some of the atom types still can be None
            # these atoms have same PDB name as standard amino acid PDB atoms, but belong to non standard residues
            if atom_type is None:
                atom_type = 'plain'

            element_symbol = mol.GetAtomWithIdx(atm_idx).GetSymbol()
            # sulphur and selen belong to one atom type, they are not distinguished by element symbol
            if atom_type != 'SSE':
                atom_type = f'{element_symbol}_{atom_type}'
            mol_atom_types[atm_idx] = atom_type

    def get_atom_types(self, mol):
        """
        Classifies atoms of molecule using 'peptide' classifier.
        :param mol: molecule (instance of class Mol)
        :return: list of assigned atom types
        """
        mol_atom_types = [None] * mol.GetNumAtoms()
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            try:
                values = self.PDB_atom_types[atom.GetPDBResidueInfo().GetName().strip()]
                for residue, assigned_atom_type in values:
                    if mol_atom_types[i] is None and (residue == atom.GetPDBResidueInfo().GetResidueName()
                            or residue == '*'):
                        mol_atom_types[i] = assigned_atom_type
            except KeyError:
                mol_atom_types[i] = 'plain'

        self.complete_classification(mol_atom_types, mol)
        return mol_atom_types

