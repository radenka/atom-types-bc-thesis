from rdkit import Chem
from ..classifier import Classifier
from ..io import load_external_atom_types


class SubstructClassifier(Classifier):

    def __init__(self):
        super().__init__(__file__)
        self.SMARTS_and_atom_types = load_external_atom_types('SMARTS')

    def classifify_remaining_H(self, hydrogen):
        """
        Assigns atom type to not classified hydrogen atom.
        :param hydrogen: hydrogen atom
        :return: assigned atom type
        """
        return '-' + str(hydrogen.GetNeighbors()[0].GetSymbol())

        # UNCOMMENT following lines and COMMENT "initial classification" above to use simplified classification
        # neigh = str(hydrogen.GetNeighbors()[0].GetSymbol())
        # if neigh == 'C' or neigh == 'N':
        #     return neigh
        # return '1bond'

    def complete_classification(self, atom_types, molecule):
        for atm_idx, atom_type in enumerate(atom_types):
            if atom_type is None:
                if molecule.GetAtomWithIdx(atm_idx).GetAtomicNum() == 1:
                    atom_type = self.classifify_remaining_H(molecule.GetAtomWithIdx(atm_idx))
                else:
                    atom_type = 'plain'
            atom_type = molecule.GetAtomWithIdx(atm_idx).GetSymbol() + '*' + atom_type
            atom_types[atm_idx] = atom_type

    def analyze_aromatic_rings(self, molecule, atom_types):
        """
        Assigns corresponding atom types to aromatic atoms.
        :param molecule: molecule (instance of class Mol)
        :param atom_types: list to store assigned atom types of molecule
        :return: None
        """
        rings = molecule.GetRingInfo().AtomRings()
        for ring in rings:
            for atm_idx in ring:
                atom = molecule.GetAtomWithIdx(atm_idx)
                if atom.GetIsAromatic() and atom_types[atm_idx] is None:
                    atom_types[atm_idx] = 'A'
                    for neigh in atom.GetNeighbors():
                        # aromatic Hs detection
                        if neigh.GetAtomicNum() == 1:
                            # initial classification
                            atom_types[neigh.GetIdx()] = f'-{atom.GetSymbol()}A'

                            # UNCOMMENT following line command and COMMENT "initial classification" above
                            # to use simplified classification
                            # atom_types[neigh.GetIdx()] = f'{atom.GetSymbol()}'

    def get_atom_types(self, mol):
        """
        Classifies atoms of molecule using 'substruct' classifier.
        :param mol: molecule (instance of class Mol)
        :return: list of assigned atom types
        """
        mol_atom_types = [None] * mol.GetNumAtoms()
        # aromatic atoms are detected first
        self.analyze_aromatic_rings(mol, mol_atom_types)
        for pattern, loaded_atom_types in self.SMARTS_and_atom_types:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                pattern_atoms = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
                for atom_tuple in pattern_atoms:
                    for atm_idx, atom_type in zip(atom_tuple, loaded_atom_types):
                        if mol_atom_types[atm_idx] is None:
                            mol_atom_types[atm_idx] = atom_type

        self.complete_classification(mol_atom_types, mol)
        return mol_atom_types
