from ..classifier import Classifier


class HybridClassifier(Classifier):
    def __init__(self):
        super().__init__(__file__)

    def get_hybridization(self, atom):
        """
        Returns assigned hybridization atom type of atom.
        """
        return f'{atom.GetSymbol()}#{str(atom.GetHybridization())}'

    def get_atom_types(self, mol):
        """
        Classifies atoms of molecule using 'hbo' classifier.
        :param mol: molecule (instance of class Mol)
        :return: list of assigned atom types
        """
        mol_atom_types = []
        for i in range(mol.GetNumAtoms()):
            mol_atom_types.append(self.get_hybridization(mol.GetAtomWithIdx(i)))
        return mol_atom_types

