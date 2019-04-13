from ..classifier import Classifier


class HybridClassifier(Classifier):
    def __init__(self):
        super().__init__(__file__)

    def get_hybridization(self, atom):
        return f'{atom.GetSymbol()}#{str(atom.GetHybridization())}'

    def get_atom_types(self, mol):
        mol_atom_types = []
        for i in range(mol.GetNumAtoms()):
            mol_atom_types.append(self.get_hybridization(mol.GetAtomWithIdx(i)))
        return mol_atom_types

