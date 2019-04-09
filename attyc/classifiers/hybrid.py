from ..classifier import Classifier


class HybridClassifier(Classifier):
    def __init__(self):
        super().__init__(__file__)

    def get_hybridization(self, atom):
        return f'{atom.GetSymbol()}#{str(atom.GetHybridization())}'

    def classify_atoms(self, supplier):
        for mol in supplier:
            atom_types = []
            for i in range(mol.GetNumAtoms()):
                atom_types.append(self.get_hybridization(mol.GetAtomWithIdx(i)))
            self.all_atom_types.append(atom_types)
