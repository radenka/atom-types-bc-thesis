from modules.classifier import Classifier


class HBOClassifier(Classifier):

    def __init__(self, input_sdf):
        super().__init__(input_sdf)

    # TODO
    def get_hbo(self, atom):
        return f'{atom.GetSymbol()}#'

    def classify_atoms(self):
        for mol in self.supplier:
            atom_types = []
            for i in range(mol.GetNumAtoms()):
                atom_types.append(self.get_hbo(mol.GetAtomWithIdx(i)))
            self.all_atom_types.append(atom_types)