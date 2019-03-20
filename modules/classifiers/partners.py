from modules.classifier import Classifier


class PartnersClassifier(Classifier):

    def __init__(self, input_sdf):
        super().__init__(input_sdf)

    def get_partners(self, atom):
        atom_type = atom.GetSymbol() + ":"
        neighs = []
        for neigh in atom.GetNeighbors():
            neighs.append(neigh.GetSymbol())
        # sorted(neighs) creates new list, I use list.sort() instead
        neighs.sort()
        atom_type += "".join(neighs)
        return atom_type

    def classify_atoms(self):
        for mol in self.supplier:
            atom_types = []
            for i in range(mol.GetNumAtoms()):
                atom_types.append(self.get_partners(mol.GetAtomWithIdx(i)))
            self.all_atom_types.append(atom_types)