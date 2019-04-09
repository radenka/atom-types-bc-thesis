from ..classifier import Classifier


class PartnersClassifier(Classifier):

    def __init__(self):
        super().__init__(__file__)

    def get_partners(self, atom):
        atom_type = atom.GetSymbol() + ":"
        neighs = []
        for neigh in atom.GetNeighbors():
            neighs.append(neigh.GetSymbol())
        # sorted(neighs) creates new list, I use list.sort() instead
        neighs.sort()
        atom_type += "".join(neighs)
        return atom_type

    def classify_atoms(self, supplier):
        for mol in supplier:
            atom_types = []
            for i in range(mol.GetNumAtoms()):
                atom_types.append(self.get_partners(mol.GetAtomWithIdx(i)))
            self.all_atom_types.append(atom_types)