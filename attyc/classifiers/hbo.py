from ..classifier import Classifier


class HBOClassifier(Classifier):

    def __init__(self):
        super().__init__(__file__)

    def create_bond_sign(self, hbo):
        bond_sign = ''
        if hbo.is_integer():
            bond_sign = str(int(hbo))
        else:
            bond_sign = 'A'
        return bond_sign

    def get_hbo(self, atom):
        hbo = 1.0
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() > hbo:
                hbo = bond.GetBondTypeAsDouble()
        atom_type = f'{atom.GetSymbol()}~{self.create_bond_sign(hbo)}'
        return atom_type

    def classify_atoms(self, supplier):
        for mol in supplier:
            atom_types = []
            for i in range(mol.GetNumAtoms()):
                atom_types.append(self.get_hbo(mol.GetAtomWithIdx(i)))
            self.all_atom_types.append(atom_types)