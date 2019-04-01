from attyc.classifier import Classifier


class HBOClassifier(Classifier):

    def __init__(self, input_sdf, SMARTS_and_atom_types):
        super().__init__(input_sdf, SMARTS_and_atom_types)

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

    def classify_atoms(self):
        for mol in self.supplier:
            atom_types = []
            for i in range(mol.GetNumAtoms()):
                atom_types.append(self.get_hbo(mol.GetAtomWithIdx(i)))
            self.all_atom_types.append(atom_types)