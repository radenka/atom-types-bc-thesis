from ..classifier import Classifier


class HBOClassifier(Classifier):

    def __init__(self):
        super().__init__(__file__)

    def create_bond_sign(self, hbo):
        """
        :param hbo: character that represents highest bond order of arbitrary atom
        :return: bond sign later used for atom sign creation
        """
        bond_sign = ''
        if hbo.is_integer():
            bond_sign = str(int(hbo))
        else:
            bond_sign = 'A'
        return bond_sign

    def get_hbo(self, atom):
        """
        Creates atom type sign.
        :param atom: classified atom
        :return: atom type sign
        """
        hbo = 1.0
        for bond in atom.GetBonds():    
            if bond.GetBondTypeAsDouble() > hbo:
                hbo = bond.GetBondTypeAsDouble()
        atom_type = f'{atom.GetSymbol()}~{self.create_bond_sign(hbo)}'
        return atom_type

    def get_atom_types(self, mol):
        """
        Classifies atoms of molecule using 'hbo' classifier.
        :param mol: molecule (instance of class Mol)
        :return: list of assigned atom types
        """
        mol_atom_types = []
        for i in range(mol.GetNumAtoms()):
            mol_atom_types.append(self.get_hbo(mol.GetAtomWithIdx(i)))
        return mol_atom_types
