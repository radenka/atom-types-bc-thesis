from abc import abstractmethod
import os


class Classifier:
    def __init__(self, classifier_path):
        self.assigned_atom_types = []
        self.name = os.path.basename(classifier_path)[:-3]

    @abstractmethod
    def get_atom_types(self, mol):
        pass

    def classify_atoms(self, supplier, is_pdb):
        if is_pdb:
            # because supplier created from PDB file is one instance of class Mol
            self.assigned_atom_types.append(
                self.get_atom_types(supplier)
            )
        else:
            for i, mol in enumerate(supplier):
                # print('Molecule no.', i)
                self.assigned_atom_types.append(
                    self.get_atom_types(mol)
                )

    def get_assigned_atom_types(self):
        return self.assigned_atom_types

    def get_name(self):
        return self.name

