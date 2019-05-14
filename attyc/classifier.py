import os
from abc import abstractmethod


class Classifier:
    def __init__(self, classifier_path):
        self.assigned_atom_types = []
        self.name = os.path.basename(classifier_path)[:-3]

    @abstractmethod
    def get_atom_types(self, mol):
        """
        Classifies atoms of molecule using specific classifier.
        """
        pass

    def classify_atoms(self, supplier, is_pdb):
        """
        Runs atom types classification.
        :param supplier: contains molecule(s) to classify
        :param is_pdb: boolean to distinguish SDF and PDB file
        :return: None
        """
        if is_pdb:
            # supplier created from PDB file is instance of class Mol
            self.assigned_atom_types.append(
                self.get_atom_types(supplier)
            )
        else:
            for i, mol in enumerate(supplier):
                self.assigned_atom_types.append(
                    self.get_atom_types(mol)
                )

    def get_assigned_atom_types(self):
        return self.assigned_atom_types

    def get_name(self):
        return self.name

