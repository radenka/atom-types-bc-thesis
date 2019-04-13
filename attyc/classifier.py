from abc import abstractmethod
import os


class Classifier:
    def __init__(self, classifier_path):
        self.assigned_atom_types = []
        self.name = os.path.basename(classifier_path)[:-3]
        # SDMolSupplier support indexing, class Mol *doesn't*
        # mol.GetNumAtoms(onlyExplicit=False)  # counts Hs, too

    @abstractmethod
    def get_atom_types(self, mol):
        pass

    def classify_atoms(self, supplier):
        for mol in supplier:
            # print('Molecule no.', i)
            self.assigned_atom_types.append(
                self.get_atom_types(mol)
            )

    def get_assigned_atom_types(self):
        return self.assigned_atom_types

    def get_name(self):
        return self.name

