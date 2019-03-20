from rdkit import Chem
from abc import abstractmethod
import sys


class Classifier:
    def __init__(self, sdf):
        suppl = Chem.SDMolSupplier(sdf, removeHs=False)
        if len(suppl) < 1:
            print("ERROR: SDF file doesn't contain any molecules. Check it and try again.\n"
                  "End of program.")
            sys.exit(1)
        self.supplier = suppl
        self.all_atom_types = []

    @abstractmethod
    def classify_atoms(self):
        pass

    def get_atom_types(self):
        return self.all_atom_types

    # DO THIS LATER
    def create_atom_types_file(self):
        # use data stored in self.atom_types, gained using classify_atoms method
        pass


