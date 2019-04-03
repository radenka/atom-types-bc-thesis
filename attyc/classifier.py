from rdkit import Chem
from abc import abstractmethod
import sys


class Classifier:
    def __init__(self):
        self.all_atom_types = []
        self.name = None
        # used in names of output files
        # self.sdf_name = sdf[:-4]   # !!! Files that don't end with '.sdf' processed incorrectly!
        # SDMolSupplier support indexing, class Mol *doesn't*
        # mol.GetNumAtoms(onlyExplicit=False)  # counts Hs, too

    @abstractmethod
    def classify_atoms(self, supplier):
        pass

    def get_atom_types(self):
        return self.all_atom_types

    def get_name(self):
        return self.name
