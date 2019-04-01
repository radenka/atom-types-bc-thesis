from rdkit import Chem
from abc import abstractmethod
import sys

class Classifier:
    def __init__(self, sdf, SMARTS_and_atom_types):
        suppl = Chem.SDMolSupplier(sdf, removeHs=False)
        if len(suppl) < 1:
            print("ERROR: SDF file doesn't contain any molecules. Check it and try again.\n"
                  "End of program.")
            sys.exit(1)
        self.supplier = suppl
        self.all_atom_types = []
        # used in names of output files
        self.sdf_name = sdf[:-4]
        # SDMolSupplier support indexing, class Mol *doesn't*
        # mol.GetNumAtoms(onlyExplicit=False)  # counts Hs, too
        self.SMARTS_and_atom_types = SMARTS_and_atom_types

    @abstractmethod
    def classify_atoms(self):
        pass

    def get_atom_types(self):
        return self.all_atom_types

    # DO THIS LATER
    def create_atom_types_file(self):
        # creates txt file of atom types assigned to atoms in given molecule set
        # use data stored in self.atom_types gained through classify_atoms() method
        pass


