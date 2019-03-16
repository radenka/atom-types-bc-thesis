from rdkit import Chem
import sys


class MoleculeSet:
    def __init__(self, sdf_file):
        suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
        self.molecules = suppl
        # SDMolSupplier support indexing, class Mol *doesn't*
        # mol.GetNumAtoms(onlyExplicit=False)  # counts Hs, too

    def classify_hybridization(self):
        all_atom_types = []
        for mol in self.molecules:
            # later version, before "removeHS" in supplier was set to False:
            # molecules in SDMolSupplier don't contain Hs implicitly; it is necesarry to add them manually
            # mol2 = Chem.AddHs(mol, addCoords=True)
            atom_types = []
            for i in range(mol.GetNumAtoms()):
                atom_types.append(f'{mol.GetAtomWithIdx(i).GetSymbol()}#{str(mol.GetAtomWithIdx(i).GetHybridization())}')
            all_atom_types.append(atom_types)
        return all_atom_types
        # less comprehensible version:
        # return  [[f'{mol.GetAtomWithIdx(i).GetSymbol()}#{str(mol.GetAtomWithIdx(i).GetHybridization())}'
        #      for i in range(mol.GetNumAtoms())] for mol in self.molecules]

    def classify_hbo(self):
        return True

    def classify_partner(self):
        return True


    def classify_group(self):
        return True

    def classify_atoms(self, classificator):
        if classificator == "hybrid":
            return self.classify_hybridization()

        elif classificator == "hbo":
            return self.classify_hbo()

        elif classificator == "group":
            return self.classify_hbo()

        elif classificator == "partner":
            return self.classify_hbo()

