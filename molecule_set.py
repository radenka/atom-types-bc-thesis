from rdkit import Chem
import sys


class MoleculeSet:
    def __init__(self, sdf_file):
        suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
        self.molecules = suppl
        # SDMolSupplier support indexing, class Mol *doesn't*
        # mol.GetNumAtoms(onlyExplicit=False)  # counts Hs, too

    def get_hybridization(self, atom):
        return f'{atom.GetSymbol()}#{str(atom.GetHybridization())}'

    def get_partners(self, atom):
        # atom.GetBonds() doesn't process bonded Hs implicitly; necessary to handle them in another way
        # if atom.GetSymbol() == 'H':
            # print((atom.GetSymbol(), atom.GetIdx()),
            #       (atom.GetNeighbors()[0].GetSymbol(), atom.GetNeighbors()[0].GetIdx()), atom.GetBonds()[0].GetBondType())
            # print((atom.GetSymbol(), atom.GetIdx()), len(atom.GetNeighbors()))

        # for bond in atom.GetBonds():
        #     print((bond.GetBeginAtom().GetSymbol(), bond.GetBeginAtom().GetIdx()),
        #           (bond.GetEndAtom().GetSymbol(), bond.GetEndAtom().GetIdx()),
        #            bond.GetBondType())

        # works fine, for each atom (Hs included) prints sth like:
        # [('N', 0), ('C', 8), ('H', 27), ('H', 28)]
        # [('C', 11), ('C', 16), ('C', 17)]
        # print((atom.GetSymbol(), atom.GetIdx()))
        # print("\t", [(neigh.GetSymbol(), neigh.GetIdx()) for neigh in atom.GetNeighbors()])
        atom_type = atom.GetSymbol() + ":"
        neighs = []
        for neigh in atom.GetNeighbors():
            neighs.append(neigh.GetSymbol())
        # sorted(neighs) creates new list, pointless, I use list.sort() instead
        neighs.sort()
        atom_type += "".join(neighs)
        return atom_type

    # former version where particular classifications weren't called in nested for cycles in "classify_atoms()" meth
    def classify_hybridization(self):
        all_atom_types = []
        for mol in self.molecules:
            # former version, before "removeHS" in supplier was set to False:
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


    def classify_group(self):
        return True

    def classify_atoms(self, classificator):
        all_atom_types = []
        for mol in self.molecules:
            # print("New molecule")
            atom_types = []
            for i in range(mol.GetNumAtoms()):
                if classificator == "hybrid":
                    atom_types.append(self.get_hybridization(mol.GetAtomWithIdx(i)))

                elif classificator == "hbo":
                    return self.classify_hbo()

                elif classificator == "group":
                    return self.classify_hbo()

                elif classificator == "partners":
                    atom_types.append(self.get_partners(mol.GetAtomWithIdx(i)))

            all_atom_types.append(atom_types)
        return all_atom_types
