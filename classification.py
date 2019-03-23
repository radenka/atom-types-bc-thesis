from rdkit import Chem
from modules import file
import sys

# make automatic saving of molecules into separate sdf files -> DONE!
# TODO:
#   create labels for atom types based on substructures
#   create substructure text file and implement its loading into ATTYC
#
#   covert numeric representation (mol.GetSubstructureMatches() output) of atoms into atom type labels
#   and assign atom types to appropriate atoms
#


# 19. 3. consultation:
# atom types output = text file with assigned atom types to specific atoms (could be based on nested list output?)
# effective atom types assigning: detect substructures from small ones to bigger ones (e. g. C(=O) -> C(=O)[OH])
# how to do it?


class Classification:
    def __init__(self, sdf_file, classificator):
        suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
        if len(suppl) < 1:
            print("ERROR: SDF file doesn't contain any molecules. Check it and try again.\n"
                  "End of program.")
            sys.exit(1)
        self.molecules = suppl
        self.classificator = classificator
        self.sdf_name = sdf_file[:-4]
        # SDMolSupplier support indexing, class Mol *doesn't*
        # mol.GetNumAtoms(onlyExplicit=False)  # counts Hs, too

    def get_hybridization(self, atom):
        return f'{atom.GetSymbol()}#{str(atom.GetHybridization())}'

    def get_partners(self, atom):
        atom_type = atom.GetSymbol() + ":"
        neighs = []
        for neigh in atom.GetNeighbors():
            neighs.append(neigh.GetSymbol())
        # sorted(neighs) creates new list, I use list.sort() instead
        neighs.sort()
        atom_type += "".join(neighs)
        return atom_type

    def get_substructure(self, mol_idx):
        # atom_types = [None] * molecule.GetNumAtoms()
        molecule = self.molecules[mol_idx]
        make_file = False
        # print(mol_idx, Chem.MolToSmiles(molecule))

        # list where the atomic type labels will be put to
        atom_types = []

        # better transform to dictionary of SMILES: labels of atom types
        # import patterns from external txt file
        # SMARTS for carboxylic acid gives the same results as appropriate SMILES but returns the first carbon only
        # carboxylic acid '[$([$([C;$(C=[$([O;D1;$(O=C)])])]);$(C[$([O;$([H1&-0,H0&-1])])]);$(C[#6,#1])])]'
        # ether [#6]O[#6] (includes esters, too)
        # disulfide [#6]SS[#6]
        patterns = ['[#6]SS[#6]'] # C(=O)[OH]', 'C(=O)', 'C[OH]'
        all_pattern_atoms = {}
        for pattern in patterns:
            pattern_atoms = molecule.GetSubstructMatches(Chem.MolFromSmarts(pattern))
            # if mol_idx == 460:
            if len(pattern_atoms) > 0:
                print('Molecule no.', mol_idx)
                # just for my use
                all_pattern_atoms[pattern] = str(pattern_atoms)
                file.create_output_paths(self.sdf_name)
                if not make_file:
                    make_file = True

                print(pattern, pattern_atoms)
        # just for my use
        # if make_file:
        #     file.create_outputs(self.sdf_name, mol_idx, molecule, all_pattern_atoms)

        return atom_types

    # former version, before "removeHS" in supplier was set to False:
    # mol2 = Chem.AddHs(mol, addCoords=True)

    def get_atom_types(self):
        all_atom_types = []

        # utils.create_mol_from_pattern('C(=O)[OH]', 'outputs/testmol.sdf')

        for i, mol in enumerate(self.molecules):
            # print('Molecule no.', i)
            if self.classificator == "substruct":

                # Mol object doesn't have .GetName() method, necessary to keep reference through its index in SDFile
                # get_substructure() returns list of atom types for all atoms in given mol
                all_atom_types.append(self.get_substructure(i))

            else:
                atom_types = []
                for i in range(mol.GetNumAtoms()):
                    if self.classificator == "hybrid":
                        atom_types.append(self.get_hybridization(mol.GetAtomWithIdx(i)))

                    elif self.classificator == "hbo":
                        continue

                    # elif self.classificator == "substruct":
                        # atom_types.append(self.get_substructure(mol.GetAtomWithIdx(i), Chem.MolToSmiles(mol)))

                    elif self.classificator == "partners":
                        atom_types.append(self.get_partners(mol.GetAtomWithIdx(i)))
                all_atom_types.append(atom_types)

        return all_atom_types