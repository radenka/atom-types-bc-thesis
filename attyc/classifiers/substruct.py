from attyc.classifier import Classifier
from rdkit import Chem
from ..io import create_output_paths, create_substruct_outputs
from ..utils import create_mol_from_pattern


class SubstructClassifier(Classifier):
    def __init__(self, input_sdf):
        super().__init__(input_sdf)

    def analyze_aromatic_rings(self, molecule, atom_types):
        rings = molecule.GetRingInfo().AtomRings()
        for ring in rings:
            # tuples = []
            for atm_idx in ring:
                atom = molecule.GetAtomWithIdx(atm_idx)
                # print(atom.GetSymbol(), atm_idx)
                if atom.GetIsAromatic():
                    atom_types[atm_idx] = f'{atom.GetSymbol()}$A'
                    # implicit and explicit Hs not detected, necessary to use brute force
                    # methods Atom.GetNumExplicitHs(), .GetNumImplicitHs(), atom.GetTotalNumHs()) return 0 for each atom
                    for neigh in atom.GetNeighbors():
                        # print('\t' + neigh.GetSymbol() + ' ' + str(neigh.GetIdx()))
                        # aromatic Hs detection
                        if neigh.GetAtomicNum() == 1:
                            atom_types[neigh.GetIdx()] = 'H$A'
                # tuples.append((atom.GetSymbol(), atm_idx, atom.GetIsAromatic()))
            # print(tuples)

    def get_substructures(self, mol_idx, molecule):
        # list where the atom type labels will be put to
        atom_types = [None] * molecule.GetNumAtoms()
        make_file = False

        # import patterns from external txt file
        # OrderedDict?
        # IDEA: adding second atom type to already assigned one (e.g. C$A_connectedX) # useful?
        # self.analyze_aromatic_rings(molecule, atom_types)
        patterns = ['N=S']    # C(=O)[OH]', 'C(=O)', 'C[OH]', 'C(=O)[H]', '[SX2]' '[O]'
        all_pattern_atoms = {}
        for pattern in patterns:
            if molecule.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                print('Molecule no.', mol_idx)
                pattern_atoms = molecule.GetSubstructMatches(Chem.MolFromSmarts(pattern))
                # just for my use
                all_pattern_atoms[pattern] = str(pattern_atoms)
                # create_output_paths(self.sdf_name)
                if not make_file:
                    make_file = True
                # connect = molecule.GetAtomWithIdx(pattern_atoms[0][0])
                # print(connect.GetSymbol())
                print(pattern, pattern_atoms)

            # just for my own use
            # if make_file:
            #     create_substruct_outputs(self.sdf_name, mol_idx, molecule, all_pattern_atoms)

        # !!! get molecule name: mol.GetProp('_Name') -> returns NSC_1000089, NSC_100992 etc.
        return atom_types

    # former version, before "removeHS" in supplier was set to False:
    # mol2 = Chem.AddHs(mol, addCoords=True)

    def classify_atoms(self):
        # create_mol_from_pattern('C[CX3](=[OX1])[OX2]C', 'outputs/testmol_modified_ester.sdf')

        for i, mol in enumerate(self.supplier):
            # print('Molecule no.', i)
                self.all_atom_types.append(
                    self.get_substructures(i, mol)
                )