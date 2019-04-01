from attyc.classifier import Classifier
from rdkit import Chem
from collections import Counter
from ..io import create_output_paths, create_substruct_outputs
from ..utils import create_mol_from_pattern


class SubstructClassifier(Classifier):
    def __init__(self, input_sdf, SMARTS_and_atom_types):
        super().__init__(input_sdf, SMARTS_and_atom_types)

    def classifify_remaining_H(self, hydrogen):
        return str(hydrogen.GetNeighbors()[0].GetSymbol())
        # another lines of code may be added, depends on search of aromatic substructures

    def complete_classification(self, atom_types, molecule, counter):
        for atm_idx, atom_type in enumerate(atom_types):
            if atom_type is None:
                if molecule.GetAtomWithIdx(atm_idx).GetAtomicNum() == 1:
                    atom_type = self.classifify_remaining_H(molecule.GetAtomWithIdx(atm_idx))
                else:
                    atom_type = 'plain'
            atom_type = molecule.GetAtomWithIdx(atm_idx).GetSymbol() + '*' + atom_type  # str(atm_idx) + '|' +
            atom_types[atm_idx] = atom_type
            counter[atom_type] += 1

    def analyze_aromatic_rings(self, molecule, atom_types, counter):
        rings = molecule.GetRingInfo().AtomRings()
        for ring in rings:
            for atm_idx in ring:
                atom = molecule.GetAtomWithIdx(atm_idx)
                # print(atom.GetSymbol(), atm_idx)
                if atom.GetIsAromatic() and atom_types[atm_idx] is None:
                    atom_types[atm_idx] = 'A'
                    # implicit and explicit Hs not detected, necessary to use brute force
                    # methods Atom.GetNumExplicitHs(), .GetNumImplicitHs(), atom.GetTotalNumHs()) return 0 for each atom
                    for neigh in atom.GetNeighbors():
                        # aromatic Hs detection
                        if neigh.GetAtomicNum() == 1:
                            atom_types[neigh.GetIdx()] = f'A{atom.GetSymbol()}'

    def get_substructures(self, mol_idx, molecule, counter):
        # list where the atom type labels will be put to
        assigned_atom_types = [None] * molecule.GetNumAtoms()
        make_file = False

        self.analyze_aromatic_rings(molecule, assigned_atom_types, counter)
        for pattern, atom_types in self.SMARTS_and_atom_types:
            if molecule.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                pattern_atoms = molecule.GetSubstructMatches(Chem.MolFromSmarts(pattern))
                for atom_tuple in pattern_atoms:
                    for atm_idx, atom_type in zip(atom_tuple, atom_types):
                        if assigned_atom_types[atm_idx] is None:
                            assigned_atom_types[atm_idx] = atom_type
                # just for my own use
                # if not make_file:
                #     make_file = True
                # a = molecule.GetAtomWithIdx(pattern_atoms[0][0])
                # b = molecule.GetAtomWithIdx(pattern_atoms[0][2])
                # c = molecule.GetAtomWithIdx(pattern_atoms[0][4])
                # counter[a.GetSymbol()] += 1
                # counter[b.GetSymbol()] += 1
                # counter[c.GetSymbol()] += 1
                # print(a.GetSymbol())   # , b.GetSymbol(), c.GetSymbol()

                # just for my own use
                # if make_file:
                #     create_substruct_outputs(self.sdf_name, mol_idx, molecule)

            # move completion one level higher - complete classification at 'molecule set' level ? better ?
        self.complete_classification(assigned_atom_types, molecule, counter)
        return assigned_atom_types
        # !!! get molecule name: mol.GetProp('_Name') -> returns NSC_1000089, NSC_100992 etc.

    # former version, before "removeHS" in supplier was set to False:
    # mol2 = Chem.AddHs(mol, addCoords=True)

    def classify_atoms(self):
        # create_mol_from_pattern('C[CX3](=[OX1])[OX2]C', 'outputs/testmol_modified_ester.sdf')
        counter = Counter()
        for i, mol in enumerate(self.supplier):
            # print('Molecule no.', i)
                self.all_atom_types.append(
                    self.get_substructures(i, mol, counter)
                )
        print(counter)
