from ..classifier import Classifier
from rdkit import Chem
from collections import Counter
import pprint
from ..io import load_external_atom_types
from ..io import create_output_paths, create_substruct_outputs


class SubstructClassifier(Classifier):
    def __init__(self):
        super().__init__(__file__)
        self.SMARTS_and_atom_types = load_external_atom_types('SMARTS')

    def classifify_remaining_H(self, hydrogen):
        return '-' + str(hydrogen.GetNeighbors()[0].GetSymbol())

    def complete_classification(self, atom_types, molecule):
        for atm_idx, atom_type in enumerate(atom_types):
            if atom_type is None:
                if molecule.GetAtomWithIdx(atm_idx).GetAtomicNum() == 1:
                    atom_type = self.classifify_remaining_H(molecule.GetAtomWithIdx(atm_idx))
                else:
                    atom_type = 'plain'
            atom_type = molecule.GetAtomWithIdx(atm_idx).GetSymbol() + '*' + atom_type
            atom_types[atm_idx] = atom_type

    def analyze_aromatic_rings(self, molecule, atom_types):
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
                            atom_types[neigh.GetIdx()] = f'-{atom.GetSymbol()}A'

    def get_atom_types(self, molecule):
        # list where the atom type labels will be put to
        mol_atom_types = [None] * molecule.GetNumAtoms()
        make_file = False
        self.analyze_aromatic_rings(molecule, mol_atom_types)
        for pattern, loaded_atom_types in self.SMARTS_and_atom_types:
            if molecule.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                pattern_atoms = molecule.GetSubstructMatches(Chem.MolFromSmarts(pattern))
                for atom_tuple in pattern_atoms:
                    for atm_idx, atom_type in zip(atom_tuple, loaded_atom_types):
                        if mol_atom_types[atm_idx] is None:
                            mol_atom_types[atm_idx] = atom_type
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

            # move completion higher - complete classification at 'molecule set' level ? better ?
        self.complete_classification(mol_atom_types, molecule)
        return mol_atom_types



