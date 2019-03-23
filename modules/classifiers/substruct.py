from modules.classifier import Classifier
from rdkit import Chem
from ..file import create_output_paths, create_substruct_outputs
from ..utils import create_mol_from_pattern


class SubstructClassifier(Classifier):
    def __init__(self, input_sdf):
        super().__init__(input_sdf)

    def get_substructure(self, mol_idx, molecule):
        # atom_types = [None] * molecule.GetNumAtoms()
        make_file = False
        # print(mol_idx, Chem.MolToSmiles(molecule))

        # list where the atom type labels will be put to
        atom_types = []

        # import patterns from external txt file
        # SMARTS for carboxylic acid gives the same results as appropriate SMILES but returns the first carbon only
        # carboxylic acid '[$([$([C;$(C=[$([O;D1;$(O=C)])])]);$(C[$([O;$([H1&-0,H0&-1])])]);$(C[#6,#1])])]'

        # carboxylic acid ver. 2 [CX3](=O)[OX2H1] doesn't detect Hs, returns tuples of 3 atoms
        # charged carboxylic acid [CX3](=O)[OX1-] results found
        # oxygen anion O- [OX1-] result s found, charged NO2 groups included
        # ether [#6]O[#6] (includes esters, too)
        # disulfide [#6]SS[#6]
        # aldehyde C(=O)[H] detects Hs, too!
        # carbonic ester C[CX3](=[OX1])[OX2]C
        #    or [CX3](=[OX1])[OX2]C -> returns substances composed of substructures as NC(=O)OC, more result than above
        # nitrile group [NX1]#[CX2]
        # nitroso group [NX2]=[OX1]; N(=O) is detected also in nitro groups - SMARTS is more precise
        # nitro group:x
        #   no outputs for [NX3](=O)=O
        #   -> use [NX3+](=O)[O-] only
        # anhydride [CX3](=[OX1])[OX2][CX3](=[OX1])

        # OrderedDict?


        patterns = ['[CX3](=O)[OX1-]']    # C(=O)[OH]', 'C(=O)', 'C[OH]', 'C(=O)[H]'
        all_pattern_atoms = {}
        for pattern in patterns:
            pattern_atoms = molecule.GetSubstructMatches(Chem.MolFromSmarts(pattern))
            # if mol_idx == 460:
            if len(pattern_atoms) > 0:
                print('Molecule no.', mol_idx)
                # just for my use
                all_pattern_atoms[pattern] = str(pattern_atoms)
                # create_output_paths(self.sdf_name)
                if not make_file:
                    make_file = True

                print(pattern, pattern_atoms)
        # just for my own use
        # if make_file:
        #     create_substruct_outputs(self.sdf_name, mol_idx, molecule, all_pattern_atoms)
        # !!! get molecule name: mol.GetProp('_Name') -> returns NSC_1000089, NSC_100992 etc.

    # former version, before "removeHS" in supplier was set to False:
    # mol2 = Chem.AddHs(mol, addCoords=True)


    def classify_atoms(self):
        # create_mol_from_pattern('C[CX3](=[OX1])[OX2]C', 'outputs/testmol_modified_ester.sdf')

        for i, mol in enumerate(self.supplier):
            # print('Molecule no.', i)
            self.all_atom_types.append(
                self.get_substructure(i, mol)
            )