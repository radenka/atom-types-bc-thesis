from rdkit import Chem

def create_mol_from_pattern(pattern, file_path):
    mol = Chem.MolFromSmiles(pattern)
    mol2 = Chem.AddHs(mol, addCoords=True)
    writer = Chem.SDWriter(file_path)
    writer.write(mol2)