from rdkit import Chem


def create_mol_from_pattern(pattern, file_path):
    mol = Chem.MolFromSmarts(pattern)
    # without UpdatePropertyCache() program falled. Source: https://github.com/rdkit/rdkit/issues/1596
    mol.UpdatePropertyCache()
    mol2 = Chem.AddHs(mol, addCoords=True)
    writer = Chem.SDWriter(file_path)
    writer.write(mol2)