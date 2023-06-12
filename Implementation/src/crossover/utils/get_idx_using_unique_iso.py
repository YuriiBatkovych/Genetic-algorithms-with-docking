def get_idx_using_unique_iso(mol, iso_val):
    for atom in mol.GetAtoms():
        if atom.GetIsotope() == iso_val:
            idx = atom.GetIdx()
            return idx
    return None