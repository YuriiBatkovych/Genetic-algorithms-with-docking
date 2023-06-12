def make_dict_all_atoms_iso_to_idx_dict(mol):
    """
    Make a dictionary of every atom in a molecule with Iso as the key and the
    Idx as its value.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: an rdkit molecule

    Return
    :returns: dict mol_iso_to_idx_dict: a dictionary of the iso-label of every
        atom in the mol as the keys and the idx of that atom in the mol object.
        ie) {1008: 7, 1009: 8, 1003: 4, 1004: 3, 1010: 9, 1006: 5, 1007: 6, 10000:
        0, 10001: 1, 10002: 2, 1005: 10}
    """

    mol_iso_to_idx_dict = {}
    for atom in mol.GetAtoms():
        iso = atom.GetIsotope()
        idx = atom.GetIdx()
        mol_iso_to_idx_dict[iso] = idx
    return mol_iso_to_idx_dict


