def get_atoms_touch_mcs(mol):
    """
    Function to find all neighbors for a set of molecules touching. Isolabeled
    core atoms.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: isolabeled with atoms in the core having
        isotope. labels set as their idx number + 10000 and atoms not shared in
        the common core isotope labels set as:
            for lig_1: atom idx number + 1000
            for lig_1: atom idx number + 2000


    Returns:
    :returns: dict mcs_touches dict:  a dictionary with keys being the isotope
        label of core atoms and the items being the idx's of all non-core atoms
        which touch it. If a core atom touch no non-core atoms it will not be
        added to the dictionary.
    """

    mcs_touches = {}
    all_atoms = mol.GetAtoms()

    for atom in all_atoms:
        iso = atom.GetIsotope()
        if iso > 9999:
            # then its a core atom
            neighbors = atom.GetNeighbors()
            values = []

            for neighbor_atom in neighbors:
                iso_neighbor_x = neighbor_atom.GetIsotope()
                if iso_neighbor_x < 9999:
                    idx_of_neighbor = neighbor_atom.GetIdx()
                    values.append(idx_of_neighbor)
                    mcs_touches[iso] = values

    return mcs_touches