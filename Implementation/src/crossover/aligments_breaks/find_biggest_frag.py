def find_biggest_frag(frag_mols_obj):
    """
    This will take a frag mol object and return the largest fragment and the
    index in frag_mols_obj.

    Inputs:
    :param tuple frag_mols_obj: A tuple containing all the fragments of an
        rdkit mol.

    Returns:
    :returns: rdkit.Chem.rdchem.Mol frag_mols_obj: The largest rdkit mol obj
        in the provided tuple
    :returns: int idx_of_max: the idx number of the largest rdkit mol obj in
        the provided tuple.
    """

    if len(frag_mols_obj) > 1:

        idx_of_max = None
        num_atoms_max = None

        for i in range(0, len(frag_mols_obj)):
            frag = frag_mols_obj[i]
            atom_count = frag.GetNumAtoms()
            if num_atoms_max is None:
                idx_of_max = i
                num_atoms_max = atom_count
            elif num_atoms_max < atom_count:
                idx_of_max = i
                num_atoms_max = atom_count
            else:
                continue

        return frag_mols_obj, idx_of_max

    else:
        return frag_mols_obj, 0