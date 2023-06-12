def make_anchor_to_bonds_and_type_for_frag(mol_frag):
    """
    Create a dictionary w anchor atoms as keys.

    for each key, the items are broken into lists of lists with the 1st number
    of each as the isotope of the atom bound and the second value as the bond
    type to recreat bonds later to merge..

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_frag: an R-group which was converted into
        a mol

    Returns:
    :returns: dict anchor_to_connection_dict: a dictionary of anchor atom
        isolabels as keys, item is a list of the isolabel of the atom the key is
        bound to and the bond type. ie. anchor_to_connection_dict[10007] =
        [1004,Chem.BondType.AROMATIC]
    """

    anchor_to_connection_dict = {}
    isos_anchors_idxs_to_remove = []

    for atom in mol_frag.GetAtoms():
        if atom.GetIsotope() > 9999:  # if atom is an anchor
            iso_anchor = atom.GetIsotope()  # isotopes of the anchor atom
            anchor_idx = atom.GetIdx()  # get that atoms idx
            isos_anchors_idxs_to_remove.append(
                anchor_idx
            )

            # empty lists for subloop
            connection_iso_idx_list = (
                []
            )  # list of isotope number for atoms connected to an anchor
            bond_type_list = (
                []
            )  # list of bond types in the same order as connection_iso_idx_list

            neighbor = (
                atom.GetNeighbors()
            )  # all neighbor atoms to the Atom from above loop
            for x in neighbor:  # Atoms which are neighbors of anchor
                iso_neighbor_atom = x.GetIsotope()
                neighbor_bond_idx = x.GetIdx()
                connection_iso_idx_list.append(iso_neighbor_atom)

                # get bond type between anchor and connected atoms
                bond_object = mol_frag.GetBondBetweenAtoms(
                    anchor_idx, neighbor_bond_idx
                )
                bond_type = bond_object.GetBondType()
                bond_type_list.append(bond_type)

            for i, j in zip(connection_iso_idx_list, bond_type_list):
                list_of_atom_and_bond = [i, j]
                if iso_anchor in list(anchor_to_connection_dict.keys()):
                    tmp = anchor_to_connection_dict[iso_anchor]
                    tmp.append(list_of_atom_and_bond)
                    anchor_to_connection_dict[iso_anchor] = tmp
                else:
                    anchor_to_connection_dict[iso_anchor] = [list_of_atom_and_bond]

    return anchor_to_connection_dict
