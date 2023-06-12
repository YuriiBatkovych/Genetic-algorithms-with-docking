def unpack_lists_of_atoms_and_bond_type(anchor_to_connection_dict, anchor_atom_iso,
                                        core_merg_iso_to_idx_dict):
    """
    Iterate through all atoms which will be bound to the anchor and unpackage
    all the bond types in a list.

    Inputs:
    :param dict anchor_to_connection_dict: a dictionary of anchor isotope
        labels as keys and a lists as the items. these lists have 2 variables, the
        1st is the atom iso-label of the atom connected to an anchor, and the
        second variable is the rdkit bond type. ie) {10004: [1007,
        rdkit.Chem.rdchem.BondType.SINGLE]}
    :param int anchor_atom_iso: the interger of an anchor atom's isotope
        label. ie) 10004
    :param dict core_merg_iso_to_idx_dict: a dictionary of atom's isotope
        labels as keys and their corresponding Idx as the items. ie) {1008: 14,
        1014: 11, 1009: 15, 1010: 7, 1007: 13, 10000: 0, 10001: 1, 10002: 2,
        10003: 3, 10004: 4}

    Returns:
    :returns: list list_of_atom_idx: a list containing the atom idx. ie) [13]
    :returns: list list_of_bond_types: a list containing the bond types of
        bonds connected to an anchor atom. ie) [rdkit.Chem.rdchem.BondType.SINGLE]
    """

    list_of_atom_idx = []
    list_of_bond_types = []
    connection_list = anchor_to_connection_dict[anchor_atom_iso]

    if type(connection_list[0]) == int:
        # get the atom iso label
        atom_iso = anchor_to_connection_dict[anchor_atom_iso][0]

        # get the atoms idx in the merged core
        atom_idx = core_merg_iso_to_idx_dict[atom_iso]

        # get bond type
        bond_type = anchor_to_connection_dict[anchor_atom_iso][1]

        # append to the lists
        list_of_atom_idx.append(atom_idx)
        list_of_bond_types.append(bond_type)

    elif type(connection_list[0]) == list:
        if type(connection_list[0][0]) == int:

            for i in range(0, len(connection_list)):
                # get the atom iso label
                atom_iso = anchor_to_connection_dict[anchor_atom_iso][i][0]

                # get the atoms idx in the merged core
                atom_idx = core_merg_iso_to_idx_dict[atom_iso]

                # get bond type
                bond_type = anchor_to_connection_dict[anchor_atom_iso][i][1]

                # append to the lists
                list_of_atom_idx.append(atom_idx)
                list_of_bond_types.append(bond_type)

    return list_of_atom_idx, list_of_bond_types