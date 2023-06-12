def add_mcs_isolabels(mol_1, mol_2, common_core, picked_alignment):
    """
    This will modify every atom in mol_1, mol_2, and the common_core to have
    the same isotope labels based on the index of the common_core atoms.

    Isotope number is set as the index number for the common_core atoms + 10,000

    Input
    :param rdkit.Chem.rdchem.Mol mol_1: an rdkit molecule
    :param rdkit.Chem.rdchem.Mol mol_2: an rdkit molecule
    :param Chem.MolFromSmarts(mcs_res_SMART) common_core: mcs_res_Mol
    :param tupple picked_alignment: a tuple with 3 subtuples for
        mol_1,mol_2,common_core. The numbers within the sublist is the atom IDx
        for a given atom in a ligand. The sublist index for each atom in for
        picked_alignment corresponds to the Idx of that atoms match in the
        commmon_core. ie picked_alignment[1][3] = 10; thus the IDx atom in mol_2
        which corresponds to the 3rd atom in the common_core is 10.

    Returns:
    :returns: tuple final_index: tuple with three sublists which are the same
        as picked_alignment[2]).
    """
    i = 0
    index_list = []
    for lig1, lig2, c1 in zip(
            picked_alignment[0], picked_alignment[1], picked_alignment[2]
    ):
        atom1 = mol_1.GetAtomWithIdx(lig1)
        atom2 = mol_2.GetAtomWithIdx(lig2)
        atom_c = common_core.GetAtomWithIdx(c1)

        atom1.SetIsotope(10000 + i)
        atom2.SetIsotope(10000 + i)
        atom_c.SetIsotope(10000 + i)
        index_list.append(i)
        i = i + 1
    final_index = index_list, index_list, index_list
    return final_index
