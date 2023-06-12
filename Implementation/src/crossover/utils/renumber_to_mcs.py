from rdkit import Chem


def renumber_to_mcs(mol, tuple_order_list):
    """
    This renumbers the indexes of the atoms in a lig to that of the MCS and
    returns a renumbered atom.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: an rdkit molecule
    :param tuple tuple_order_list: a tuple with 3 subtuples for
        mol_1,mol_2,common_core. The numbers within the sublist is the atom IDx
        for a given atom in a ligand. The sublist index for each atom in for
        tuple_order_list corresponds to the Idx of that atoms match in the
        commmon_core. ie tuple_order_list[1][3] = 10; thus the IDx atom in mol_2
        which corresponds to the 3rd atom in the common_core is 10.

    :returns:
    :returns: rdkit.Chem.rdchem.Mol mol: the same rdkit molecule but with the
        atom idx's renumbered to be consistent with the common core, as provided
        by the tuple_order_list
    """

    full_tuple_order = [x for x in tuple_order_list]
    num = 0
    while num < len(mol.GetAtoms()):
        if num in full_tuple_order:
            num = num + 1
        else:
            full_tuple_order.append(num)
            num = num + 1

    mol = Chem.RenumberAtoms(mol, full_tuple_order)
    return mol