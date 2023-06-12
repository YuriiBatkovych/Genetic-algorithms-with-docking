from src.crossover.groups_dict_r.r_group_list import r_group_list


def check_replace_mol(mol_1, mol_2, mcs_mol):
    """
    Confirm that mcs_mol can be replaced in mol_1 and mol_2 around 0.8% of the
    time this function fails so we will filter this 1st

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_1: an rdkit mol
    :param rdkit.Chem.rdchem.Mol mol_2: an rdkit mol
    :param rdkit.Chem.rdchem.Mol mcs_mol: rdkit mol for shared common core
        between mol_1 and mol_2

    Returns:
    :returns: bool True/False: Returns True if it passes for both mol_1 and
        mol_2 returns False if either fails.
    """

    temp = r_group_list(mol_1, mcs_mol)
    if temp is None:
        return False
    temp = r_group_list(mol_2, mcs_mol)
    if temp is None:
        return False
    return True