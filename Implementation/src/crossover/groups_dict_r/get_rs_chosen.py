def get_rs_chosen_smiles(rs_chosen, r_smiles_dict_1, r_smiles_dict_2):
    """
    This function returns a list of SMILES strings for every R-group chosen.
    It requires the R_smile_dictionary for both ligands to function.

    Inputs:
    :param list rs_chosen: A list of the chosen R-groups which will be used to
        generate a new mol. ie) ['2R2', '1R1']
    :param dict r_smiles_dict_1: A dictionary which has can find the SMILES
        string for each R-group of Ligand 1. ie) {'1R1':
        '[10006*][1009N]=[1008N+]=[1007N-]'}
    :param dict r_smiles_dict_2: A dictionary which has can find the SMILES
        string for each R-group of Ligand 2. ie) {'2R2': '[10006*][2009OH]',
        '2R1': '[10003*][2007CH2][2008OH]'}


    Returns:
    :returns: list rs_chosen_smiles: A list of all the SMILES string which are
        to be added to make the child ligand. Each SMILES is a sublist.
        ie)[['[10006*][1009N]=[1008N+]=[1007N-]'],['[10006*][2009OH]']]
    """

    rs_chosen_smiles = []
    for R in rs_chosen:
        Rs_for_the_R = []
        lig_number = R[0]
        if lig_number == str(1):
            Rs_for_the_R.append(r_smiles_dict_1[R])
        elif lig_number == str(2):
            Rs_for_the_R.append(r_smiles_dict_2[R])

        rs_chosen_smiles.append(Rs_for_the_R)

    return rs_chosen_smiles


def get_rs_chosen_from_bs(bs_chosen, b_to_r_master_dict_1, b_to_r_master_dict_2):
    """
    this function returns a list of R-groups chosen based on the list of
    chosen B's. It requires the b_to_r_master_dict_1 for both ligands to
    function.

    Inputs:
    :param list bs_chosen: A list of the chosen B-groups. ie) ['1B1', 1B2',
        '2B3']
    :param dict b_to_r_master_dict_1: a Dictionary to reference B and R-groups
        from mol_1. keys are names of B-groups; items are R-groups that a B-group
        represents. ie) {'1B1':['1R1'],'1B2':['1R2','1R3','1R4'],'1B3': ['1R5']}
    :param dict b_to_r_master_dict_2: a Dictionary to reference B and R-groups
        from mol_2. keys are names of B-groups; items are R-groups that a B-group
        represents. ie) {'2B1':['2R1'],'2B2':['2R2','2R3','2R4'],'2B3':
        ['2R5','2R6]}

    Returns:
    :returns: list rs_chosen: a list containing all the R-groups represented
        by the chosen B-groups. ie) ['1R1', '1R2', '1R3','1R4', '2R5', '2R6']
    """

    rs_chosen = []
    for B in bs_chosen:
        Rs_for_the_B = []
        lig_number = B[0]
        if lig_number == str(1):
            for i in b_to_r_master_dict_1[B]:
                Rs_for_the_B.append(i)

        elif lig_number == str(2):
            for i in b_to_r_master_dict_2[B]:
                Rs_for_the_B.append(i)
        for i in Rs_for_the_B:
            rs_chosen.append(i)

    return rs_chosen