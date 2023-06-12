from rdkit import Chem


def r_groups_dict(mol_frags, lig_number_for_multiplier):
    """
    given a set of mol_frags and the ligand_number (ie. 1 for mol_1 and 2 for
    mol_2) this will make dictionaries of all the Rgroup and all the smiles
    for each Rgroup

    Input
    :param rdkit.Chem.rdchem.Mol mol_frags: a rdkit mol containing fragments
    :param int lig_number_for_multiplier: an int either 1 for mol_1 or 2 for
        mol_2, used to make labels which are traceable to the ligand being used

    Returns:
    :returns: dict r_chain_dictionary: a dictionary with the R-groups and the
        anchor atoms they connect to ie) {'1R1':[13,14],'1R2':[21,22],'1R3':[25]}
    :returns: dict r_smiles_dictionary: a dictionary with the R-groups and the
        SMILES strings of those groups ie
        {'1R1':'[1*]:[1013c]([1020H])[1014c]([1019H])[1015c]([1018H])[1016c](:[2*])[1017H]',
        '1R2':'[3*][1024C]([1026H])([1027H])[1023N] = [1022N+] = [1021N-]',
        '1R3':'[4*][1025O][1029H]'}
    """

    num_frags = len(mol_frags)
    i = 0
    r_chain_dictionary = {}
    r_smiles_dictionary = {}
    k = int(lig_number_for_multiplier)
    while i < num_frags:
        frag = mol_frags[i]
        r_list_temp = []
        r_list_smiles = Chem.MolToSmiles(frag, isomericSmiles=True)
        for atoms in frag.GetAtoms():
            iso = atoms.GetIsotope()
            if 3000 > iso > 100:
                r_list_temp.append(iso - (1000 * k))
                atoms.SetIsotope(0)
            if iso > 3000:
                r_list_temp.append(iso)
            lig_num_r_r_num = "{}R{}".format(k, i + 1)
            r_chain_dictionary[lig_num_r_r_num] = r_list_temp
            r_smiles_dictionary[lig_num_r_r_num] = r_list_smiles
        i = i + 1

    return r_chain_dictionary, r_smiles_dictionary
