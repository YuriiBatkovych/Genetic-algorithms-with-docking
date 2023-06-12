def get_r_dict(r_chain_dict, lig_r_atom_touch_mcs):
    """
    This will take the r_chain_dict and the dict of all the atoms which touch
    the core and return a dict of Rs groups as keys and their nodes as values

    Inputs:
    :param dict r_chain_dict: dict of all the atom isolabels for in an
        R-group. keys are R-groups;  items are iso-labels of atoms in the R-group.
        ie) {'1R1': [3, 4, 5, 6, 7, 8, 9, 10, 11, 10000]}
    :param dict lig_r_atom_touch_mcs: dict of all the atoms which directly
        touch the core and what anchor they touch. keys are atom isolabels of
        atoms touching an anchor; items are iso-labels of anchor atoms. ie) {3:
        [10000]}

    Returns:
    :returns: dict r_s_dict:  dictionary of R-groups and anchor atoms they are
        connected to. keys are R-groups. items are isolabel of anchor atoms. ie)
        {'1R1': [10000]}
    """

    r_s_dict = {}
    for key in list(r_chain_dict.keys()):
        node_list = []
        for atom in r_chain_dict[key]:
            for key_id in list(lig_r_atom_touch_mcs.keys()):
                if atom == key_id:
                    for x in lig_r_atom_touch_mcs[key_id]:
                        node_list.append(x)
                    r_s_dict[key] = node_list

    return r_s_dict