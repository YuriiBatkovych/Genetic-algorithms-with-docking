def make_b_dic(i_dictionary, r_dict_num, lig_number):
    """
    This generates the dictionaries for the B-groups. one is to track the
    R-groups which a B-group represents (this is the b_to_r_master_dict). one
    is to track the anchor atoms a B-group branches from (this is the
    b_to_anchor_master_dict).

    Inputs:
    :param dict i_dictionary:dictionary for R groups bound to nodes (aka I's).
        ie) {'10008':[1R1,1R2],'10009':[1R2,1R3]}
    :param dict r_dict_num: dictionary for anchors which are attached to an R
        group. ie) {'1R1':[10008],'1R2':[10008,10009],'1R3':[10009]}
    :param int lig_number: an int either 1 or 2 for (mol_1 or mol_2
        respectively)

    Returns:
    :returns: dict b_to_r_master_dict: key is unique B-name and the R-groups
        it represents. example {'1B1':['1R1'],'1B2':['1R2','1R3','1R4'],'1B3':
        ['1R5']}
    :returns: dict b_to_anchor_master_dict: key is unique B-name and items are
        anchors that B connects to. example
        {'1B1':[10008,10007],'1B2':[10000],'1B3':[10006]}
    """

    k = lig_number
    b_to_r_master_dict = {}
    b_to_anchor_master_dict = {}
    counter = 1
    anchor_list = list(i_dictionary.keys())

    while len(anchor_list) > 0:
        anchor = anchor_list[0]
        B_key = "{}B{}".format(k, counter)
        temp_r_list = []
        temp_anchor_list = []

        for Rs in i_dictionary[anchor]:
            temp_r_list.append(Rs)
            r_dict_i = r_dict_num[Rs]
            for I in r_dict_i:
                temp_anchor_list.append(I)

        temp_anchor_list = list(set(temp_anchor_list))
        temp_r_list = list(set(temp_r_list))

        # make new B-group entry in the dictionaries
        b_to_r_master_dict[B_key] = temp_r_list  # This B-represents these R-groups
        b_to_anchor_master_dict[
            B_key
        ] = temp_anchor_list  # This B connects to these anchor atoms

        counter = counter + 1

        # make a list of atoms to remove in the next iteration if they are in
        # both the temp_anchor_list and anchor_list
        for i in temp_anchor_list:
            if i in anchor_list:
                anchor_list.remove(i)

    return b_to_r_master_dict, b_to_anchor_master_dict
