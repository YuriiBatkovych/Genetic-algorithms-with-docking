def invert_dictionary(old_dic):
    """
    This will invert any dictionary so that the keys are the values and the
    values are the keys.

    Inputs:
    :param dict old_dic: a dictionary to invert

    Returns:
    :returns: dict inverted_dic: old_dict dict inverted so the keys are the
        items and the items are the keys
    """

    values = set([a for b in list(old_dic.values()) for a in b])
    values = list(values)
    inverted_dic = dict(
        (new_key, [key for key, value in list(old_dic.items()) if new_key in value])
        for new_key in values
    )

    return inverted_dic
