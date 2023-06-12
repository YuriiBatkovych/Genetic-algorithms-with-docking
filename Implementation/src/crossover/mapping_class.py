"""
This script holds the Mapping class.
This is used when mapping most common substructure (MCS)
    to combine two molecules.
"""
import __future__

import random
import copy


class Mapping(object):

    def __init__(self, b_to_is, i_to_bs):
        self.b_to_is = copy.deepcopy(b_to_is)
        self.i_to_bs = copy.deepcopy(i_to_bs)

    def locate_b(self, i):
        return self.i_to_bs[i]

    def locate_i(self, b):
        return self.b_to_is[b]

    def delete_b(self, b):
        i_list_to_modify = self.locate_i(b)
        for i in i_list_to_modify:
            blank = self.i_to_bs[i].remove(b)
        del self.b_to_is[b]

    def delete_i(self, i):
        bs_to_modify = self.locate_b(i)
        for b in bs_to_modify:
            self.b_to_is[b].remove(i)
        del self.i_to_bs[i]

    def chose_b_from_i(self, i):
        if i in list(self.i_to_bs.keys()):
            options = self.locate_b(i)
            if len(options) > 1:
                b_x = random.choice(options)
            elif len(options) == 1:
                b_x = options[0]
            else:
                return "None"

            list_is = self.locate_i(b_x)
            list_bs = []
            for x in list_is:
                list_bs.append(self.locate_b(x))

            flattened = [val for sublist in list_bs for val in sublist]
            unique_bs = list(set(flattened))
            for b in unique_bs:
                self.delete_b(b)

            for x in list_is:
                self.delete_i(x)
            return b_x

        return "None"

    def testing_function_return_self_dicts(self):
        return self.b_to_is, self.i_to_bs


def run_mapping(b_dict, i_dict):
    a_mapping_object = Mapping(b_dict, i_dict)
    bs_chosen = []
    for i in i_dict:
        b_choice = a_mapping_object.chose_b_from_i(i)
        bs_chosen.append(b_choice)

    bs_chosen = list(set(bs_chosen))

    for i in bs_chosen:
        if i == "None":
            bs_chosen.remove(i)

    return bs_chosen
