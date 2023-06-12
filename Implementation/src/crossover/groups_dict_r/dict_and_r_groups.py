"""
Dictionary and Dictionary handling functions
"""
import __future__

import rdkit
from rdkit import Chem

from src.crossover import mapping_class
from src.crossover.groups_dict_r.check_replace_mol import check_replace_mol
from src.crossover.groups_dict_r.get_atoms_touch_mcs import get_atoms_touch_mcs
from src.crossover.groups_dict_r.get_r_dict import get_r_dict
from src.crossover.groups_dict_r.get_rs_chosen import get_rs_chosen_from_bs, get_rs_chosen_smiles
from src.crossover.groups_dict_r.make_b_dict import make_b_dic
from src.crossover.groups_dict_r.r_group_list import r_group_list
from src.crossover.groups_dict_r.r_groups_dict import r_groups_dict
from src.crossover.groups_dict_r.replace_core_mol_dummy_atoms import replace_core_mol_dummy_atoms
from src.crossover.utils.inver_dict import invert_dictionary

rdkit.RDLogger.DisableLog("rdApp.*")

def handle_dicts_and_select_b_groups(mol_1, mol_2, mcs_mol):
    """
    this takes 3 rdkit.Chem.rdchem.Mol objects 1 for lig_1,lig_2, and the
    common core(mcs_mol). It creates all the necessary dictionaries, mapping,
    and selects the ligands that will be added to make the final molecule.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_1: rdkit mol for ligand 1
    :param rdkit.Chem.rdchem.Mol mol_2: rdkit mol for ligand 2
    :param rdkit.Chem.rdchem.Mol mcs_mol: rdkit mol for shared common core
        between mol_1 and mol_2

    Returns:
    :returns: list rs_chosen_smiles: smiles strings for the the R groups which
        correspond to the chosen B's returns None if it fails
    """

    they_pass = check_replace_mol(mol_1, mol_2, mcs_mol)
    if they_pass is False:
        return None

    r_smiles_dict_1, b_to_r_master_dict_1, b_to_anchor_master_dict_1 = \
        mol_handling_of_fragmenting_labeling_and_indexing(mol_1, mcs_mol, 1)

    if r_smiles_dict_1 is None:
        return None

    r_smiles_dict_2, b_to_r_master_dict_2, b_to_anchor_master_dict_2 = mol_handling_of_fragmenting_labeling_and_indexing(
        mol_2, mcs_mol, 2
    )

    if r_smiles_dict_2 is None:
        return None

    b_to_anchor_master = b_to_anchor_master_dict_1
    for i in list(b_to_anchor_master_dict_2.keys()):
        b_to_anchor_master[i] = b_to_anchor_master_dict_2[i]

    anchor_to_b_master = invert_dictionary(b_to_anchor_master)

    bs_chosen = mapping_class.run_mapping(b_to_anchor_master, anchor_to_b_master)

    rs_chosen = get_rs_chosen_from_bs(
        bs_chosen, b_to_r_master_dict_1, b_to_r_master_dict_2
    )

    rs_chosen_smiles = get_rs_chosen_smiles(rs_chosen, r_smiles_dict_1, r_smiles_dict_2)

    return rs_chosen_smiles


def mol_handling_of_fragmenting_labeling_and_indexing(mol, mcs_mol, lig_number):
    """
    This takes an rdkit mol for a ligand and 1 for the mcs_mol. It fragments
    the ligand by replacing the MCS. and it determines which anchors are in
    each fragment. These fragments are our R-groups and the assignment of
    anchors. is how we determine which R-group goes where relative to the MCS.

    lig_number  int    is the number of the ligand that is mol
                       ie if mol is mol_1 lig_number = 1
                       ie if mol is mol_2 lig_number = 2


    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: an rdkit mol (either mol_1 or mol_2)
    :param rdkit.Chem.rdchem.Mol mcs_mol: rdkit mol for shared common core
        between mol_1 and mol_2
    :param int lig_number: an int either 1 or 2 for (mol_1 or mol_2
        respectively)

    Returns:
    :returns: dict r_smiles_dictionary: a dictionary of the R-groups which
        branch off the common core keys are the R-groups; items are the SMILES
        strings of that R-groups returns None if it fails. Example: {'1R1':
        '[10003*][1007N]=[1013O]', '1R2': '[10000*][1011CH2]=[1008O]'}
    :returns: dict b_to_r_master_dict: A dictionary which tracks the R groups
        which belong to a B-group keys are the B-groups; items are the R-groups
        which belong to the B-group. returns None if it fails. Example: {'1B1':
        ['1R2'], '1B2': ['1R1']}
    :returns: dict b_to_anchor_master_dict: A dictionary which tracks the iso
        label of the anchor atoms for B-group. keys are the B-groups; items are
        the iso label of the anchor atoms for B-group returns None if it fails.
        Example:{'1B1': [10000], '1B2': [10003]}
    """

    mcs_touches = get_atoms_touch_mcs(mol)
    lig_r_atoms_touch_mcs = invert_dictionary(mcs_touches)

    replace_core = r_group_list(mol, mcs_mol)
    if replace_core is None:
        return None, None, None

    replace_core = replace_core_mol_dummy_atoms(replace_core)
    if replace_core is None:
        return None, None, None

    mol_frags = Chem.GetMolFrags(replace_core, asMols=True, sanitizeFrags=False)
    list_r_groups = []
    i = 0
    while i < len(mol_frags):
        val = Chem.MolToSmiles(mol_frags[i], isomericSmiles=True)
        list_r_groups.append(val)
        i = i + 1

    r_chain_dictionary, r_smiles_dictionary = r_groups_dict(mol_frags, lig_number)

    r_dict = get_r_dict(r_chain_dictionary, lig_r_atoms_touch_mcs)
    i_dict = invert_dictionary(r_dict)

    b_to_r_master_dict, b_to_anchor_master_dict = make_b_dic(i_dict, r_dict, lig_number)

    return r_smiles_dictionary, b_to_r_master_dict, b_to_anchor_master_dict




