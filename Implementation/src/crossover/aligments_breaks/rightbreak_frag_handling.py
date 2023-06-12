import copy

import rdkit
from rdkit import Chem

from src.crossover.aligments_breaks.find_biggest_frag import find_biggest_frag

rdkit.RDLogger.DisableLog("rdApp.*")

def ringbreak_frag_handling(new_core, mcs_ringbreak_idx):
    """
    This takes a rdkit mol of the core (after atoms were removed to resolve
    cyclic breaks). It then tests and handles ringbreaks and fragmentation.

    if an atom in the core was deleted, if it had additional atoms attached to
    it it can cause fragmenation. so this step handles that and takes the
    largest fragment as the new common core.

    Inputs:
    :param rdkit.Chem.rdchem.Mol new_core: common core mol object
    :param list mcs_ringbreak_idx: list of the idx's of the common core for
        iso labels and later adjustment iso labels and idx numbers

    Returns:
    :returns: list iso_core_frag_list: list of the idx's of common core; same
        as mcs_ringbreak_idx unless there was fragmentation that needed to be
        handled.
    """

    ##########################
    # Check for fragmentation, if core is fragmented than additional
    # processing is required to find the largest frag and to then reassign the
    # indexes in the ligs and core

    check_fragmentation = Chem.GetMolFrags(new_core, asMols=True, sanitizeFrags=False)
    num_frag_len = len(check_fragmentation)
    iso_core_frag_list = copy.deepcopy(mcs_ringbreak_idx)

    if num_frag_len > 1:
        # determine the largest fragment in the list of frags
        largest_frag, largest_frag_index_num = find_biggest_frag(check_fragmentation)

        # make a list without the largest fragment
        list_frag_mols = []
        list_of_frag_idxs = range(0, len(check_fragmentation))

        for i in list_of_frag_idxs:
            if i == largest_frag_index_num:
                continue

            frag = check_fragmentation[int(i)]
            list_frag_mols.append(frag)

        # get the idx for all atoms in all frags EXCEPT THE LARGEST FRAG.
        # these will be the idx's of the original common core, before deleting
        # things which will be identified using the Isolabels we added before.
        # We will be deleting these atoms shortly
        for frag in list_frag_mols:

            # get all atom idx's (for the original unaltered common_core) in
            # the frag based on the Iso-labels we added before
            for atoms in frag.GetAtoms():
                index_val = atoms.GetIsotope() - 10000
                iso_core_frag_list.append(index_val)

        # Remove redundancy
        iso_core_frag_list = list(set(iso_core_frag_list))

        return iso_core_frag_list

    # if no fragmentation occured
    return iso_core_frag_list