import copy

import rdkit

from src.crossover.aligments_breaks.rightbreak_frag_handling import ringbreak_frag_handling
from src.crossover.utils.remove_iso_labels import remove_iso_labels
from src.gypsum_dl.gypsum_dl.MolObjectHandling import check_sanitization, remove_atoms

rdkit.RDLogger.DisableLog("rdApp.*")

def check_cyclic_breaks(alignment_tuple, mol_1, mol_2, core):
    """
    Check for cyclic breaks and fixes them.

    Fixing cyclic breaks is handled by removing the atoms from the core (MCS
    substructure). Removing atoms will often cause fragmentation (ie. removing
    a carbon from a methyl will leave 3 fragmented hydrogens). Fragmented
    atoms need to be removed. The largest fragment is assumed to be the core
    of interest which we don't want to remove.

    Inputs:
    :param tuple alignment_tuple: a tuple with the atoms with IDX which match
        in the same order ie. alignment_tuple[0][0] is the IDX for the atom in
        mol_1 which is matched to alignment_tuple[1][0] and alignment_tuple[2][0]
        (for mol_2) alignment_tuple[3] (for core) is the reference index which
        ranges from 0 to the number of atoms in the order (0,1,2,3,4...n)
    :param rdkit.Chem.rdchem.Mol mol_1: rdkit mol for ligand 1
    :param rdkit.Chem.rdchem.Mol mol_2: rdkit mol for ligand 2
    :param rdkit.Chem.rdchem.Mol core: rdkit mol for shared common core
        between mol_1 and mol_2

    Returns:
    :returns: rdkit.Chem.rdchem.Mol core: the original core rdkit mol returned
        if did_a_ring_break is False
    :returns: tuple alignment_tuple: the unaltered input param alignment_tuple
        returned if did_a_ring_break is False
    :returns: bool did_a_ring_break: True if the ring broke and was fixed;
        False if there were no breaks and required no modifications to be made to
        the alignment_tuple or core
    :returns: rdkit.Chem.rdchem.Mol new_core: the modified core rdkit mol
        returned if did_a_ring_break is True
    :returns: tuple new_align_tuple: the modified alignment_tuple returned if
        did_a_ring_break is True
    :returns: bool None: returns 3 Nones if it failed to fix the cyclic breaks
    """

    if type(mol_1) is not rdkit.Chem.rdchem.Mol:
        return None, None, None
    if type(mol_2) is not rdkit.Chem.rdchem.Mol:
        return None, None, None
    if type(core) is not rdkit.Chem.rdchem.Mol:
        return None, None, None

    mcs_ringbreak_idx = []
    mol_1_ringbreak_idx = []
    mol_2_ringbreak_idx = []
    for l1, l2, c1 in zip(alignment_tuple[0], alignment_tuple[1], alignment_tuple[2]):
        atom1 = mol_1.GetAtomWithIdx(l1)
        atom2 = mol_2.GetAtomWithIdx(l2)
        atom_c = core.GetAtomWithIdx(c1)

        # ring breaks can occur when an atom in either lig is a ring atom
        # but the common substructure has that as a non-ring atom
        if atom_c.IsInRing() is False and (
                atom1.IsInRing() is True or atom2.IsInRing() is True
        ):
            mcs_ringbreak_idx.append(l1)
            mol_1_ringbreak_idx.append(l2)
            mol_2_ringbreak_idx.append(c1)

    if len(mcs_ringbreak_idx) > 0:
        new_align_list_l1 = []
        new_align_list_l2 = []
        new_align_list_c1 = []

        temp_core_removed_breaks = copy.deepcopy(core)

        temp_core_removed_breaks = remove_atoms(temp_core_removed_breaks, mcs_ringbreak_idx)
        if temp_core_removed_breaks is None:
            return None, None, None

        all_atoms_to_delete = ringbreak_frag_handling(
            temp_core_removed_breaks, mcs_ringbreak_idx
        )
        if all_atoms_to_delete is None:
            return None, None, None

        # Now work on the original core. THIS WILL BE THE OFFICIAL NEW CORE.
        # delete any cyclic breaks or anything connected to a cyclic break
        # which would fragment only delete from core mol
        new_core = remove_atoms(core, all_atoms_to_delete)
        if new_core is None:
            return None, None, None
        new_core = check_sanitization(new_core)
        if new_core is None:
            return None, None, None

        # now that we've made a new core, the idx's are different so we need
        # to relabel mol_1 and mol_2

        # remove the Iso-labels from lig 1 and 2 for anything deleted
        remove_iso_labels(mol_1, all_atoms_to_delete)
        remove_iso_labels(mol_2, all_atoms_to_delete)

        # make a new index series for comparing mol_1 and mol_2 to the core.
        # this is done using the original indexing and the atoms which were
        # removed from mcs.
        count = 0
        for l1, l2, c1 in zip(
                alignment_tuple[0], alignment_tuple[1], alignment_tuple[2]
        ):
            if c1 not in all_atoms_to_delete:
                new_align_list_l1.append(l1)
                new_align_list_l2.append(l2)
                new_align_list_c1.append(count)
                count = count + 1
        new_align_tuple = (new_align_list_l1, new_align_list_l2, new_align_list_c1)
        did_a_ring_break = True
        return new_core, new_align_tuple, did_a_ring_break

    # len(mcs_ringbreak_idx) less than or equal to 0
    did_a_ring_break = False
    return core, alignment_tuple, did_a_ring_break