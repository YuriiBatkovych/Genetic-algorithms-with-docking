import rdkit
from rdkit import Chem

from src.crossover.aligments_breaks.check_cyclic_breaks import check_cyclic_breaks
from src.crossover.utils.add_atom_isolabels import add_r_atom_isolabels
from src.crossover.utils.add_mcs_isolabels import add_mcs_isolabels
from src.crossover.aligments_breaks.pick_mcs_alignment import pick_mcs_alignment
from src.crossover.utils.renumber_to_mcs import renumber_to_mcs


def handle_mcs_align_labeling_and_cyclicbreaks(mol_1, mol_2, mcs_mol):
    """
    This will take 2 aligned molecules and pick a specific alignment, check
    for any ring and cyclic breaks and fragmentation removing any atoms from
    the common core which caused those issues, renumber and isotope atoms in
    mol_1, mol_2, mcs_mol to be tractable, and consistent amongst the three.

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_1: an rdkit molecule
    :param rdkit.Chem.rdchem.Mol mol_2: an rdkit molecule
    :param Chem.MolFromSmarts(mcs_res_SMART) mcs_mol: a result object from
        rdkits MCS function

    Returns:
    :returns: rdkit.Chem.rdchem.Mol mol_1: an rdkit molecule isolabeled and
        renumbered or None if it fails
    :returns: rdkit.Chem.rdchem.Mol mol_2: an rdkit molecule isolabeled  and
        renumbered or None if it fails
    :returns: mcs_res_Mol mcs_mol: an MCS result isolabeled and no breaks or
        None if it fails
    """

    if type(mol_1) is not Chem.rdchem.Mol:
        return None, None, None
    if type(mol_2) is not Chem.rdchem.Mol:
        return None, None, None
    if type(mcs_mol) is not Chem.rdchem.Mol:
        return None, None, None

    # Set Alignment, Isotope label, and Handle breaks
    picked_alignment = pick_mcs_alignment(mol_1, mol_2, mcs_mol)
    if picked_alignment is None:
        return None, None, None

    # Isotope label the MCS core
    index_tuple = add_mcs_isolabels(mol_1, mol_2, mcs_mol, picked_alignment)

    mol_1 = renumber_to_mcs(mol_1, picked_alignment[0])
    mol_2 = renumber_to_mcs(mol_2, picked_alignment[1])

    mcs_mol, new_index, are_there_breaks = check_cyclic_breaks(index_tuple, mol_1, mol_2, mcs_mol)

    if mcs_mol is None:
        return None, None, None

    if are_there_breaks is True:
        mcs_m3, new_index_2, are_there_breaks_2 = check_cyclic_breaks(new_index, mol_1, mol_2, mcs_mol)
        if (are_there_breaks_2 is True) or (mcs_mol is None):
            return None, None, None

    # Check for any additional fragmentation issues
    confirm_no_frag = Chem.GetMolFrags(mcs_mol, asMols=True, sanitizeFrags=False)
    if len(confirm_no_frag) != 1:
        return None, None, None

    mol_1 = renumber_to_mcs(mol_1, new_index[0])
    mol_2 = renumber_to_mcs(mol_2, new_index[1])

    add_r_atom_isolabels(mol_1, mol_2)
    return mol_1, mol_2, mcs_mol

