import copy

from rdkit import Chem

from src.crossover.core.make_anchor_to_bonds_and_type_for_frag import make_anchor_to_bonds_and_type_for_frag
from src.crossover.core.make_dict_all_atoms_iso_to_idx_dict import make_dict_all_atoms_iso_to_idx_dict
from src.crossover.core.unpack_list_of_atoms_and_bond_types import unpack_lists_of_atoms_and_bond_type
from src.gypsum_dl.gypsum_dl.MolObjectHandling import check_sanitization, remove_atoms


def merge_smiles_with_core(rs_chosen_smiles, mcs_mol):
    """
    This function runs most of the ligand merger portion of the script. It
    merges the chosen R-groups (rs_chosen_smiles) with the common core
    (mcs_mol) by bonding anchors in the core to the atoms bound to the
    respective anchor in the R-group fragment.

    Variables used in this function
        anchor_to_connection_dict: dict of anchor to connected atoms and bond
        types.

        example anchor_to_connection_dict = {10008: [1016, rdkit.Chem.rdchem.BondType.AROMATIC],
                                            10007: [1013, rdkit.Chem.rdchem.BondType.AROMATIC]})

        mol_frag, keys are isotope label, value is Idx.

        example mol_frag_iso_to_idx_dict = {10007: 0, 10008: 8, 1013: 1, 1014: 3, 1015: 5, 1016: 7,
                                           1017: 9, 1018: 6, 1019: 4, 1020: 2})

    Inputs:
    :param list rs_chosen_smiles: A list containing the  SMILES strings for
        the chosen R-groups to add
    :param rdkit.Chem.rdchem.Mol mcs_mol: an rdkit molecule representing the
        Most common substructure (MCS) which will be expanded by adding R-groups
        to make the child molecule

    Returns:
    :returns: rdkit.Chem.rdchem.Mol rw_core_merg: The child molecule with the
        added R-groups built onto the mcs_mol. returns None if the process fails
        or if a None-type makes it through.
    """

    # convert to RWMOL class of molecule which are able to add and remove
    # bonds. RWMOL class is the Read and Write-Mol Class in rdkit.
    rw_core_merg = Chem.RWMol(mcs_mol)

    # sanitize the mol_frag
    rw_core_merg = check_sanitization(rw_core_merg)
    if rw_core_merg is None:
        # ("rw_core_merg failed to be sanitizable (merge_smiles_with_core)")
        return None

    for r_groups in rs_chosen_smiles:
        for frag in r_groups:
            # make a rdkit mol out of the smiles string of the R-group frag
            mol_frag = Chem.MolFromSmiles(frag, sanitize=False)

            mol_frag_copy = copy.deepcopy(mol_frag)
            mol_frag = check_sanitization(mol_frag)
            if mol_frag is None:
                mol_frag = mol_frag_copy

            anchor_to_connection_dict = make_anchor_to_bonds_and_type_for_frag(mol_frag)
            mol_frag_iso_to_idx_dict = make_dict_all_atoms_iso_to_idx_dict(mol_frag)

            anchors_idxs_to_remove = []
            for anchors in list(anchor_to_connection_dict.keys()):
                # anchors are iso numbers of 10,000 or higher
                idx_val = mol_frag_iso_to_idx_dict[anchors]
                anchors_idxs_to_remove.append(idx_val)

            # remove the anchor atoms from mol_frag
            mol_frag = remove_atoms(mol_frag, anchors_idxs_to_remove)
            rw_core_merg = Chem.CombineMols(rw_core_merg, mol_frag)

            # convert to RWMOL class of molecule which are able to add and
            # remove bonds
            rw_core_merg = Chem.RWMol(rw_core_merg)

            # make a dictionary of every atom in rw_core_merg with Iso as the
            # key and the Idx as its value
            core_merg_iso_to_idx_dict = make_dict_all_atoms_iso_to_idx_dict(rw_core_merg)

            for anchor_atom_iso in list(anchor_to_connection_dict.keys()):
                # Idx of the anchor in merged core
                idx_for_anchor = core_merg_iso_to_idx_dict[anchor_atom_iso]
                list_of_atom_idx, list_of_bond_types = unpack_lists_of_atoms_and_bond_type(
                    anchor_to_connection_dict, anchor_atom_iso, core_merg_iso_to_idx_dict
                )

                for ai, bt in zip(list_of_atom_idx, list_of_bond_types):
                    try:
                        rw_core_merg.AddBond(idx_for_anchor, ai, bt)
                    except:
                        return None

    return rw_core_merg
