from rdkit import Chem
from rdkit.Chem import rdFMCS

from src.crossover.aligments_breaks.alignments_and_break_handling import handle_mcs_align_labeling_and_cyclicbreaks
from src.crossover.core.merge_smiles_with_core import merge_smiles_with_core
from src.crossover.groups_dict_r.dict_and_r_groups import handle_dicts_and_select_b_groups
from src.crossover.utils.handleHs import handleHs
from src.crossover.utils.remove_iso_labels import remove_all_isolabels
from src.gypsum_dl.gypsum_dl.MolObjectHandling import check_sanitization, handle_frag_check, check_for_unassigned_atom


def process_ligand_new_mol(ligand_new_mol):
    """
    This function processes the ligand_new_mol.
    It either returns the SMILES string of ligand_new_mol (ligand_new_smiles)
    or None if it failed at any step.

    Inputs:
    :param str lig_string_1: smile string for lig 1

    Returns:
    :returns: str ligand_new_smiles: either returns the SMILES
        string of ligand_new_mol or None if it failed at any step.
    """

    ligand_new_mol = check_sanitization(ligand_new_mol)
    if ligand_new_mol is None:
        return None

    # REMOVE ALL THE ISOTOPES IN THE NEW MOLECULE
    ligand_new_mol_final = remove_all_isolabels(ligand_new_mol)

    # Remove any fragments incase 1 made it through
    ligand_new_mol_final = handle_frag_check(ligand_new_mol_final)
    if ligand_new_mol_final is None:
        return None

    # Make sure there are no unassigned atoms which made it through. These are
    # very unlikely but possible
    ligand_new_mol_final = check_for_unassigned_atom(ligand_new_mol_final)
    if ligand_new_mol_final is None:
        return None

    ligand_new_mol = check_sanitization(ligand_new_mol_final)
    if ligand_new_mol is None:
        return None

    ligand_new_smiles = Chem.MolToSmiles(ligand_new_mol, isomericSmiles=True)

    return ligand_new_smiles


def run_main_smiles_merge(lig_string_1, lig_string_2):
    """
    This runs the main script for SmileMerge.

    Inputs:
    :param dict vars: User variables which will govern how the programs runs

    :param str lig_string_1: smile string for lig 1
    :param str lig_string_2: smile string for lig 2. example: lig_string_1 =
        "[N-] = [N+] = NCC(O)COc1cccc2ccccc12"; example: lig_string_2 = "C#
        CCOc1ccc2ccccc2c1CO"

    Returns:
    :returns: str ligand_new_smiles: smile string for the child ligand derived
        from lig_1 and lig_2. Returns None if it failed at any point.
    """

    # Sanitize
    lig_smile_1 = Chem.MolFromSmiles(lig_string_1, sanitize=False)
    lig_smile_2 = Chem.MolFromSmiles(lig_string_2, sanitize=False)

    # Sanitize, deprotanate, and reprotanate both molecules
    mol_1 = check_sanitization(lig_smile_1)
    mol_2 = check_sanitization(lig_smile_2)
    if mol_1 is None or mol_2 is None:
        return False

    mol_1 = handleHs(lig_smile_1)
    mol_2 = handleHs(lig_smile_2)

    if mol_1 is None or mol_2 is None:
        return None

    # make a list of the two rdkit.Chem.rdchem.Mol objects
    mols = [mol_1, mol_2]

    # Use the below mcs_H function for Most Common Substructure searching.
    # This will prevent broken rings.
    mcs_results = rdFMCS.FindMCS(
        mols,
        matchValences=False,
        ringMatchesRingOnly=True,
        completeRingsOnly=False,
        timeout=1,
    )

    if mcs_results.canceled is True:
        return None

    # confirm that this meets the minimum number of matching atoms
    if mcs_results.numAtoms < 4:
        return None

    ### Convert mcs_res from into usable and referable forms
    mcs_mol = Chem.MolFromSmarts(mcs_results.smartsString)

    # handle_mcs_align_labeling_and_cyclicbreaks
    mol_1, mol_2, mcs_mol = handle_mcs_align_labeling_and_cyclicbreaks(mol_1, mol_2, mcs_mol)

    if mol_1 is None or mol_2 is None or mcs_mol is None:
        return None

    rs_chosen_smiles = handle_dicts_and_select_b_groups(mol_1, mol_2, mcs_mol)
    if rs_chosen_smiles is None:
        return None

    ligand_new_mol = merge_smiles_with_core(rs_chosen_smiles, mcs_mol)
    if ligand_new_mol is None:
        return None

    ligand_new_smiles = process_ligand_new_mol(ligand_new_mol)
    return ligand_new_smiles
