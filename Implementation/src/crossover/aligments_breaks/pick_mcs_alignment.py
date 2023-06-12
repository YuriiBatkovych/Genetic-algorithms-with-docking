import random


def pick_mcs_alignment(mol_1, mol_2, common_core):
    """
    This will take the common substructure (aka the common_core) and find
    every match to it within each of the two ligands and produce to list of
    tuples for the atom index (Idx) relative to the indexing of the
    common_core substrure.

    It will then randomly pick combination of an alignment for mol_1 and
    mol_2.

    This should pick the alignments for future numbering.

    This should return picked_alignment (tuple)
        picked_alignment[0] is alignment for mol_1
        picked_alignment[1] is alignment for mol_2
        picked_alignment[2] is the atom idx of the common_core

        ie. picked_alignment[0][0] is the atom IDx for the atom in mol_1 which
        corresponds to the atom in the common substructure index as 0

        picked_alignment[0][1] is the atom IDx for the atom in mol_1 which
        corresponds to the atom in the common substructure index as 1

        picked_alignment[1][0] is the atom IDx for the atom in mol_2 which
        corresponds to the atom in the common substructure index as 0

        picked_alignment[1][1] is the atom IDx for the atom in mol_2 which
        corresponds to the atom in the common substructure index as 1

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_1: rdkit mol for ligand 1
    :param rdkit.Chem.rdchem.Mol mol_2: rdkit mol for ligand 2
    :param rdkit.Chem.rdchem.Mol common_core: rdkit mol for the shared core
        between mol_1 and mol_2

    Returns:
    :returns: tuple picked_alignment: tuple with 3 sub lists of the atoms IDx
        which correspond to their respective mol object in the order that they
        match the atom in the same index in the sublist for all 3 mols
        (mol_1,mol_2,common_core)
    """

    # Get the substructure match for the MCS within each ligand
    mol_1_match_idx = mol_1.GetSubstructMatches(
        common_core, uniquify=False, maxMatches=10
    )
    mol_2_match_idx = mol_2.GetSubstructMatches(
        common_core, uniquify=False, maxMatches=10
    )

    all_drug_pairings = []
    for mol_1_match in mol_1_match_idx:
        for mol_2_match in mol_2_match_idx:
            all_drug_pairings.append((mol_1_match, mol_2_match))

    if type(all_drug_pairings) != list:
        return None
    else:
        if len(all_drug_pairings) == 0 or type(all_drug_pairings[0]) != tuple:
            return None
        if len(all_drug_pairings[0]) == 0:
            return None

    alignment_choice = random.choice(all_drug_pairings)

    substruc_idx = tuple(list(range(0, len(common_core.GetAtoms()))))
    picked_alignment = (alignment_choice[0], alignment_choice[1], substruc_idx)
    return picked_alignment