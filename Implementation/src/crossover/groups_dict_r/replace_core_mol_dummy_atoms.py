def replace_core_mol_dummy_atoms(replace_core_mol):
    """
    This function will replace the dummy atoms (*) with the isotope label from
    the core atoms. example:
        mol = Chem.MolFromSmiles("[10000N-]=[10001N+]=[10002N][10003CH2][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")
        mcs = Chem.MolFromSmiles("[10003CH3][10002N]=[10001N+]=[10000NH]")
        replace_core = Chem.MolFromSmiles("[3*][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]")

        resulting replace_core = '[10003*][2004CH]1[2005NH2+][2006CH2][2007CH]([2008OH])[2009CH]([2010OH])[2011CH]1[2012OH]'

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol: an rdkit molecule
    :param rdkit.Chem.rdchem.Mol mcs: an rdkit molecule for the shared common
        core
    :param rdkit.Chem.rdchem.Mol replace_core_mol: the mol with the MCS
        anchors labeled with * and an isotope label of the idx of the core anchor
        atom

    Returns:
    :returns: rdkit.Chem.rdchem.Mol replace_core_mol: an rdkit molecule with
        the common core removed from a ligand fragments the mol which can be used
        to make lists of R-groups. The * atoms will be isotope labeled with the
        isotope label from the core.
    """

    anchor_dict = {}
    anchor_to_set_dict = {}
    for atom in replace_core_mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            anchor_iso = atom.GetIsotope() + 10000
            neighbors = atom.GetNeighbors()
            tmp = []
            for n_atom in neighbors:
                tmp.append(n_atom.GetIsotope())
            anchor_dict[anchor_iso] = tmp

            anchor_to_set_dict[atom.GetIdx()] = anchor_iso

    for idx in list(anchor_to_set_dict.keys()):
        atom = replace_core_mol.GetAtomWithIdx(idx)
        anchor_iso = anchor_to_set_dict[idx]
        atom.SetIsotope(anchor_iso)

    return replace_core_mol