def add_r_atom_isolabels(mol_1, mol_2):
    """
    Label the 1st atom in R-group with its idx + (1000 for lig1 and 2000 for
    lig 2) the 1000 vs 2000 label will be used to deconvelute later. These
    will be the iso values assigned to the 1st atom in an R-group

    Inputs:
    :param rdkit.Chem.rdchem.Mol mol_1: rdkit mol for ligand 1
    :param rdkit.Chem.rdchem.Mol mol_2: rdkit mol for ligand 2
    """

    # isotope label mol_1
    for atom in mol_1.GetAtoms():
        if atom.GetIsotope() < 9999:
            atom_iso = atom.GetIdx() + 1000
            atom.SetIsotope(atom_iso)

    # isotope label mol_2
    for atom in mol_2.GetAtoms():
        if atom.GetIsotope() < 9999:
            atom_iso = atom.GetIdx() + 2000
            atom.SetIsotope(atom_iso)