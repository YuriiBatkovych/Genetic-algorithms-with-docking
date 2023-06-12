import rdkit


def remove_iso_labels(mol, list_of_idx_to_remove):
    for i in list_of_idx_to_remove:
        atom = mol.GetAtomWithIdx(i)
        atom.SetIsotope(0)


def remove_all_isolabels(rw_core_merg):
    """
    Remove all the isotope labels from a molecule. One of the finalizing steps
    used when LigSmiles is nearly complete.

    Inputs:
    :param rdkit.Chem.rdchem.RWMol rw_core_merg: a read-write rdkit molecule
        of the child molecule (after R-groups have been added) to have the
        isotopes to be removed

    Returns:
    :returns: rdkit.Chem.rdchem.RWMol rw_core_merg: the read-write rdkit
        molecule of the child molecule with all the isotope labels removed
    """

    if rw_core_merg is None:
        return None

    # If mol is wrong data type (excluding None) raise TypeError
    if (
            type(rw_core_merg) != rdkit.Chem.rdchem.Mol
            and type(rw_core_merg) != rdkit.Chem.rdchem.RWMol
    ):
        printout = "rw_core_merg is the wrong data type. \n"
        printout = (
            printout
            + "Input should be a rdkit.Chem.rdchem.Mol or rdkit.Chem.rdchem.RWMol\n"
        )
        printout = printout + "Input mol was {} type.".format(type(rw_core_merg))
        raise TypeError(printout)

    for atom in rw_core_merg.GetAtoms():
        atom.SetIsotope(0)

    return rw_core_merg