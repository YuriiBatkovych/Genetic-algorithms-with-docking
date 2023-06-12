from src.gypsum_dl.gypsum_dl.MolObjectHandling import check_sanitization, try_deprotanation


def handleHs(mol):
    """
    Given a rdkit.Chem.rdchem.Mol this script will sanitize the molecule, remove all non-explicit H's
    and add back on all implicit H's. This is to control for any discrepencies in the smiles strings or presence and
    absense of H's.
    If it fails it will return a None rather than causing the outer script to fail. Handled here so there are no problems later.

    Inputs:
    :param rdkit.Chem.rdchem.Mol sanitized_deprotanated_mol: an rdkit molecule already sanitized and deprotanated.

    Returns:
    :returns: rdkit.Chem.rdchem.Mol mol: an rdkit molecule with H's handled (either added or removed) and sanitized.
                                            it returns None if H's can't be added or if sanitation fails
    """
    mol = check_sanitization(mol)
    if mol is None:
        return None

    mol = try_deprotanation(mol)
    if mol is None:
        return None

    return mol
