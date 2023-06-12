from rdkit import Chem


def r_group_list(mol, core_mol):
    replace_core_mol = Chem.ReplaceCore(
        mol, core_mol, labelByIndex=True, replaceDummies=True, requireDummyMatch=False
    )

    if len(replace_core_mol.GetAtoms()) == 0:
        return None

    return replace_core_mol
