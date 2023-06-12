import numpy as np
from rdkit import Chem

import src.gypsum_dl.gypsum_dl.MolObjectHandling as MOH


def preprocess_mol(mol: Chem.rdchem.Mol):
    mol = MOH.check_sanitization(mol)
    if mol is None:
        return None

    mol = MOH.handle_frag_check(mol)
    if mol is None:
        return None

    mol = MOH.check_for_unassigned_atom(mol)
    if mol is None:
        return None

    mol = MOH.try_reprotanation(mol)
    if mol is None:
        return None

    mol = MOH.try_deprotanation(mol)
    if mol is None:
        return None

    mol = MOH.check_sanitization(mol)
    if mol is None:
        return None

    return mol


def are_rings_ok(mol: Chem.rdchem.Mol) -> bool:
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('[R]')):
        return True

    ring_allene = mol.HasSubstructMatch(Chem.MolFromSmarts('[R]=[R]=[R]'))

    cycle_list = mol.GetRingInfo().AtomRings()
    max_cycle_length = max([len(j) for j in cycle_list])
    macro_cycle = max_cycle_length > 6

    double_bond_in_small_ring = mol.HasSubstructMatch(Chem.MolFromSmarts('[r3,r4]=[r3,r4]'))

    return not ring_allene and not macro_cycle and not double_bond_in_small_ring


average_size = 39.15
size_stdev = 3.50


def is_mol_ok(mol: Chem.rdchem.Mol) -> bool:
    try:
        Chem.SanitizeMol(mol)
        target_size = size_stdev * np.random.randn() + average_size  # parameters set in GA_mol
        if 5 < mol.GetNumAtoms() < target_size:
            return True
        else:
            return False
    except ValueError:
        return False
