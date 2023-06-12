import logging as log
import random

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from src.mutations.validation import preprocess_mol, are_rings_ok, is_mol_ok

rdBase.DisableLog('rdApp.error')


def delete_atom():
    choices = ['[*:1]~[D1:2]>>[*:1]', '[*:1]~[D2:2]~[*:3]>>[*:1]-[*:3]',
               '[*:1]~[D3:2](~[*;!H0:3])~[*:4]>>[*:1]-[*:3]-[*:4]',
               '[*:1]~[D4:2](~[*;!H0:3])(~[*;!H0:4])~[*:5]>>[*:1]-[*:3]-[*:4]-[*:5]',
               '[*:1]~[D4:2](~[*;!H0;!H1:3])(~[*:4])~[*:5]>>[*:1]-[*:3](-[*:4])-[*:5]']
    p = [0.25, 0.25, 0.25, 0.1875, 0.0625]

    return np.random.choice(choices, p=p)


def append_atom():
    choices = [['single', ['C', 'N', 'O', 'F', 'S', 'Cl', 'Br'], 7 * [1.0 / 7.0]],
               ['double', ['C', 'N', 'O'], 3 * [1.0 / 3.0]],
               ['triple', ['C', 'N'], 2 * [1.0 / 2.0]]]
    p_BO = [0.60, 0.35, 0.05]

    index = np.random.choice(list(range(3)), p=p_BO)

    BO, atom_list, p = choices[index]
    new_atom = np.random.choice(atom_list, p=p)

    if BO == 'single':
        return '[*;!H0:1]>>[*:1]X'.replace('X', '-' + new_atom)
    if BO == 'double':
        return '[*;!H0;!H1:1]>>[*:1]X'.replace('X', '=' + new_atom)

    return '[*;H3:1]>>[*:1]X'.replace('X', '#' + new_atom)


def insert_atom():
    choices = [['single', ['C', 'N', 'O', 'S'], 4 * [1.0 / 4.0]],
               ['double', ['C', 'N'], 2 * [1.0 / 2.0]],
               ['triple', ['C'], [1.0]]]
    p_BO = [0.60, 0.35, 0.05]

    index = np.random.choice(list(range(3)), p=p_BO)

    BO, atom_list, p = choices[index]
    new_atom = np.random.choice(atom_list, p=p)

    if BO == 'single':
        return '[*:1]~[*:2]>>[*:1]X[*:2]'.replace('X', new_atom)
    if BO == 'double':
        return '[*;!H0:1]~[*:2]>>[*:1]=X-[*:2]'.replace('X', new_atom)

    return '[*;!R;!H1;!H0:1]~[*:2]>>[*:1]#X-[*:2]'.replace('X', new_atom)


def change_bond_order():
    choices = ['[*:1]!-[*:2]>>[*:1]-[*:2]', '[*;!H0:1]-[*;!H0:2]>>[*:1]=[*:2]',
               '[*:1]#[*:2]>>[*:1]=[*:2]', '[*;!R;!H1;!H0:1]~[*:2]>>[*:1]#[*:2]']
    p = [0.45, 0.45, 0.05, 0.05]

    return np.random.choice(choices, p=p)


def delete_cyclic_bond():
    return '[*:1]@[*:2]>>([*:1].[*:2])'


def add_ring():
    choices = ['[*;!r;!H0:1]~[*;!r:2]~[*;!r;!H0:3]>>[*:1]1~[*:2]~[*:3]1',
               '[*;!r;!H0:1]~[*!r:2]~[*!r:3]~[*;!r;!H0:4]>>[*:1]1~[*:2]~[*:3]~[*:4]1',
               '[*;!r;!H0:1]~[*!r:2]~[*:3]~[*:4]~[*;!r;!H0:5]>>[*:1]1~[*:2]~[*:3]~[*:4]~[*:5]1',
               '[*;!r;!H0:1]~[*!r:2]~[*:3]~[*:4]~[*!r:5]~[*;!r;!H0:6]>>[*:1]1~[*:2]~[*:3]~[*:4]~[*:5]~[*:6]1']
    p = [0.05, 0.05, 0.45, 0.45]

    return np.random.choice(choices, p=p)


def change_atom(mol):
    choices = ['#6', '#7', '#8', '#9', '#16', '#17', '#35']
    p = [0.15, 0.15, 0.14, 0.14, 0.14, 0.14, 0.14]

    X = np.random.choice(choices, p=p)
    while not mol.HasSubstructMatch(Chem.MolFromSmarts('[' + X + ']')):
        X = np.random.choice(choices, p=p)
    Y = np.random.choice(choices, p=p)
    while Y == X:
        Y = np.random.choice(choices, p=p)

    return '[X:1]>>[Y:1]'.replace('X', X).replace('Y', Y)


def generate_single_mutation(mol: Chem.rdchem.Mol):
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except ValueError:
        return None

    for i in range(10):
        rxn_smarts_list = [
            insert_atom(),
            change_bond_order(),
            delete_cyclic_bond(),
            add_ring(),
            delete_atom(),
            change_atom(mol),
            append_atom(),
        ]
        rxn_smarts = np.random.choice(rxn_smarts_list)

        rxn = AllChem.ReactionFromSmarts(rxn_smarts)

        new_mol_trial = rxn.RunReactants((mol,))

        new_mols = []
        for m in new_mol_trial:
            m = m[0]

            m = preprocess_mol(m)
            if m is None:
                continue

            mol_ok = is_mol_ok(m)
            rings_ok = are_rings_ok(m)

            if mol_ok and rings_ok:
                new_mols.append(m)

        if len(new_mols) > 0:
            return random.choice(new_mols)

    log.warning(f'Failed to mutate molecule: {mol}')
    return None


def generate_list_of_mutations(mol: Chem.rdchem.Mol, list_size: int):
    mutations_list = []
    counter = 0

    while len(mutations_list) < list_size:

        if counter > 5 * list_size:
            log.warning(f'Failed to generate list of mutations for molecule {mol}')
            return None

        mutation = generate_single_mutation(mol)
        if mutation is not None:
            mutations_list.append(mutation)

        counter += 1

    return mutations_list


if __name__ == '__main__':
    mol = Chem.MolFromSmiles('CC1(C)CC(=O)N(CCCCN2Cc3ccccc3C2)C(=O)C1')
    mol = generate_single_mutation(mol)
    smiles = Chem.MolToSmiles(mol)
    print(smiles)

    # original CC1(C)CC(=O)N(CCCCN2Cc3ccccc3C2)C(=O)C1
    # insert atom CNC1(C)CC(=O)N(CCCCN2Cc3ccccc3C2)C(=O)C1
    # change bond order CC1(C)CC(=O)N(CC=CCN2Cc3ccccc3C2)C(=O)C1
    # delete cyclic bond C=CC=C1CN(CCCCN2C(=O)CC(C)(C)CC2=O)CC1=C
    # add ring CC1(C)CC(=O)N(CC2CC2N2Cc3ccccc3C2)C(=O)C1
    # delete atom CC1(C)CC(=O)N(CCCCN2Cc3ccccc3C2)C1O
    # change atom CC1(C)CC(=O)N(CCCC[SH]2Cc3ccccc3C2)C(=O)C1
    # append atom CC1C(=O)N(CCCCN2Cc3ccccc3C2)C(=O)CC1(C)C
