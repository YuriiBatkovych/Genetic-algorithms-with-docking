import random

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from src.mutations.validation import is_mol_ok, are_rings_ok, preprocess_mol

rdBase.DisableLog('rdApp.error')


def cut(mol):
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('[*]-;!@[*]')):
        return None

    bis = random.choice(mol.GetSubstructMatches(Chem.MolFromSmarts('[*]-;!@[*]')))  # single bond not in ring

    bs = [mol.GetBondBetweenAtoms(bis[0], bis[1]).GetIdx()]

    fragments_mol = Chem.FragmentOnBonds(mol, bs, addDummies=True, dummyLabels=[(1, 1)])

    try:
        return Chem.GetMolFrags(fragments_mol, asMols=True, sanitizeFrags=True)
    except ValueError:
        return None

    return None


def cut_ring(mol):
    for i in range(10):
        if random.random() < 0.5:
            if not mol.HasSubstructMatch(Chem.MolFromSmarts('[R]@[R]@[R]@[R]')):
                return None
            bis = random.choice(mol.GetSubstructMatches(Chem.MolFromSmarts('[R]@[R]@[R]@[R]')))
            bis = ((bis[0], bis[1]), (bis[2], bis[3]),)
        else:
            if not mol.HasSubstructMatch(Chem.MolFromSmarts('[R]@[R;!D2]@[R]')):
                return None
            bis = random.choice(mol.GetSubstructMatches(Chem.MolFromSmarts('[R]@[R;!D2]@[R]')))
            bis = ((bis[0], bis[1]), (bis[1], bis[2]),)

        bs = [mol.GetBondBetweenAtoms(x, y).GetIdx() for x, y in bis]

        fragments_mol = Chem.FragmentOnBonds(mol, bs, addDummies=True, dummyLabels=[(1, 1), (1, 1)])

        try:
            fragments = Chem.GetMolFrags(fragments_mol, asMols=True, sanitizeFrags=True)
            if len(fragments) == 2:
                return fragments
        except ValueError:
            return None

    return None


average_size = 39.15
size_stdev = 3.50


def crossover_ring(parent_A, parent_B):
    ring_smarts = Chem.MolFromSmarts('[R]')
    if not parent_A.HasSubstructMatch(ring_smarts) and not parent_B.HasSubstructMatch(ring_smarts):
        return None

    rxn_smarts1 = ['[*:1]~[1*].[1*]~[*:2]>>[*:1]-[*:2]', '[*:1]~[1*].[1*]~[*:2]>>[*:1]=[*:2]']
    rxn_smarts2 = ['([*:1]~[1*].[1*]~[*:2])>>[*:1]-[*:2]', '([*:1]~[1*].[1*]~[*:2])>>[*:1]=[*:2]']

    for i in range(10):
        fragments_A = cut_ring(parent_A)
        fragments_B = cut_ring(parent_B)

        if fragments_A is None or fragments_B is None:
            return None

        new_mol_trial = []
        for rs in rxn_smarts1:
            rxn1 = AllChem.ReactionFromSmarts(rs)
            new_mol_trial = []
            for fa in fragments_A:
                for fb in fragments_B:
                    new_mol_trial.append(rxn1.RunReactants((fa, fb))[0])

        new_mols = []
        for rs in rxn_smarts2:
            rxn2 = AllChem.ReactionFromSmarts(rs)
            for m in new_mol_trial:
                m = m[0]
                if is_mol_ok(m):
                    new_mols += list(rxn2.RunReactants((m,)))

        new_mols2 = []
        for m in new_mols:
            m = m[0]
            if is_mol_ok(m) and are_rings_ok(m):
                new_mols2.append(m)

        if len(new_mols2) > 0:
            return random.choice(new_mols2)

    return None


def crossover_non_ring(parent_A, parent_B):
    for i in range(10):
        fragments_A = cut(parent_A)
        fragments_B = cut(parent_B)
        if fragments_A is None or fragments_B is None:
            return None
        rxn = AllChem.ReactionFromSmarts('[*:1]-[1*].[1*]-[*:2]>>[*:1]-[*:2]')
        new_mol_trial = []
        for fa in fragments_A:
            for fb in fragments_B:
                new_mol_trial.append(rxn.RunReactants((fa, fb))[0])

        new_mols = []
        for mol in new_mol_trial:
            mol = mol[0]
            if is_mol_ok(mol):
                new_mols.append(mol)

        if len(new_mols) > 0:
            return random.choice(new_mols)

    return None


def crossover(mol1_smiles, mol2_smiles):
    parent_A = Chem.MolFromSmiles(mol1_smiles)
    parent_B = Chem.MolFromSmiles(mol2_smiles)
    parent_smiles = [mol1_smiles, mol2_smiles]
    try:
        Chem.Kekulize(parent_A, clearAromaticFlags=True)
        Chem.Kekulize(parent_B, clearAromaticFlags=True)

    except ValueError:
        pass

    number_of_attempts = 10
    for i in range(number_of_attempts):
        if random.random() <= 0.5:
            # print 'non-ring crossover'
            new_mol = crossover_non_ring(parent_A, parent_B)
            new_mol = preprocess_mol(new_mol)
            if new_mol is not None:
                new_smiles = Chem.MolToSmiles(new_mol)
                if new_smiles is not None and new_smiles not in parent_smiles:
                    return new_smiles
        else:
            # print 'ring crossover'
            new_mol = crossover_ring(parent_A, parent_B)
            new_mol = preprocess_mol(new_mol)
            if new_mol is not None:
                new_smiles = Chem.MolToSmiles(new_mol)
                if new_smiles is not None and new_smiles not in parent_smiles:
                    return new_smiles

    return None
