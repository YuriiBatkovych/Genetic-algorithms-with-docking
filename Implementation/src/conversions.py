import logging as log
from typing import List

import numpy as np
import torch
from Bio.PDB import PDBParser
import warnings
from Bio import BiopythonWarning
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

warnings.simplefilter('ignore', BiopythonWarning)

def smiles_to_fingerprints(smiles_string: str, n_bits: int = 2048, device: str = 'cpu') -> torch.Tensor:
    mol = Chem.MolFromSmiles(smiles_string)
    Chem.SanitizeMol(mol)
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
    features = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fingerprint, features)
    fingerprint = torch.from_numpy(features).to(device).float()
    return fingerprint


def convert_pdbqt_to_vectors(pdbqt_file):
    parser = PDBParser()
    structure = parser.get_structure('structure', pdbqt_file)

    symbol_vector = []
    position_vector = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    symbol_vector.append(atom.get_name())
                    position_vector.append(atom.get_coord())

    return symbol_vector, position_vector


class AtomSymbolConverter:
    def __init__(self, symbols_list: List):
        self.symbols_list = symbols_list

    def atom_to_int(self, atom: str):
        if atom in self.symbols_list:
            return self.symbols_list.index(atom)
        return len(self.symbols_list) - 1


def atoms_to_ints(atoms: List, symbols_list: List) -> List:
    converter = AtomSymbolConverter(symbols_list)
    return [converter.atom_to_int(atom) for atom in atoms]


def convert_smiles_to_mols(smiles_list: List):
    updated_smiles_list = []
    mols_list = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)

        if mol is not None:
            updated_smiles_list.append(smiles)
            mols_list.append(mol)

    num_failures = len(smiles_list) - len(mols_list)
    if num_failures > 0:
        log.warning(f'{num_failures} SMILES failed to convert to mols.')

    return updated_smiles_list, mols_list
