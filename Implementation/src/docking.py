import logging as log
import os
import subprocess

from rdkit import Chem
from rdkit.Chem import AllChem


def optimize_conformation(mol):
    try:
        mol = Chem.AddHs(mol)  # Adds hydrogens to make optimization more accurate
        AllChem.EmbedMolecule(mol)  # Adds 3D positions
        AllChem.MMFFOptimizeMolecule(mol)  # Improves the 3D positions using a force-field method
        return mol
    except ValueError as e:
        log.warning(e)
        return None


def dock_molecule(target_smile: str, protein_pdb: str, res_file: str, center_x: float, center_y: float, center_z: float,
                  size_x: float = 30, size_y: float = 30, size_z: float = 30, exhaustiveness: int = 8) -> float or None:
    mol = Chem.MolFromSmiles(target_smile)
    mol = optimize_conformation(mol)
    if mol is None:
        return None
    Chem.MolToMolFile(mol, 'molecule.mol')
    babel = 'obabel -imol molecule.mol -omol2 -O molecule.mol2'
    subprocess.run(babel, shell=True, capture_output=True, text=True)
    os.remove('molecule.mol')

    smina = f'smina -r {protein_pdb} -l molecule.mol2 --center_x {center_x} --center_y {center_y} --center_z {center_z}' \
            f' --size_x {size_x} --size_y {size_y} --size_z {size_z} --exhaustiveness {exhaustiveness} --out {res_file}'
    result = subprocess.run(smina, shell=True, capture_output=True, text=True)

    return float(result.stdout.split('\n')[29].split()[1])
