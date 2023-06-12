from rdkit import Chem


def convert_sdf_file_to_smiles_file(input_sdf_filename, output_smiles_filename):
    suppl = Chem.SDMolSupplier(input_sdf_filename)
    with open(output_smiles_filename, 'w') as f:
        for mol in suppl:
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                f.write(smiles + '\n')


def load_n_smiles_from_file(smiles_filename, num_smiles):
    smiles_list = []

    with open(smiles_filename, 'r') as file:
        for i, line in enumerate(file):
            if i >= num_smiles:
                break

            smiles = line.strip().split()[0]
            smiles_list.append(smiles)

    return smiles_list
