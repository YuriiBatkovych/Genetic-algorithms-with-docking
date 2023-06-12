import numpy as np
import torch
import torch.nn.functional as F
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from torch import nn


def smiles_to_fingerprint(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    Chem.SanitizeMol(mol)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    features = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, features)
    fingerprint = torch.from_numpy(features).float().view(1, -1)
    return fingerprint


class LigandSMILES(nn.Module):
    def __init__(self, device):
        super(LigandSMILES, self).__init__()
        self.device = device
        self.ligand_mlp = nn.Linear(2048, 100).to(device)
        self.product_mlp = nn.Linear(2048, 100).to(device)
        self.output_mlp = nn.Linear(200, 1).to(device)

    def forward(self, ligand_smiles, product_smiles_list):
        n = len(product_smiles_list)
        ligand_fp = smiles_to_fingerprint(ligand_smiles).to(self.device)
        ligand_embedding = F.relu(self.ligand_mlp(ligand_fp))
        ligand_embedding = ligand_embedding.repeat(n, 1)

        product_fps = [smiles_to_fingerprint(smiles).to(self.device) for smiles in product_smiles_list]
        product_fps = torch.cat(product_fps, 0)
        product_embeddings = F.relu(self.product_mlp(product_fps))

        latent_variable = torch.cat([ligand_embedding, product_embeddings], 1)
        output = self.output_mlp(latent_variable).view(-1)
        log_output = F.log_softmax(output, dim=0)
        prob_output = F.softmax(output, dim=0)
        return log_output, prob_output
