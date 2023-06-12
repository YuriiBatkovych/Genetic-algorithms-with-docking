from torch import nn
import torch
from torch.nn.functional import one_hot, softmax

from src.constants import ATOM_CATEGORIES


class ENN(nn.Module):
    def __init__(self, device, latent_size=50, vocab_size=len(ATOM_CATEGORIES)):
        super(ENN, self).__init__()
        self.device = device
        self.latent_size = latent_size
        self.vocab_size = vocab_size
        self.learnable_matrix = nn.Parameter(torch.randn(vocab_size, latent_size)).to(device)

        self.mlp_e = self.create_mlp(2 * self.latent_size + 1, self.latent_size)
        self.mlp_x = self.create_mlp(self.latent_size, 1)
        self.mlp_h = self.create_mlp(self.latent_size + 1, self.latent_size)
        self.output_mlp = nn.Linear(self.latent_size, 1).to(device)

    def create_mlp(self, in_size, out_size):
        return nn.Sequential(
            nn.Linear(in_size, self.latent_size),
            nn.Tanh(),
            nn.Linear(self.latent_size, out_size),
            nn.Tanh()).to(self.device)

    def forward(self, atom_categories, atom_coordinates):
        if len(atom_categories) != len(atom_coordinates):
            raise ValueError('Atom categories and atom coordinates lists must be of same size')

        # initialization
        Z = atom_coordinates
        H = self.init_H(atom_categories)

        # single iteration
        W = self.calculate_W(H, Z)
        V = self.calculate_V(W)
        H = self.calculate_H(H, V)

        H = torch.sum(H, dim=0)
        return self.output_mlp(H)[0]

    def init_H(self, atom_categories):
        category_lines = []
        for category in atom_categories:
            category_one_hot = one_hot(category, len(ATOM_CATEGORIES)).float()
            category_line = torch.matmul(self.learnable_matrix.T, category_one_hot)
            category_lines.append(category_line)
        return torch.stack(category_lines)

    def calculate_W(self, H, Z):
        W = torch.zeros((H.shape[0], H.shape[0], self.latent_size), device=self.device)
        for i in range(H.shape[0]):
            for j in range(H.shape[0]):
                z_diff_norm_squared = torch.pow(torch.norm(Z[i] - Z[j]), 2).reshape(-1)
                hi_hj_z_diff_concat = torch.cat((H[i], H[j], z_diff_norm_squared))
                W[i][j] = self.mlp_e(hi_hj_z_diff_concat)
        return W

    def calculate_V(self, W):
        V = torch.zeros((W.shape[0]), device=self.device)
        for i in range(W.shape[0]):
            V[i] = torch.sum(torch.stack([W[i][j] for j in range(W.shape[0]) if i != j]))
        return V

    def calculate_H(self, H, V):
        new_H = torch.zeros(H.shape, device=self.device)
        for i in range(H.shape[0]):
            h_v_concat = torch.cat((H[i], V[i].reshape(1,)))
            new_H[i] = self.mlp_h(h_v_concat)
        return new_H


if __name__ == '__main__':
    enn = ENN(torch.device('cpu'))
    atom_categories = torch.tensor([4, 1, 0])
    atom_coordinates = torch.tensor([[0.2, 0.3, 0.1], [0.5, 0.31, 0.3], [0.12, 0.32, 0.1]])

    print('---- LEARNABLE MATRIX TEST ----')
    a = one_hot(atom_categories, len(ATOM_CATEGORIES)).float()
    print(a)
    H = [torch.matmul(enn.learnable_matrix.T, one_hot(a, len(ATOM_CATEGORIES)).float()) for a in atom_categories]
    print(H)

    print('---- ENN TEST ----')
    distribution = softmax(enn(atom_categories, atom_coordinates), dim=0)
    print(distribution)
