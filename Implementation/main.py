import logging as log
import os.path
import random
from collections import defaultdict, namedtuple

import torch
from rdkit import Chem
from torch.nn.functional import softmax, log_softmax
from tqdm import tqdm

from src.constants import ATOM_CATEGORIES
from src.conversions import convert_smiles_to_mols, convert_pdbqt_to_vectors, atoms_to_ints
from src.docking import dock_molecule
from src.files import load_n_smiles_from_file
from src.mutations.mutations import generate_list_of_mutations
from src.my_crossover.create_cross import crossover
from src.networks.ENN import ENN
from src.networks.LigandSMILES import LigandSMILES

log.basicConfig(level=log.INFO,
                format='%(asctime)s [%(threadName)-10.10s][%(levelname)-5.5s] %(message)s',
                datefmt='%H:%M:%S')

CrossoverGradData = namedtuple('CrossoverGradData', ['parent_1_idx', 'parent_2_idx',
                                                     'parent_1_docking_score', 'parent_2_docking_score',
                                                     'result_docking_score'])

MutationGradData = namedtuple('MutationGradData', ['original_idx', 'mutation_idx',
                                                   'original_docking_score', 'mutation_docking_score'])

receptor_info_list = [
    # ('5ht1a', '5ht1a_chembl_Ki_data.smi', '5ht1a_7E2X_cleaned.pdb', 101.897, 108.541, 111.118, 15, 15, 15),
    # ('5ht7r', '5ht7_chembl_Ki_data.smi', '5ht7r_7XTC_cleaned.pdb', 89.854, 102.530, 82.123, 15, 15, 15),
    ('beta2', 'beta2AR_chembl_Ki_data.smi', 'beta2_5D5A_cleaned.pdb', -35.132, 11.562, 3.035, 15, 15, 15),
    # ('d2',    'D2_chembl_Ki_data.smi', 'd2_6CM4_cleaned.pdb', 4.598, 9.999, -7.445, 15, 15, 15),
    # ('h1',    'H1_chembl_Ki_data.smi', 'h1_3RZE_cleaned.pdb', 9.705, 35.243, 21.005, 15, 15, 15)
]

# Variables to move to input parameters
ligands_dir_name = 'ligands'
proteins_dir_name = 'proteins'
outputs_dir_name = 'outputs'
docking_file_name = 'docking_result'
num_of_initial_smiles = 30
top_n_smiles = 10
exhaustiveness = 3
num_of_epochs = 5
num_of_crossovers = 0
num_of_mutations = 10
lr_crossover = 0.01
lr_mutations = 0.01


if __name__ == '__main__':

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    log.info(f'Using device {device}.')

    if not os.path.exists(outputs_dir_name):
        os.makedirs(outputs_dir_name)

    for receptor_info in receptor_info_list:
        protein_name, input_smiles_filename, target_protein_filename, \
        center_x, center_y, center_z, \
        size_x, size_y, size_z = receptor_info

        log.info(f'Current receptor: {protein_name}')

        input_smiles_path = os.path.join(ligands_dir_name, input_smiles_filename)
        initial_smiles = load_n_smiles_from_file(input_smiles_path, num_of_initial_smiles)
        log.info(f'Successfully loaded {len(initial_smiles)} smiles from file.')
        initial_smiles, initial_mols = convert_smiles_to_mols(initial_smiles)
        log.info(f'Successfully converted {len(initial_smiles)} smiles to mols.')

        log.info(f'Docking {len(initial_smiles)} smiles. This might take a while.')
        target_protein_path = os.path.join(proteins_dir_name, target_protein_filename)
        res_directory = os.path.join(outputs_dir_name, f'{outputs_dir_name}_{protein_name}_0')
        if not os.path.exists(res_directory):
            os.makedirs(res_directory)
        docking_result = [(dock_molecule(smiles,
                                         target_protein_path,
                                         os.path.join(res_directory, f'{docking_file_name}_{i}.pdbqt'),
                                         center_x, center_y, center_z,
                                         size_x, size_y, size_z,
                                         exhaustiveness),
                           f'{docking_file_name}_{i}.pdbqt')
                          for i, smiles in tqdm(enumerate(initial_smiles), total=len(initial_smiles))]

        docked_smiles_population = list(zip(initial_smiles,
                                            map(lambda x: x[0], docking_result),
                                            map(lambda x: x[1], docking_result)))
        docked_smiles_population = list(filter(lambda ds: ds[1] is not None, docked_smiles_population))
        log.info(f'Successfully docked {len(docked_smiles_population)} smiles.')

        docked_smiles_population = sorted(docked_smiles_population, key=lambda x: x[1])[:top_n_smiles]
        log.info(f'Taking {len(docked_smiles_population)} (no more than {top_n_smiles}) '
                 f'best smiles based on docking scores. '
                 f'Best score: {docked_smiles_population[0][1]}, worst score: {docked_smiles_population[-1][1]}')

        log.info(f'Saving top {top_n_smiles} smiles from epoch 0')
        with open(os.path.join(res_directory, 'best_smiles.smi'), 'w') as file:
            for docked_smiles in docked_smiles_population:
                file.write(f'{docked_smiles[0]} {docked_smiles[1]}\n')

        crossover_policy_net_1 = ENN(device)
        crossover_policy_net_2 = ENN(device)
        crossover_optimizer = torch.optim.Adam(list(crossover_policy_net_1.parameters())
                                               + list(crossover_policy_net_2.parameters()), lr=lr_crossover)

        mutation_policy_net_1 = ENN(device)
        mutation_policy_net_2 = LigandSMILES(device)
        mutation_optimizer = torch.optim.Adam(list(mutation_policy_net_1.parameters())
                                              + list(mutation_policy_net_2.parameters()), lr=lr_mutations)

        for e in range(num_of_epochs):
            current_epoch = e + 1
            log.info(f'Starting epoch number {current_epoch}/{num_of_epochs} for {protein_name}')

            res_directory = os.path.join(outputs_dir_name, f'{outputs_dir_name}_{protein_name}_{current_epoch}')
            if not os.path.exists(res_directory):
                os.makedirs(res_directory)

            ligands_atoms_information = []
            for docked_smiles in docked_smiles_population:
                docked_filename = os.path.join(outputs_dir_name,
                                               f'{outputs_dir_name}_{protein_name}_{current_epoch - 1}',
                                               docked_smiles[2])
                ligands_atoms_information.append(convert_pdbqt_to_vectors(docked_filename))

            # CROSSOVER
            crossover_grad_data = []

            probs_1 = []
            probs_2 = []
            for ligand_atom_information in ligands_atoms_information:
                encoded_atom_symbols = atoms_to_ints(list(map(lambda x: x[0], ligand_atom_information[0])),
                                                     ATOM_CATEGORIES)
                encoded_atom_symbols = torch.tensor(encoded_atom_symbols).to(device)
                spatial_positions = list(map(lambda x: x[1], ligand_atom_information[1]))
                spatial_positions = torch.tensor(spatial_positions).to(device)
                probs_1.append(crossover_policy_net_1(encoded_atom_symbols, spatial_positions))
                probs_2.append(crossover_policy_net_2(encoded_atom_symbols, spatial_positions))

            probs_1 = softmax(torch.stack(probs_1), dim=0)
            probs_2 = softmax(torch.stack(probs_2), dim=0)

            parents_1_indexes = random.choices(list(range(len(probs_1))), weights=probs_1, k=num_of_crossovers)
            parents_2_indexes = random.choices(list(range(len(probs_2))), weights=probs_2, k=num_of_crossovers)

            log.info(f'Performing {num_of_crossovers}*{num_of_crossovers} crossovers and dockings.')
            pbar = tqdm(total=num_of_crossovers*num_of_crossovers)

            docked_smiles_population_temp = []
            for parent_1_idx in parents_1_indexes:
                for parent_2_idx in parents_2_indexes:
                    if parent_1_idx == parent_2_idx:
                        pbar.update(1)
                        continue

                    parent_1 = docked_smiles_population[parent_1_idx]
                    parent_2 = docked_smiles_population[parent_2_idx]

                    parent_1_prob = probs_1[parent_1_idx]
                    parent_2_prob = probs_2[parent_2_idx]

                    parent_1_smiles = parent_1[0]
                    parent_2_smiles = parent_2[0]
                    result_smiles = crossover(parent_1_smiles, parent_2_smiles)

                    if result_smiles is None:
                        pbar.update(1)
                        continue

                    parent_1_docking_score = parent_1[1]
                    parent_2_docking_score = parent_2[1]

                    result_docking_score = dock_molecule(result_smiles,
                                                         target_protein_path,
                                                         os.path.join(res_directory,
                                                                      f'{docking_file_name}_{len(docked_smiles_population_temp)}.pdbqt'),
                                                         center_x, center_y, center_z,
                                                         size_x, size_y, size_z,
                                                         exhaustiveness)

                    if result_docking_score is None:
                        pbar.update(1)
                        continue

                    docked_smiles_population_temp.append((result_smiles,
                                                          result_docking_score,
                                                          f'{docking_file_name}_{len(docked_smiles_population_temp)}.pdbqt'))

                    crossover_grad_data.append(CrossoverGradData(parent_1_idx, parent_2_idx, parent_1_docking_score,
                                                                 parent_2_docking_score, result_docking_score))

                    pbar.update(1)

            pbar.close()
            log.info(f'Generated {len(crossover_grad_data)} crossovers. '
                     f'{num_of_crossovers*num_of_crossovers - len(crossover_grad_data)} crossover pairs failed.')

            # MUTATIONS

            mutation_grad_data = []

            probs_1 = []
            for ligand_atom_information in ligands_atoms_information:
                encoded_atom_symbols = atoms_to_ints(list(map(lambda x: x[0], ligand_atom_information[0])),
                                                     ATOM_CATEGORIES)
                encoded_atom_symbols = torch.tensor(encoded_atom_symbols).to(device)
                spatial_positions = list(map(lambda x: x[1], ligand_atom_information[1]))
                spatial_positions = torch.tensor(spatial_positions).to(device)
                probs_1.append(mutation_policy_net_1(encoded_atom_symbols, spatial_positions))

            probs_1 = softmax(torch.tensor(probs_1, device=device), dim=0)

            original_ligands_indexes = random.choices(list(range(len(probs_1))), weights=probs_1, k=num_of_mutations)

            log.info(f'Performing {num_of_mutations}*2 mutations and dockings.')
            pbar = tqdm(total=num_of_mutations*2)

            mutations_dict = {}
            for original_ligand_index in original_ligands_indexes:
                original_ligand = docked_smiles_population[original_ligand_index]

                mutations = generate_list_of_mutations(Chem.MolFromSmiles(original_ligand[0]), num_of_mutations)
                mutations = [Chem.MolToSmiles(mol) for mol in mutations]

                _, probs_2 = mutation_policy_net_2(original_ligand[0], mutations)

                mutations_indexes = random.choices(list(range(len(mutations))), weights=probs_2, k=2)

                for i in mutations_indexes:
                    result_docking_score = dock_molecule(mutations[i],
                                                         target_protein_path,
                                                         os.path.join(res_directory,
                                                                      f'{docking_file_name}_{len(docked_smiles_population_temp)}.pdbqt'),
                                                         center_x, center_y, center_z,
                                                         size_x, size_y, size_z,
                                                         exhaustiveness)

                    if result_docking_score is None:
                        pbar.update(1)
                        continue

                    docked_smiles_population_temp.append((mutations[i],
                                                          result_docking_score,
                                                          f'{docking_file_name}_{len(docked_smiles_population_temp)}.pdbqt'))

                    mutation_grad_data.append(MutationGradData(original_ligand_index, i, original_ligand[1], result_docking_score))
                    pbar.update(1)

                mutations_dict[original_ligand_index] = mutations

            pbar.close()
            log.info(f'Generated {len(mutation_grad_data)} mutations. '
                     f'{num_of_mutations*2 - len(mutation_grad_data)} mutations failed.')

            # CROSSOVER OPTIMIZATION

            if len(crossover_grad_data) > 0:
                log.info('Optimizing crossover network.')

            for grad_data in crossover_grad_data:
                probs_1 = []
                probs_2 = []
                for ligand_atom_information in ligands_atoms_information:
                    encoded_atom_symbols = atoms_to_ints(list(map(lambda x: x[0], ligand_atom_information[0])),
                                                         ATOM_CATEGORIES)
                    encoded_atom_symbols = torch.tensor(encoded_atom_symbols).to(device)
                    spatial_positions = list(map(lambda x: x[1], ligand_atom_information[1]))
                    spatial_positions = torch.tensor(spatial_positions).to(device)
                    probs_1.append(crossover_policy_net_1(encoded_atom_symbols, spatial_positions))
                    probs_2.append(crossover_policy_net_2(encoded_atom_symbols, spatial_positions))

                probs_1 = log_softmax(torch.stack(probs_1), dim=0)
                probs_2 = log_softmax(torch.stack(probs_2), dim=0)

                parent_1_log_prob = probs_1[grad_data.parent_1_idx]
                parent_2_log_prob = probs_2[grad_data.parent_2_idx]

                _, _, parent_1_docking_score, parent_2_docking_score, result_docking_score = grad_data

                reward = -result_docking_score - max(-parent_1_docking_score, -parent_2_docking_score)
                log_prob = (parent_1_log_prob + parent_2_log_prob) * reward

                crossover_optimizer.zero_grad()
                log_prob.backward()
                crossover_optimizer.step()

            # MUTATION OPTIMIZATION

            if len(mutation_grad_data) > 0:
                log.info('Optimizing mutation network.')

            for grad_data in mutation_grad_data:
                original_index, mutation_index, original_docking_score, mutation_docking_score = grad_data

                probs_1 = []
                for ligand_atom_information in ligands_atoms_information:
                    encoded_atom_symbols = atoms_to_ints(list(map(lambda x: x[0], ligand_atom_information[0])),
                                                         ATOM_CATEGORIES)
                    encoded_atom_symbols = torch.tensor(encoded_atom_symbols).to(device)
                    spatial_positions = list(map(lambda x: x[1], ligand_atom_information[1]))
                    spatial_positions = torch.tensor(spatial_positions).to(device)
                    probs_1.append(mutation_policy_net_1(encoded_atom_symbols, spatial_positions))

                probs_1 = log_softmax(torch.stack(probs_1), dim=0)
                probs_2, _ = mutation_policy_net_2(docked_smiles_population[original_index][0],
                                                   mutations_dict[original_index])

                parent_1_log_prob = probs_1[original_index]
                parent_2_log_prob = probs_2[mutation_index]

                reward = original_docking_score - mutation_docking_score
                log_prob = (parent_1_log_prob + parent_2_log_prob) * reward

                mutation_optimizer.zero_grad()
                log_prob.backward()
                mutation_optimizer.step()

            docked_smiles_population = docked_smiles_population_temp
            docked_smiles_population = sorted(docked_smiles_population, key=lambda x: x[1])[:top_n_smiles]
            log.info(f'Taking {len(docked_smiles_population)} (no more than {top_n_smiles}) '
                     f'best smiles based on docking scores. '
                     f'Best score: {docked_smiles_population[0][1]}, worst score: {docked_smiles_population[-1][1]}')

            log.info(f'Saving top {top_n_smiles} smiles from epoch {num_of_epochs}')
            with open(os.path.join(res_directory, 'best_smiles.smi'), 'w') as file:
                for docked_smiles in docked_smiles_population:
                    file.write(f'{docked_smiles[0]} {docked_smiles[1]}\n')
