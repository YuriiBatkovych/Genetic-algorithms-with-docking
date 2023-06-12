# Genetic algorithms with docking
Project carried out as part of the Machine Learning in Drug Design course at the Jagiellonian University in Krakow under the mentorship of dr Sabina Podlewska and mgr Tomasz Danel. The group implementing the project: Kamil Harłacz, Artur Górak, Piotr Faron, Yurii Titov.

## Tasks
* Implement a genetic algorithm generating chemical compounds with improved affinity (better results of
docking).
* Compare different methods of compound mutation.

## Research hypotheses
- The applied genetic algorithm improves the function
docking assessments from population to population.
- The docking algorithm is optimized to
ensure a satisfactory docking speed at
high accuracy.
- Mutations give results comparable to crossovers.

## Genetic algorithm
The implemented genetic algorithm is based on the article Reinforced Genetic Algorithm for
Structure-based Drug Design (https://arxiv.org/pdf/2211.16508.pdf). Our modifications of the mentioned algorithm are the following:
- different mutations
- different crossovers
- use of a different docking program (Smina)

## How to run our code?
To run the algorithm you could use the `mldd23` conda environment.
The list of dependencies for `mldd23` can be found in `environment-gpu.yml` file.
After activating the env, run `main.py` file.

The global parameters can be found at the top of `main.py` file. Currently configurable global params are:

```
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
```

You should be able to test the algorithm for different ligands and target receptors.
To do that you need to:
- add your targets to the directory specified in `proteins_dir_name` directory (default `proteins`) in .pdb file format 
- add the initial populations of ligands active to the targets in .smi files format to the directory specified in `ligands_dir_name` (default `ligands`) 
- modify receptor_info_list in `main.py` file. The fields are, from the left: 
  - name of the receptor,
  - name of the .smi file with the initial population,
  - name of the .pdb file describing the biological target,
  - 3D coordinates of the pocket center,
  - 3D size of the pocket
  
  An example of the receptor_info_list tuple element:
  `('5ht1a', '5ht1a_chembl_Ki_data.smi', '5ht1a_7E2X_cleaned.pdb', 101.897, 108.541, 111.118, 15, 15, 15)`

Best SMILES generated for each epoch can be found in `outputs/outputs_receptor_epoch/best_smiles.smi`.
Each line contains a generated SMILES sequence and its docking score.

## Dependencies
In our project we are basing our code on Reinforced Genetic Algorithm for
Structure-based Drug Design article and use open-source library gypsum_dl for checking validity of generated SMILES.

## Results
The description and visualisations of the final results can be found in the `final_presentation.pdf` file.