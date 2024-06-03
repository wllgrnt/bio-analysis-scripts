"""
Reads the CIF files output by alphafold and gets the average confidence score for each residue.
"""
import gemmi
import numpy as np
import matplotlib.pyplot as plt

def extract_cif_info_by_chain(cif_file):
    # Read the CIF file
    structure = gemmi.read_structure(cif_file)
    
    # Dictionary to store extracted information by chain
    chain_info = {}

    # Iterate over all models, chains, residues, and atoms in the structure
    for model in structure:
        for chain in model:
            chain_name = chain.name
            if chain_name not in chain_info:
                chain_info[chain_name] = {"residue_info": {}, "alpha_carbons": []}
            for residue in chain:
                residue_number = residue.seqid.num
                b_factors = []
                for atom in residue:
                    b_factor = atom.b_iso
                    b_factors.append(b_factor)
                    if atom.name == 'CA':  # Check for alpha carbon
                        chain_info[chain_name]["alpha_carbons"].append((residue_number, b_factor))
                if b_factors:
                    avg_b_factor = np.mean(b_factors)
                    std_b_factor = np.std(b_factors)
                    chain_info[chain_name]["residue_info"][residue_number] = (avg_b_factor, std_b_factor)

    return chain_info

def plot_b_factors_by_chain(chain_info):
    for chain_name, data in chain_info.items():
        residues = sorted(data["residue_info"].keys())
        avg_b_factors = [data["residue_info"][res][0] for res in residues]
        std_b_factors = [data["residue_info"][res][1] for res in residues]

        alpha_carbon_residues = [ac[0] for ac in data["alpha_carbons"]]
        alpha_carbon_b_factors = [ac[1] for ac in data["alpha_carbons"]]

        fig, axes = plt.subplots(3, 1, figsize=(10, 15))

        # Plot average B factors
        axes[0].plot(residues, avg_b_factors, marker='o')
        axes[0].set_title(f'Average B Factor per Residue (Chain {chain_name})')
        axes[0].set_xlabel('Residue Number')
        axes[0].set_ylabel('Average B Factor')

        # Plot standard deviation of B factors
        axes[1].plot(residues, std_b_factors, marker='o')
        axes[1].set_title(f'Standard Deviation of B Factor per Residue (Chain {chain_name})')
        axes[1].set_xlabel('Residue Number')
        axes[1].set_ylabel('Standard Deviation B Factor')

        # Plot alpha carbon B factors
        axes[2].plot(alpha_carbon_residues, alpha_carbon_b_factors, marker='o', linestyle='None')
        axes[2].set_title(f'Alpha Carbon B Factor per Residue (Chain {chain_name})')
        axes[2].set_xlabel('Residue Number')
        axes[2].set_ylabel('Alpha Carbon B Factor')

        plt.tight_layout()
        plt.show()

# Example usage
cif_file = "input_folder/fold_trak1_dimer_model_0.cif"
chain_info = extract_cif_info_by_chain(cif_file)
plot_b_factors_by_chain(chain_info)

# def analyze_cif_structure(cif_file):
#     # Read the CIF file
#     structure = gemmi.read_structure(cif_file)
    
#     # Sets to store unique models, chains, residues, and alpha carbons
#     models_set = set()
#     chains_set = set()
#     residues_set = set()
#     alpha_carbons_set = set()

#     # Iterate over all models, chains, residues, and atoms in the structure
#     for model in structure:
#         models_set.add(model.name)
#         for chain in model:
#             chains_set.add(chain.name)
#             for residue in chain:
#                 residue_number = residue.seqid.num
#                 residues_set.add((chain.name, residue_number))
#                 for atom in residue:
#                     if atom.name == 'CA':  # Check for alpha carbon
#                         alpha_carbons_set.add((chain.name, residue_number, atom.name, atom.b_iso))

#     print(f"Number of unique models: {len(models_set)}")
#     print(f"Number of unique chains: {len(chains_set)}")
#     print(f"Number of unique residues: {len(residues_set)}")
#     print(f"Number of unique alpha carbons: {len(alpha_carbons_set)}")

#     return {
#         "models": models_set,
#         "chains": chains_set,
#         "residues": residues_set,
#         "alpha_carbons": alpha_carbons_set
#     }

# from pprint import pprint
# pprint(analysis_result)