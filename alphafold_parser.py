import gemmi
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import Dark2
import os

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
                    std_b_factors = np.std(b_factors)
                    chain_info[chain_name]["residue_info"][residue_number] = (avg_b_factor, std_b_factors)

    return chain_info

def plot_b_factors_by_chain(chain_info, output_folder, base_filename):
    # cmap = get_cmap('Dark2')
    colors = Dark2.colors

    chain_names = list(chain_info.keys())
    num_chains = len(chain_names)
    
    fig, axes = plt.subplots(2, num_chains, figsize=(15, 10), sharey='row')

    for idx, chain_name in enumerate(chain_names):
        data = chain_info[chain_name]
        residues = sorted(data["residue_info"].keys())
        avg_b_factors = [data["residue_info"][res][0] for res in residues]
        std_b_factors = [data["residue_info"][res][1] for res in residues]

        alpha_carbon_residues = [ac[0] for ac in data["alpha_carbons"]]
        alpha_carbon_b_factors = [ac[1] for ac in data["alpha_carbons"]]

        color = colors[idx % len(colors)]

        # Plot average B factors with standard deviation band
        axes[0, idx].plot(residues, avg_b_factors, linewidth=2, label='Average B Factor', color=color)
        axes[0, idx].fill_between(residues, 
                                  np.array(avg_b_factors) - np.array(std_b_factors), 
                                  np.array(avg_b_factors) + np.array(std_b_factors), 
                                  color=color, alpha=0.2, label='Standard Deviation')
        axes[0, idx].set_title(f'Chain {chain_name}')
        axes[0, idx].set_xlabel('Residue Number')
        if idx == 0:
            axes[0, idx].set_ylabel('B Factor')
        axes[0, idx].legend()

        # Plot alpha carbon B factors
        axes[1, idx].plot(alpha_carbon_residues, alpha_carbon_b_factors, linewidth=2, color=color)
        axes[1, idx].set_xlabel('Residue Number')
        if idx == 0:
            axes[1, idx].set_ylabel('Alpha Carbon B Factor')

    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"{base_filename}.png"))
    plt.close()

def output_csv_by_chain(chain_info, output_folder, base_filename):
    for chain_name, data in chain_info.items():
        csv_file = os.path.join(output_folder, f"{base_filename}_{chain_name}.csv")
        with open(csv_file, "w") as f:
            f.write("residue_number,b_factor\n")
            for residue_number, b_factor in data["alpha_carbons"]:
                f.write(f"{residue_number},{b_factor}\n")

def process_cif_files(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    cif_files = [f for f in os.listdir(input_folder) if f.endswith('.cif')]
    
    for cif_file in cif_files:
        cif_path = os.path.join(input_folder, cif_file)
        base_filename = os.path.splitext(cif_file)[0]

        chain_info = extract_cif_info_by_chain(cif_path)
        output_csv_by_chain(chain_info, output_folder, base_filename)
        plot_b_factors_by_chain(chain_info, output_folder, base_filename)

# Example usage
input_folder = "input_folder/cifs"
output_folder = "output_folder/cifs"

process_cif_files(input_folder, output_folder)
