
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import sys

input_folder = "input_folder/colabfold_json"
output_folder = "output_folder"
plot = True

path_to_special_lines = {"2ca5d": [369, 474, 487, 729, 834, 847],
                         "6d779": [369, 474, 483, 729, 834, 843]}


# the above code reveals that pae is a square array. Output to a csv.
for filename in os.listdir(input_folder):
    if filename.endswith(".json"):
        file_path = os.path.join(input_folder, filename)
        with open(file_path, "r") as f:
            colabfold_data = json.load(f)
        
        pae_array = np.array(colabfold_data['pae'])
        csv_filename = os.path.splitext(filename)[0] + "_pae.csv"
        output_path = os.path.join(output_folder, csv_filename)
        print(pae_array.shape)
        if plot:
            fig, ax = plt.subplots(figsize=(10,10))
            cmap = sns.color_palette("coolwarm", as_cmap=True)
            # cmap = sns.diverging_palette(220, 20, l=65, as_cmap=True)

            sns.heatmap(pae_array,
                         cmap=cmap, 
                         xticklabels=False,
                         yticklabels=False,
                         square=True,
                         vmin=0,
                         vmax=30,
                         cbar=False
                           )
            
            ax.vlines([143, 286, 646, 1006], *ax.get_ylim(), colors='black')
            ax.hlines([143, 286, 646, 1006], *ax.get_xlim(), colors='black')
            
            for subpath, extra_lines in path_to_special_lines.items():
                if subpath in file_path:
                    ax.vlines(extra_lines, *ax.get_ylim(), colors='black', linestyles='dotted')
                    ax.hlines(extra_lines, *ax.get_xlim(), colors='black', linestyles='dotted')

            output_image_path = os.path.join(output_folder, os.path.splitext(filename)[0] + "_pae.png") 
            plt.savefig(output_image_path, dpi=300, bbox_inches='tight', pad_inches=0)
            plt.close(fig)

            # separate cbar fig
            fig_cmap, ax_cmap = plt.subplots(figsize=(1, 5))
            plt.imshow(np.arange(30).reshape((30, 1)), cmap=cmap)
            ax_cmap.set_axis_off()

            output_cmap_path = os.path.join(output_folder, os.path.splitext(filename)[0] + "_cmap.png") 
            plt.savefig(output_cmap_path, dpi=300, bbox_inches='tight', pad_inches=0)
            plt.close(fig_cmap)  # Close the colormap figure

        np.savetxt(output_path, pae_array, delimiter=",", fmt='%.2f')
