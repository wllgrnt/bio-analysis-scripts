
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns

input_folder = "input_folder/colabfold_json"

# for filename in os.listdir(input_folder):
#     if filename.endswith(".json"):
#         file_path = os.path.join(input_folder, filename)
#         with open(file_path, "r") as f:
#             colabfold_data = json.load(f)
        
#         print(f"Keys in {filename}:")
#         for key in colabfold_data.keys():
#             print(key)
#             print(type(colabfold_data[key]))
#             if type(colabfold_data[key]) is list:
#                 arr = np.array(colabfold_data[key])
#                 print(arr.shape)
#                 if len(arr.shape) == 2:
#                     fig, ax = plt.subplots(figsize=(10,10))
#                     sns.heatmap(arr)
#                     plt.show()
#         print()


# the above code reveals that pae is a square array. Output to a csv.
for filename in os.listdir(input_folder):
    if filename.endswith(".json"):
        file_path = os.path.join(input_folder, filename)
        with open(file_path, "r") as f:
            colabfold_data = json.load(f)
        
        pae_array = np.array(colabfold_data['pae'])
        csv_filename = os.path.splitext(filename)[0] + "_pae.csv"
        np.savetxt(csv_filename, pae_array, delimiter=",")
