"""
count_number_of_edge_spots.py

Given a csv file, extract the well number and xy from the file name, and count
the number of nuclei and 'edge spots' in each xy.

"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

# set up your variable names

sns.set_style("whitegrid")

FILENAME_COLUMN = "FileName_Hoechst"  # from which we extract the well number and xy
NUCLEI_COUNT_COLUMN = "ObjectNumber"
EDGE_SPOT_COUNT_COLUMN = "Number_Object_Number"
INPUT_PATHS = [
    "input_files/C06/All_measurements.csv",
    "input_files/E06/All_measurements.csv",
    "input_files/B03/All_measurements.csv",
    "input_files/E03/All_measurements.csv",
]


# read the csv
def extract_xy(input_string):
    return int(re.search(r"(?<=_XY)\d+", input_string).group(0))


def extract_well_number(input_string):
    return int(re.search(r"(?<=_T)\d+", input_string).group(0))


for input_path in INPUT_PATHS:
    # header=0,1 would get you a multiindex but we can drop for now
    cellprofiler_output_df = pd.read_csv(input_path, header=[0, 1])
    # cellprofiler_output_df = pd.read_csv(input_path, header=1)

    # print(cellprofiler_output_df['Image'])
    # print(cellprofiler_output_df['Image'].columns)
    # print(cellprofiler_output_df['Image']['Filename_Hoechst'])
    # extract the XY and T values from the filename
    cellprofiler_output_df["XY"] = cellprofiler_output_df["Image"][
        FILENAME_COLUMN
    ].apply(extract_xy)
    cellprofiler_output_df["T"] = cellprofiler_output_df["Image"][
        FILENAME_COLUMN
    ].apply(extract_well_number)
    cellprofiler_output_df["nuclei_count"] = cellprofiler_output_df["Nuclei"][
        NUCLEI_COUNT_COLUMN
    ]
    cellprofiler_output_df["edge_spot_count"] = cellprofiler_output_df["edge_spots"][
        EDGE_SPOT_COUNT_COLUMN
    ]
    # group by XY and well number, and count the number of nuclei and edge spots
    cellprofiler_output_df = cellprofiler_output_df[
        ["T", "XY", "nuclei_count", "edge_spot_count"]
    ].reset_index(drop=True)
    # cellprofiler_output_df.columns = ['well_number', 'XY', 'nuclei_count', 'edge_spot_count']

    counts = cellprofiler_output_df.groupby(["T", "XY"]).count().reset_index()
    counts["edge_spot_fraction"] = counts["edge_spot_count"] / counts["nuclei_count"]

    print(counts)
    fig, ax = plt.subplots()
    counts.plot.scatter(x="T", y="edge_spot_fraction", ax=ax)
    ax.set_title(input_path)
    plt.show()
