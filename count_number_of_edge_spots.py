"""
count_number_of_edge_spots.py

Given a csv file, extract the well number and xy from the file name, and count
the number of nuclei and 'edge spots' in each xy.

"""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

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


def extract_xy(input_string):
    """Christina is a criminal, and her xys are not consistently formatted."""
    try:
        return int(re.search(r"(?<=_XY)\d+", input_string).group(0))
    except AttributeError:
        return int(re.search(r"(?<=_)\d{4}(?=_)", input_string).group(0))


def extract_well_number(input_string):
    return int(re.search(r"(?<=_T)\d+", input_string).group(0))


for input_path in INPUT_PATHS:
    cellprofiler_output_df = pd.read_csv(input_path, header=[0, 1])
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
    counts = cellprofiler_output_df.groupby(["T", "XY"]).count().reset_index()
    counts["edge_spot_fraction"] = counts["edge_spot_count"] / counts["nuclei_count"]

    # plot the results
    fig, ax = plt.subplots()
    counts.plot.scatter(x="T", y="edge_spot_fraction", ax=ax)
    ax.set_title(input_path)
    # get the name of the input file without its parent folder
    output_path = (
        "output_files/"
        + os.path.basename(os.path.dirname(input_path))
        + "_"
        + os.path.basename(input_path)[:-4]
        + ".png"
    )
    plt.savefig(output_path)
    counts.to_csv(output_path[:-4] + ".csv", index=False)
    # plt.show()
