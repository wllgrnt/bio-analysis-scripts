"""
count_number_of_edge_spots.py

Given a csv file, extract the well number and xy from the file name, and count
the number of nuclei and 'edge spots' in each xy.

"""

import glob
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

sns.set_style("whitegrid")

FILENAME_COLUMN = "FileName_Hoechst"  # from which we extract the well number and xy
NUCLEI_COUNT_COLUMN = "Number_Object_Number"
EDGE_SPOT_COUNT_COLUMN = "Number_Object_Number"
INPUT_PATHS = glob.glob("input_files/*/All_measurements.csv")
OUTPUT_PATH = "output_files/"
SHOW_FIG = False
VARIABLES = ["WellNumber", "XY"]


def extract_xy(input_string: str) -> int:
    """(Christina is a criminal, and her xys are not consistently formatted)"""
    try:
        return int(re.search(r"(?<=_XY)\d+", input_string).group(0))
    except AttributeError:
        return int(re.search(r"(?<=_)\d{4}(?=_)", input_string).group(0))


def extract_timestamp(input_string: str) -> int:
    return int(re.search(r"(?<=_T)\d+", input_string).group(0))


def extract_wellnumber(input_string: str) -> str:
    return re.search(r"(?<=_Well)[A-Z]\d+", input_string).group(0)


def process_cellprofiler_df(
    input_df: pd.DataFrame, get_timestamp=False, get_xy=False, get_wellnumber=False
) -> pd.DataFrame:
    """Given a cellprofiler dataframe, extract the data from the filename column."""
    use_cols = []
    if get_timestamp:
        input_df["T"] = input_df["Image"][FILENAME_COLUMN].apply(extract_timestamp)
        use_cols.append("T")
    if get_xy:
        input_df["XY"] = input_df["Image"][FILENAME_COLUMN].apply(extract_xy)
        use_cols.append("XY")
    if get_wellnumber:
        input_df["WellNumber"] = input_df["Image"][FILENAME_COLUMN].apply(
            extract_wellnumber
        )
        use_cols.append("WellNumber")
    input_df["nuclei_count"] = input_df["Nuclei"][NUCLEI_COUNT_COLUMN]
    input_df["edge_spot_count"] = input_df["edge_spots"][EDGE_SPOT_COUNT_COLUMN]
    use_cols.extend(["nuclei_count", "edge_spot_count"])
    return input_df[use_cols].reset_index(drop=True)


for input_path in INPUT_PATHS:
    print(input_path)
    cellprofiler_df = pd.read_csv(input_path, header=[0, 1])

    processed_df = process_cellprofiler_df(
        cellprofiler_df,
        get_timestamp="T" in VARIABLES,
        get_xy="XY" in VARIABLES,
        get_wellnumber="WellNumber" in VARIABLES,
    )

    counts = processed_df.groupby(VARIABLES).count().reset_index()
    counts["edge_spot_fraction"] = counts["edge_spot_count"] / counts["nuclei_count"]

    # plot the results
    fig, ax = plt.subplots()
    counts.plot.scatter(x=VARIABLES[0], y="edge_spot_fraction", ax=ax)
    ax.set_title(input_path)
    # get the name of the input file without its parent folder
    output_path = (
        OUTPUT_PATH
        + os.path.basename(os.path.dirname(input_path))
        + "_"
        + os.path.basename(input_path)[:-4]
        + ".png"
    )
    plt.savefig(output_path)
    counts.to_csv(output_path[:-4] + ".csv", index=False)
    if SHOW_FIG:
        plt.show()
