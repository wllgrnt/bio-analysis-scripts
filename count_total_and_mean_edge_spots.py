"""
count_number_of_edge_spots.py

Given a csv file, extract the well number and xy from the file name, and count
the number of nuclei and 'edge spots' in each time/wellnumber/xy. Also grab
the edge spot median intensity, and average (mean + median) over time/wellnumber/xy.
"""

import glob
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

sns.set_style("whitegrid")

FILENAME_COLUMN = "FileName_Hoechst"  # from which we extract the well number and xy
NUCLEI_COUNT_COLUMN = "ObjectNumber"
EDGE_SPOT_COUNT_COLUMN = "Number_Object_Number"
EDGE_SPOT_MEDIAN_INTENSITY_COLUMN = (
    "Intensity_MedianIntensity_Subtract_perinucleus_Miro160mer"
)
EDGE_SPOT_MEAN_INTENSITY_COLUMN = (
    "Intensity_MeanIntensity_Subtract_perinucleus_Miro160mer"
)
INPUT_PATHS = glob.glob("input_files/for_will/*/All_measurements.csv")
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
    input_df["edge_spot_median_intensity"] = input_df["edge_spots"][
        EDGE_SPOT_MEDIAN_INTENSITY_COLUMN
    ]
    input_df["edge_spot_mean_intensity"] = input_df["edge_spots"][
        EDGE_SPOT_MEAN_INTENSITY_COLUMN
    ]
    use_cols.extend(
        [
            "nuclei_count",
            "edge_spot_count",
            "edge_spot_median_intensity",
            "edge_spot_mean_intensity",
        ]
    )
    processed_df = input_df[use_cols].reset_index(drop=True)
    # drop the bottom of the multiindex to prevent future headaches
    processed_df.columns = processed_df.columns.droplevel(1)
    return processed_df


for input_path in INPUT_PATHS:
    print(input_path)
    cellprofiler_df = pd.read_csv(input_path, header=[0, 1])

    processed_df = process_cellprofiler_df(
        cellprofiler_df,
        get_timestamp="T" in VARIABLES,
        get_xy="XY" in VARIABLES,
        get_wellnumber="WellNumber" in VARIABLES,
    )
    # count the non-nan nuclei and edge spots, and mean/median the intensity
    vals_per_fov = (
        processed_df.groupby(VARIABLES)
        .agg(
            edge_spot_count=("edge_spot_count", "count"),
            nuclei_count=("nuclei_count", "count"),
            edge_spot_median_intensity_median=("edge_spot_median_intensity", "median"),
            edge_spot_median_intensity_mean=("edge_spot_median_intensity", "mean"),
            edge_spot_mean_intensity_median=("edge_spot_mean_intensity", "median"),
            edge_spot_mean_intensity_mean=("edge_spot_mean_intensity", "mean"),
        )
        .reset_index()
    )

    # vals_per_fov = processed_df.groupby(VARIABLES).count().reset_index()
    vals_per_fov["edge_spot_fraction"] = (
        vals_per_fov["edge_spot_count"] / vals_per_fov["nuclei_count"]
    )

    # plot the results
    fig, axes = plt.subplots(nrows=2, figsize=(7, 10))
    vals_per_fov.plot.scatter(x=VARIABLES[0], y="edge_spot_fraction", ax=axes[0])
    vals_per_fov.plot.scatter(
        x=VARIABLES[0], y="edge_spot_median_intensity_median", ax=axes[1]
    )
    for ax in axes:
        ax.set_title(input_path)
    # get the name of the input file without its parent folder
    output_path = (
        OUTPUT_PATH
        + os.path.basename(os.path.dirname(input_path))
        + "_"
        + os.path.basename(input_path)[:-4]
        + ".png"
    )
    plt.tight_layout()
    plt.savefig(output_path)
    processed_df.to_csv(output_path[:-4] + "_raw.csv", index=False)
    vals_per_fov.to_csv(output_path[:-4] + ".csv", index=False)
    if SHOW_FIG:
        plt.show()
