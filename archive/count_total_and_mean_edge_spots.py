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
NUCLEI_COUNT_COLUMN = "Number_Object_Number"
EDGE_SPOT_COUNT_COLUMN = "Number_Object_Number"
EDGE_SPOT_MEDIAN_INTENSITY_COLUMN = (
    "Intensity_MedianIntensity_Subtract_perinucleus_Miro160mer"
)
EDGE_SPOT_MEAN_INTENSITY_COLUMN = (
    "Intensity_MeanIntensity_Subtract_perinucleus_Miro160mer"
)
INPUT_PATHS = glob.glob("input_files/for_will_7/All_measurements.csv")
OUTPUT_PATH = "output_files/"
SHOW_FIG = False
VARIABLES = ["XY", "WellNumber"]


def extract_xy(input_string: str) -> int:
    """(Christina is a criminal, and her xys are not consistently formatted)"""
    try:
        return int(re.search(r"(?<=_XY)\d+", input_string).group(0))
    except AttributeError:
        return int(re.search(r"(?<=_)\d{4}(?=_)", input_string).group(0))


def extract_timestamp(input_string: str) -> int:
    return int(re.search(r"(?<=_T)\d+", input_string).group(0))


def extract_wellnumber(input_string: str) -> str:
    return re.search(r"(?<=Well)[A-Z]\d+", input_string).group(0)


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


def aggregate_df(input_df: pd.DataFrame) -> pd.DataFrame:
    """Given a cleaned df, aggregate over VARIABLES."""
    global VARIABLES

    output_df = (
        input_df.groupby(VARIABLES)
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

    output_df["edge_spot_fraction"] = (
        output_df["edge_spot_count"] / output_df["nuclei_count"]
    )

    return output_df


if __name__ == "__main__":
    for input_path in INPUT_PATHS:
        cellprofiler_df = pd.read_csv(input_path, header=[0, 1])

        processed_df = process_cellprofiler_df(
            cellprofiler_df,
            get_timestamp="T" in VARIABLES,
            get_xy="XY" in VARIABLES,
            get_wellnumber="WellNumber" in VARIABLES,
        )

        aggregated_by_field_of_view = aggregate_df(processed_df)

        # get the name of the input file without its parent folder
        output_path = (
            OUTPUT_PATH
            + os.path.basename(os.path.dirname(input_path))
            + "_"
            + os.path.basename(input_path)[:-4]
        )
        processed_df.to_csv(output_path + "_raw.csv", index=False)
        aggregated_by_field_of_view.to_csv(output_path + ".csv", index=False)

        # make pivot tables (there is probably a fancier way than the unique() iterator)
        VAL_VARIABLES = ["edge_spot_mean_intensity_mean", "edge_spot_fraction"]
        X_VARIABLE = "XY"
        Y_VARIABLE = "WellNumber"
        FILENAME_VARIABLE = None

        if FILENAME_VARIABLE is not None:
            for filename_var in aggregated_by_field_of_view[FILENAME_VARIABLE].unique():
                filename_var_vals = aggregated_by_field_of_view.query(
                    f"{FILENAME_VARIABLE} == @filename_var"
                )
                for val_variable in VAL_VARIABLES:
                    pivot = pd.pivot_table(
                        data=filename_var_vals,
                        values=val_variable,
                        index=X_VARIABLE,
                        columns=Y_VARIABLE,
                    )
                    pivot_filename = (
                        output_path
                        + f"_{val_variable}_{FILENAME_VARIABLE}{filename_var}.csv"
                    )
                    pivot.to_csv(pivot_filename, index=False)

        else:
            for val_variable in VAL_VARIABLES:
                pivot = pd.pivot_table(
                    data=aggregated_by_field_of_view,
                    values=val_variable,
                    index=X_VARIABLE,
                    columns=Y_VARIABLE,
                )
                pivot_filename = output_path + f"_{val_variable}.csv"
                pivot.to_csv(pivot_filename, index=False)

        if SHOW_FIG:
            # plot the results
            fig, axes = plt.subplots(nrows=2, figsize=(7, 10))
            aggregated_by_field_of_view.plot.scatter(
                x=VARIABLES[0], y="edge_spot_fraction", ax=axes[0]
            )
            aggregated_by_field_of_view.plot.scatter(
                x=VARIABLES[0], y="edge_spot_median_intensity_median", ax=axes[1]
            )
            for ax in axes:
                ax.set_title(input_path)
            plt.tight_layout()
            plt.show()
