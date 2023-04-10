"""
cov_and_edge_spots.py

Given the output of cellprofiler, we get a set of folders with All_measurements.csv, Expand_Nuclei.csv,
Nuclei.csv, Perinuclear_region.csv. We process all_measurements.csv using the count_total_and_mean_edge_spots
code to average the intensity of the edge spots over our vars (XY, wellnumber, and time) into either a single 2d
with the average for each field of view, or a set of 2d arrays in different files, depending on whether we have
2 or 3 features.

We do a similar thing with the CoV, ending up with the normalised cov over time for each field of view.

This is getting too complicated to be one big script.
"""
# imports
import glob
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


from utils import (
    FILENAME_COLUMN,
    process_cellprofiler_edgespot_df,
    aggregate_edgespot_df,
    extract_timestamp,
    extract_xy,
    extract_wellnumber,
)

sns.set_style("whitegrid")

# constants
INPUT_FOLDERS = glob.glob("input_files/analysis/*/")
OUTPUT_PATH = "output_files/analysis/"
SHOW_FIG = False
VARIABLES = ["XY", "WellNumber"]


def generate_edge_spot_pivots(
    input_path: str,
    val_variables: list,
    x_variable: str,
    y_variable: str,
    filename_variable: str,
):
    """generate edge spot data from the All_measurements.csv file.

    Horrendous encapsulation here, should refactor.
    """
    global VARIABLES
    global OUTPUT_PATH
    global SHOW_FIG
    print(f"Processing {input_path}")
    cellprofiler_df = pd.read_csv(input_path, header=[0, 1])

    processed_df = process_cellprofiler_edgespot_df(
        cellprofiler_df,
        get_timestamp="T" in VARIABLES,
        get_xy="XY" in VARIABLES,
        get_wellnumber="WellNumber" in VARIABLES,
    )

    aggregated_by_field_of_view = aggregate_edgespot_df(processed_df, VARIABLES)

    output_path = (
        OUTPUT_PATH
        + os.path.basename(os.path.dirname(input_path))
        + "_"
        + os.path.basename(input_path)[:-4]
    )
    processed_df.to_csv(output_path + "_raw.csv", index=False)
    aggregated_by_field_of_view.to_csv(output_path + ".csv", index=False)

    if filename_variable is not None:
        for filename_var in aggregated_by_field_of_view[filename_variable].unique():
            filename_var_vals = aggregated_by_field_of_view.query(
                f"{filename_variable} == @filename_var"
            )
            for val_variable in val_variables:
                pivot = pd.pivot_table(
                    data=filename_var_vals,
                    values=val_variable,
                    index=x_variable,
                    columns=y_variable,
                )
                pivot_filename = (
                    output_path
                    + f"_{val_variable}_{filename_variable}{filename_var}.csv"
                )
                pivot.to_csv(pivot_filename, index=False)

    else:
        for val_variable in val_variables:
            pivot = pd.pivot_table(
                data=aggregated_by_field_of_view,
                values=val_variable,
                index=x_variable,
                columns=y_variable,
            )
            pivot_filename = output_path + f"_{val_variable}.csv"
            pivot.to_csv(pivot_filename, index=False)

    if SHOW_FIG:
        # plot the results
        _, axes = plt.subplots(nrows=2, figsize=(7, 10))
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


if __name__ == "__main__":
    for input_folder in INPUT_FOLDERS:
        # input_path = os.path.join(input_folder, "All_measurements.csv")
        VAL_VARIABLES = ["edge_spot_mean_intensity_mean", "edge_spot_fraction"]
        X_VARIABLE = "XY"
        Y_VARIABLE = "WellNumber"
        FILENAME_VARIABLE = None
        # generate_edge_spot_pivots(input_path,
        #                           val_variables=VAL_VARIABLES,
        #                           x_variable=X_VARIABLE,
        #                           y_variable=Y_VARIABLE,
        #                           filename_variable=FILENAME_VARIABLE)
        # generate CoV data from the Expand_Nuclei.csv file.

        STD_COLUMN = "Intensity_StdIntensity_MIRO160mer"
        MEAN_COLUMN = "Intensity_MeanIntensity_MIRO160mer"
        input_path = os.path.join(input_folder, "Expand_Nuclei.csv")
        print(input_path)
        mito_df = pd.read_csv(input_path)
        mito_df["T"] = mito_df[FILENAME_COLUMN].apply(extract_timestamp)
        mito_df["XY"] = mito_df[FILENAME_COLUMN].apply(extract_xy)
        mito_df["WellNumber"] = mito_df[FILENAME_COLUMN].apply(extract_wellnumber)
        print(mito_df.columns)
        mito_df["CoV"] = mito_df[STD_COLUMN] / mito_df[MEAN_COLUMN]
        mito_df = mito_df[["T", "XY", "WellNumber", "CoV"]]
        print(mito_df)
        # sodijfsoidf
