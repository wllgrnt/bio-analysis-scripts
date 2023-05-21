"""
cov_and_edge_spots.py

Given the output of cellprofiler, we get a set of folders with All_measurements.csv, Expand_Nuclei.csv,
Nuclei.csv, Perinuclear_region.csv. We process all_measurements.csv using the count_total_and_mean_edge_spots
code to average the intensity of the edge spots over our vars (XY, wellnumber, and time) into either a single 2d
with the average for each field of view, or a set of 2d arrays in different files, depending on whether we have
2 or 3 features.

ALL_measurements: edge spots per cell
Expand nuclei: mass displacement of 60mer or mitochondria
Nuclei: nada
Pericentriolar region: CoV

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


def create_output_path(input_path: str, suffix: str):
    global OUTPUT_PATH
    return (
        OUTPUT_PATH
        + os.path.basename(os.path.dirname(input_path))
        + "_"
        + os.path.basename(input_path)[:-4]
        + suffix
    )


def save_pivot_table(data, values, index, columns, output_filename: str):
    pivot = pd.pivot_table(data=data, values=values, index=index, columns=columns)
    print("writing csv: ", output_filename)
    pivot.to_csv(output_filename, index=False)


def plot_edge_spot_pivots(input_path, aggregated_by_field_of_view):
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


def generate_edge_spot_pivots(
    input_path: str,
    val_variables: list,
    x_variable: str,
    y_variable: str,
    filename_variable: str | None,
):
    """generate edge spot data from the input_path.

    We process the input_path into an x, y, (f), val dataframe, then aggregate and pivot.
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

    output_path = create_output_path(input_path, "_raw.csv")
    processed_df.to_csv(output_path, index=False)

    output_path = create_output_path(input_path, ".csv")
    aggregated_by_field_of_view.to_csv(output_path, index=False)

    for val_variable in val_variables:
        if filename_variable is not None:
            for filename_var in aggregated_by_field_of_view[filename_variable].unique():
                filename_var_vals = aggregated_by_field_of_view.query(
                    f"{filename_variable} == @filename_var"
                )
                pivot_filename = create_output_path(
                    input_path, f"_{val_variable}_{filename_variable}{filename_var}.csv"
                )
                save_pivot_table(
                    filename_var_vals,
                    val_variable,
                    x_variable,
                    y_variable,
                    pivot_filename,
                )
        else:
            pivot_filename = create_output_path(input_path, f"_{val_variable}.csv")
            save_pivot_table(
                aggregated_by_field_of_view,
                val_variable,
                x_variable,
                y_variable,
                pivot_filename,
            )

    if SHOW_FIG:
        plot_edge_spot_pivots(input_path, aggregated_by_field_of_view)


def generate_mass_displacement_files(
    input_path: str,
    val_variables: list,
    x_variable: str,
    y_variable: str,
    filename_variable: str | None,
):
    """generate mass displacement files (not pivots at all, they aren't square) from input_path.

    TODO refactor common logic.
    """
    global VARIABLES
    cellprofiler_df = pd.read_csv(input_path)
    # convert to T, XY, WellNumber, vals
    use_cols = []
    if "T" in VARIABLES:
        cellprofiler_df["T"] = cellprofiler_df[FILENAME_COLUMN].apply(extract_timestamp)
        use_cols.append("T")
    if "XY" in VARIABLES:
        cellprofiler_df["XY"] = cellprofiler_df[FILENAME_COLUMN].apply(extract_xy)
        use_cols.append("XY")
    if "WellNumber" in VARIABLES:
        cellprofiler_df["WellNumber"] = cellprofiler_df[FILENAME_COLUMN].apply(
            extract_wellnumber
        )
        use_cols.append("WellNumber")
    use_cols.extend(val_variables)
    processed_df = cellprofiler_df[use_cols]

    aggregated_by_field_of_view = (
        processed_df.groupby(VARIABLES)
        .agg(**{f"{val_var}_mean": (val_var, "mean") for val_var in val_variables})
        .reset_index()
    )

    # save intermediates
    output_path = create_output_path(input_path, "_raw.csv")
    processed_df.to_csv(output_path, index=False)
    output_path = create_output_path(input_path, ".csv")
    aggregated_by_field_of_view.to_csv(output_path, index=False)
    # pivot
    pivot_variables = [f"{val_var}_mean" for val_var in val_variables]
    for pivot_var in pivot_variables:
        if filename_variable is not None:
            for filename_var in aggregated_by_field_of_view[filename_variable].unique():
                filename_var_vals = aggregated_by_field_of_view.query(
                    f"{filename_variable} == @filename_var"
                )
                pivot_filename = create_output_path(
                    input_path, f"_{pivot_var}_{filename_variable}{filename_var}.csv"
                )
                save_pivot_table(
                    filename_var_vals, pivot_var, x_variable, y_variable, pivot_filename
                )
        else:
            pivot_filename = create_output_path(input_path, f"_{pivot_var}.csv")
            save_pivot_table(
                aggregated_by_field_of_view,
                pivot_var,
                x_variable,
                y_variable,
                pivot_filename,
            )


def generate_cov_files(
    input_path: str, x_variable: str, y_variable: str, filename_variable: str | None
):
    """generate cov files from input_path."""
    global VARIABLES
    cellprofiler_df = pd.read_csv(input_path)
    # convert to T, XY, WellNumber, vals
    use_cols = []
    if "T" in VARIABLES:
        cellprofiler_df["T"] = cellprofiler_df[FILENAME_COLUMN].apply(extract_timestamp)
        use_cols.append("T")
    if "XY" in VARIABLES:
        cellprofiler_df["XY"] = cellprofiler_df[FILENAME_COLUMN].apply(extract_xy)
        use_cols.append("XY")
    if "WellNumber" in VARIABLES:
        cellprofiler_df["WellNumber"] = cellprofiler_df[FILENAME_COLUMN].apply(
            extract_wellnumber
        )
        use_cols.append("WellNumber")

    # generate the CoV
    STD_COLUMN = "Intensity_StdIntensity_MIRO160mer"
    MEAN_COLUMN = "Intensity_MeanIntensity_MIRO160mer"
    cellprofiler_df["CoV"] = cellprofiler_df[STD_COLUMN] / cellprofiler_df[MEAN_COLUMN]
    use_cols.extend(["CoV"])
    processed_df = cellprofiler_df[use_cols]
    # print(cellprofiler_df)
    print(processed_df.head())
    val_variables = [["CoV"]]
    aggregated_by_field_of_view = (
        processed_df.groupby(VARIABLES)
        .agg(**{f"{val_var}_mean": (val_var, "mean") for val_var in val_variables})
        .reset_index()
    )

    # save intermediates
    output_path = create_output_path(input_path, "_raw.csv")
    processed_df.to_csv(output_path, index=False)
    output_path = create_output_path(input_path, ".csv")
    aggregated_by_field_of_view.to_csv(output_path, index=False)
    print(aggregated_by_field_of_view.head())
    # pivot
    pivot_variables = [f"{val_var}_mean" for val_var in val_variables]
    for pivot_var in pivot_variables:
        if filename_variable is not None:
            for filename_var in aggregated_by_field_of_view[filename_variable].unique():
                filename_var_vals = aggregated_by_field_of_view.query(
                    f"{filename_variable} == @filename_var"
                )
                pivot_filename = create_output_path(
                    input_path, f"_{pivot_var}_{filename_variable}{filename_var}.csv"
                )
                save_pivot_table(
                    filename_var_vals, pivot_var, x_variable, y_variable, pivot_filename
                )
        else:
            pivot_filename = create_output_path(input_path, f"_{pivot_var}.csv")
            save_pivot_table(
                aggregated_by_field_of_view,
                pivot_var,
                x_variable,
                y_variable,
                pivot_filename,
            )


if __name__ == "__main__":
    # constants
    INPUT_FOLDERS = glob.glob("input_files/230517")
    OUTPUT_PATH = "output_files/230517/"
    SHOW_FIG = False
    VARIABLES = ["XY", "WellNumber"]
    for input_folder in INPUT_FOLDERS:
        # edge spot fractions from All_mesurements.csv
        all_measurements_path = os.path.join(input_folder, "All_measurements.csv")
        VAL_VARIABLES = ["edge_spot_mean_intensity_mean", "edge_spot_fraction"]
        X_VARIABLE = "XY"
        Y_VARIABLE = "WellNumber"
        FILENAME_VARIABLE = None
        generate_edge_spot_pivots(
            all_measurements_path,
            val_variables=VAL_VARIABLES,
            x_variable=X_VARIABLE,
            y_variable=Y_VARIABLE,
            filename_variable=FILENAME_VARIABLE,
        )
        # # mass displacement from Expand_Nuclei.csv
        # expand_nuclei_path = os.path.join(input_folder, "Expand_Nuclei.csv")
        # MITO_DISPLACEMENT = "Intensity_MassDisplacement_mito"
        # CAGE_DISPLACEMENT = "Intensity_MassDisplacement_MIRO160mer"
        # X_VARIABLE = "XY"
        # Y_VARIABLE = "WellNumber"
        # FILENAME_VARIABLE = None
        # generate_mass_displacement_pivots(
        #     expand_nuclei_path,
        #     val_variables=[MITO_DISPLACEMENT, CAGE_DISPLACEMENT],
        #     x_variable=X_VARIABLE,
        #     y_variable=Y_VARIABLE,
        #     filename_variable=FILENAME_VARIABLE,
        # )
        # CoV from Perinuclear_region.csv
        # peri_path = os.path.join(input_folder, "Perinuclear_region.csv")
        # generate the CoV from the standard deviation and mean.
        # VAL_VARIABLE = "CoV"
        # X_VARIABLE = "XY"
        # Y_VARIABLE = "WellNumber"
        # FILENAME_VARIABLE = None
        # generate_cov_pivots(
        #     peri_path,
        #     x_variable=X_VARIABLE,
        #     y_variable=Y_VARIABLE,
        #     filename_variable=FILENAME_VARIABLE,
        # )
