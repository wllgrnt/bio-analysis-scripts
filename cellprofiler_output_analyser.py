"""
Given a folder with the following structure:

INPUT_FOLDER
    - subfolder1
      - All_measurements.csv
      - Expand_nuclei.csv
      - Perinuclear_region.csv

In which either wellnumber and XY vary, or wellnumber, XY, and T vary, do the following:

- from All_measurements, extract the number of edge spots as a fraction of the number of nuclei per XY.
  Generate intermediate with Well, XY, T, and fraction, one row per (well, XY, T) combination.
  Then:

  If T doesn't vary, generate a file with:

    |  Well_1  | Well_2 | Well_3   ....
    | -------------------------
XY1 |  fraction | fraction | fraction
XY2 ....

 If T varies, generate a file for each well number with:

     |  T1      |     T2   |    T3 ....
     | -------------------------
XY1  | fraction | fraction | fraction
XY2  ....

Then average over XYs and plot:

    |    Well_1      |     Well_2    | Well_3   ....
    | -------------------------
T1  |  mean_fraction | mean_fraction | mean_fraction
...

- from Expand_nuclei, extract the mass displacement of the 60mer/mito for each cell. As this is a per-cell measure,
  we don't really care about XY, but track anyway. Generate intermediate with Well, XY, T, and
  mass displacement, where there are N rows per (well, XY, T) combination, where N is the number
  of cells in that field of view. Then:

  If T doesn't vary, generate a file with:

    |   Well_1  |  Well_2   | Well_3   ....
    | XY1 | XY2 | XY1 | XY2 | XY1 | XY2
    | -------------------------
    |  mass | mass | mass | mass | mass | mass
    |  mass | mass | mass | mass | mass | mass
    ...


And
   |  Well_1 |  Well_2   | Well_3   ....
   | -------------------------
   | mass    |   mass    |   mass

with Well_1 contains all XY1 values, stacked.


NB the above is not square since each FoV may have a different number of cells. Fill with NaNs.

If T varies, generate a file for each well number with:

    | T1        |    T2     |    T3 ....
    | XY1 | XY2 | XY1 | XY2 | XY1 | XY2
    | ------------
    |  mass | mass | mass | mass | mass | mass
    ....

Plus the stacked variant as above.

Then average over values and generate:

    | Well_1 | Well_2 | Well_3 ...

 T1 | mean_mass | mean_mass | mean_mass
 ...

- from Perinuclear_region, extract the CoV for each cell, and then generate files as for the mass
  displacement.


Script parameters:
    - input_folder: path to the folder containing the input files
    - output_folder: path to the folder where the output files will be written
    - t_varies: boolean, whether the T parameter varies or not
    - plot: boolean, whether to plot the output.

"""
import glob
import logging
import os
import re
from numpy import average

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")

INPUT_FOLDER = "input_folder"  # the path to all the input files
INPUT_SUBFOLDER = (
    "analysis"  # the path to the specific files we're analysing right now.
)
OUTPUT_FOLDER = "output_folder"  # output will be saved to OUTPUT_FOLDER/INPUT_SUBFOLDER
T_VARIES = True
EDGE_SPOT_FILE = "All_measurements.csv"
LOGGING_LEVEL = logging.INFO  # logging.INFO or logging.ERROR  normally
PLOT = False

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def generate_edge_spot_files(
    input_path: str, output_folder: str, t_varies: bool, do_plot=True
):
    """
    Generate the edge_spot_count / nuclei_count for each field of view.

    Assumes two header columns.

    Args:
        - input_path: the path to the input file. Must have the following columns:
        - output_folder: where to output the processed files.
        - t_varies. If we have multiple timepoints and must therefore process the file differently.
            We check this with asserts.
    """
    logger.info(
        f"Generating edge spot data for {input_path}, output to {output_folder}, t_varies={t_varies}"
    )
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    raw_input_df = pd.read_csv(input_path, header=[0, 1])

    processed_df = extract_edgespot_cols(raw_input_df, t_varies)
    intermediate_filepath = os.path.join(output_folder, "edge_spot_fraction_raw.csv")
    logger.info(f"Writing edge spot intermediate to {intermediate_filepath}")
    processed_df.to_csv(intermediate_filepath, index=False)

    if t_varies:
        # File for each well number, T as columns, XY as rows
        for well_number in processed_df.WellNumber.unique():
            well_number_subdf = processed_df[processed_df.WellNumber == well_number]
            output_filename = os.path.join(
                output_folder, f"edge_spot_fraction_over_time_{well_number}.csv"
            )
            save_pivot_table(
                data=well_number_subdf,
                values="edge_spot_fraction",
                index="XY",
                columns="T",
                output_filename=output_filename,
            )

        # average over XYs, normalise and plot:
        average_over_time = (
            processed_df.groupby(["WellNumber", "T"])
            .edge_spot_fraction.mean()
            .reset_index()
        )
        pivot = pd.pivot_table(
            data=average_over_time,
            values="edge_spot_fraction",
            index="T",
            columns="WellNumber",
        )
        pivot = pivot.divide(pivot.iloc[0])
        output_path = os.path.join(
            output_folder, "edge_spot_fraction_mean_over_time_normalised.csv"
        )
        logger.info(
            f"Normalised pivot table has shape {pivot.shape}, writing to {output_path}"
        )
        pivot.to_csv(output_path)
        if do_plot:
            pivot.plot(title="Edge spot fraction over time, normalised to T0")
            plt.show()

    else:
        # single file - well number vs xy.
        output_filename = os.path.join(output_folder, "edge_spot_fraction_static.csv")
        save_pivot_table(
            data=processed_df,
            values="edge_spot_fraction",
            index="XY",
            columns="WellNumber",
            output_filename=output_filename,
        )


def extract_edgespot_cols(cellprofiler_df: pd.DataFrame, t_varies) -> pd.DataFrame:
    """
    Extract the (well_number, xy, t) columns and aggregate over W, XY, T
    to generate the edge spot data.

    Makes a number of assumptions about the struture of the input df.

    """
    logger.info(f"Extracting columns, input df has shape {cellprofiler_df.shape}")
    filename_column = "FileName_Hoechst"  # from which we extract the well number and xy
    nuclei_count_column = "Number_Object_Number"
    edge_spot_count_column = "Number_Object_Number"

    output_df = cellprofiler_df.copy()
    if t_varies:
        output_df["T"] = output_df["Image"][filename_column].apply(extract_timestamp)
    output_df["XY"] = output_df["Image"][filename_column].apply(extract_xy)
    output_df["WellNumber"] = output_df["Image"][filename_column].apply(
        extract_wellnumber
    )
    output_df["nuclei_count"] = output_df["Nuclei"][nuclei_count_column]
    output_df["edge_spot_count"] = output_df["edge_spots"][edge_spot_count_column]
    cols = ["WellNumber", "XY", "nuclei_count", "edge_spot_count"]
    if t_varies:
        cols.append("T")
    output_df = output_df[cols].reset_index(drop=True)
    output_df.columns = output_df.columns.droplevel(
        1
    )  # drop the second level of the multiindex

    aggregate_df = (
        output_df.groupby(
            ["WellNumber", "XY", "T"] if t_varies else ["WellNumber", "XY"]
        )
        .agg({"nuclei_count": "count", "edge_spot_count": "count"})
        .reset_index()
    )
    aggregate_df["edge_spot_fraction"] = (
        aggregate_df["edge_spot_count"] / aggregate_df["nuclei_count"]
    )

    logger.info(
        f"Extracted columns and aggregated, output df has shape {aggregate_df.shape}"
    )
    return aggregate_df


def extract_xy(input_string: str) -> int:
    try:
        return int(re.search(r"(?<=_XY)\d+", input_string).group(0))
    except AttributeError:
        return int(re.search(r"(?<=_)\d{4}(?=_)", input_string).group(0))


def extract_timestamp(input_string: str) -> int:
    return int(re.search(r"(?<=_T)\d+", input_string).group(0))


def extract_wellnumber(input_string: str) -> str:
    return re.search(r"(?<=Well)[A-Z]\d+", input_string).group(0)


def save_pivot_table(data, values, index, columns, output_filename: str):
    pivot = pd.pivot_table(data=data, values=values, index=index, columns=columns)
    logger.info(f"writing pivot table: {output_filename}")
    pivot.to_csv(output_filename)


if __name__ == "__main__":
    input_folders = glob.glob(os.path.join(INPUT_FOLDER, INPUT_SUBFOLDER, "*/"))
    for input_folder in input_folders:
        logger.info(f"Processing {input_folder}")
        edge_spot_file_path = os.path.join(input_folder, EDGE_SPOT_FILE)
        output_subfolder = input_folder.replace(INPUT_FOLDER, OUTPUT_FOLDER)
        generate_edge_spot_files(edge_spot_file_path, output_subfolder, T_VARIES, PLOT)
        logger.info("\n")
