"""
Utilities for working with cellprofiler outputs.
"""

import re
import pandas as pd


FILENAME_COLUMN = "FileName_Hoechst"  # from which we extract the well number and xy
NUCLEI_COUNT_COLUMN = "Number_Object_Number"
EDGE_SPOT_COUNT_COLUMN = "Number_Object_Number"
EDGE_SPOT_MEDIAN_INTENSITY_COLUMN = (
    "Intensity_MedianIntensity_Subtract_perinucleus_Miro160mer"
)
EDGE_SPOT_MEAN_INTENSITY_COLUMN = (
    "Intensity_MeanIntensity_Subtract_perinucleus_Miro160mer"
)


def extract_xy(input_string: str) -> int:
    try:
        return int(re.search(r"(?<=_XY)\d+", input_string).group(0))
    except AttributeError:
        return int(re.search(r"(?<=_)\d{4}(?=_)", input_string).group(0))


def extract_timestamp(input_string: str) -> int:
    return int(re.search(r"(?<=_T)\d+", input_string).group(0))


def extract_wellnumber(input_string: str) -> str:
    return re.search(r"(?<=Well)[A-Z]\d+", input_string).group(0)


def process_cellprofiler_edgespot_df(
    input_df: pd.DataFrame,
    get_timestamp=False,
    get_xy=False,
    get_wellnumber=False,
    filename_column=FILENAME_COLUMN,
    nuclei_count_column=NUCLEI_COUNT_COLUMN,
    edge_spot_count_column=EDGE_SPOT_COUNT_COLUMN,
    edge_spot_median_intensity_column=EDGE_SPOT_MEDIAN_INTENSITY_COLUMN,
    edge_spot_mean_intensity_column=EDGE_SPOT_MEAN_INTENSITY_COLUMN,
) -> pd.DataFrame:
    """Given a cellprofiler dataframe, extract the data from the filename column."""
    use_cols = []
    if get_timestamp:
        input_df["T"] = input_df["Image"][filename_column].apply(extract_timestamp)
        use_cols.append("T")
    if get_xy:
        input_df["XY"] = input_df["Image"][filename_column].apply(extract_xy)
        use_cols.append("XY")
    if get_wellnumber:
        input_df["WellNumber"] = input_df["Image"][filename_column].apply(
            extract_wellnumber
        )
        use_cols.append("WellNumber")
    input_df["nuclei_count"] = input_df["Nuclei"][nuclei_count_column]
    input_df["edge_spot_count"] = input_df["edge_spots"][edge_spot_count_column]
    input_df["edge_spot_median_intensity"] = input_df["edge_spots"][
        edge_spot_median_intensity_column
    ]
    input_df["edge_spot_mean_intensity"] = input_df["edge_spots"][
        edge_spot_mean_intensity_column
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


def aggregate_edgespot_df(input_df: pd.DataFrame, groupby_vars) -> pd.DataFrame:
    """Given a cleaned df, aggregate over <input_df>

    the df must contain groupby_vars + ["edge_spot_count", "nuclei_count", "edge_spot_median_intensity", "edge_spot_mean_intensity"]

    """

    output_df = (
        input_df.groupby(groupby_vars)
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
