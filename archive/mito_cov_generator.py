import pandas as pd  # for reading/writing and table manipulation
import seaborn as sns  # for plotting
import matplotlib.pyplot as plt  # also for plotting
import glob  # for wildcard matching

from count_total_and_mean_edge_spots import *

# set up your variable names

# INPUT_PATH = 'tests/Perinuclear_region.csv'
# INPUT_PATH = '2DG_mito/5mM2DGNOwash/Expand_Nuclei.csv'
# OUTPUT_PATH = 'output.csv'
DO_PLOT = False
SLICE_INDEX_START = 45  # 13
SLICE_INDEX_STOP = 48  # 16
IMPORTANT_COLUMNS = [
    "Intensity_MassDisplacement_MIRO160mer",
    "Intensity_MassDisplacement_mito",
]

INPUT_PATHS = glob.glob("input_files/for_will_8/Expand_Nuclei.csv")

stitched_dict = {x: {} for x in IMPORTANT_COLUMNS}

FILENAME_COLUMN = "FileName_MIRO160mer"  # from which we extract the well number and xy
STD_COLUMN = "Intensity_StdIntensity_MIRO160mer"
MEAN_COLUMN = "Intensity_MeanIntensity_MIRO160mer"

for input_path in INPUT_PATHS:
    # read the csv
    cellprofiler_output_df = pd.read_csv(input_path)

    # make the CoV
    cellprofiler_output_df["CoV"] = (
        cellprofiler_output_df[STD_COLUMN] / cellprofiler_output_df[MEAN_COLUMN]
    )

    # grab the well number
    cellprofiler_output_df["well_number"] = cellprofiler_output_df[
        FILENAME_COLUMN
    ].apply(
        extract_wellnumber  # noqa: F405
    )

    print(cellprofiler_output_df)

    for important_col in IMPORTANT_COLUMNS:
        # for each well number, take the values of <important_column>, and put in a dict
        new_array_dict = {}
        median_dict = {}
        for well_number in cellprofiler_output_df[
            "well_number"
        ].unique():  # iterate over well numbers
            # grab the section of the dataframe corresponding to the well number
            slice_of_array = cellprofiler_output_df[
                cellprofiler_output_df.well_number == well_number
            ]
            # grab the column we want from that section
            column_slice = slice_of_array[important_col]
            # do pandas magic with index so we can build a new dataframe at the end
            column_slice = column_slice.reset_index(drop=True)
            # assign to our dictionary
            new_array_dict[well_number] = column_slice
            median_dict[well_number] = column_slice.median()
        # generate new dataframe
        pivot_df = pd.DataFrame.from_dict(new_array_dict)

        # relabel the median dictionary and normalise by the first value
        median_dict_normalised = {}
        print(median_dict.keys())
        for key, value in median_dict.items():
            minute = (int(key[1:]) - 1) * 40  # this 15 should be changed
            median_dict_normalised[minute] = median_dict[key] / median_dict["T01"]

        # put in our stitched dictionary
        condition = input_path.split("/")[1]
        stitched_dict[important_col][condition] = median_dict_normalised

        # write to output_path
        output_path = input_path[:-4] + "_" + important_col + "_transform.csv"
        pivot_df.to_csv(output_path)

        if DO_PLOT:
            fig, ax = plt.subplots(figsize=(20, 10))
            sns.boxplot(data=pivot_df)
            ax.set_ylabel(important_col)
            ax.set_xlabel("well_number")
            ax.set_title(input_path)
            plt.savefig(output_path[:-4] + ".png")
            # plt.show()


# convert our big dict into a dataframe so we can write it to csv
for column_name, array_dictionary in stitched_dict.items():
    print(column_name)
    df = pd.DataFrame(array_dictionary)
    print(df.head())
    df.to_csv(column_name + ".csv")


# repeat the above for CoV
IMPORTANT_COLUMNS = ["CoV"]

stitched_dict = {x: {} for x in IMPORTANT_COLUMNS}

for input_path in glob.glob("230118/*/Perinuclear_region.csv"):
    # read the csv
    cellprofiler_output_df = pd.read_csv(input_path)

    # make the CoV
    cellprofiler_output_df["CoV"] = (
        cellprofiler_output_df["Intensity_StdIntensity_MIRO160mer"]
        / cellprofiler_output_df["Intensity_MeanIntensity_MIRO160mer"]
    )

    # grab the well number
    cellprofiler_output_df["well_number"] = cellprofiler_output_df[
        "FileName_MIRO160mer"
    ].str.slice(SLICE_INDEX_START, SLICE_INDEX_STOP)

    for important_col in IMPORTANT_COLUMNS:
        # for each well number, take the values of <important_column>, and put in a dict
        new_array_dict = {}
        median_dict = {}
        for well_number in cellprofiler_output_df[
            "well_number"
        ].unique():  # iterate over well numbers
            # grab the section of the dataframe corresponding to the well number
            slice_of_array = cellprofiler_output_df[
                cellprofiler_output_df.well_number == well_number
            ]
            # grab the column we want from that section
            column_slice = slice_of_array[important_col]
            # do pandas magic with index so we can build a new dataframe at the end
            column_slice = column_slice.reset_index(drop=True)
            # assign to our dictionary
            new_array_dict[well_number] = column_slice
            median_dict[well_number] = column_slice.median()
        # generate new dataframe
        pivot_df = pd.DataFrame.from_dict(new_array_dict)

        # relabel the median dictionary and normalise by the first value
        median_dict_normalised = {}
        for key, value in median_dict.items():
            minute = (int(key[1:]) - 1) * 40  # this 15 should be changed
            median_dict_normalised[minute] = median_dict[key] / median_dict["T01"]

        # put in our stitched dictionary
        condition = input_path.split("/")[1]
        stitched_dict[important_col][condition] = median_dict_normalised

        # write to output_path
        output_path = input_path[:-4] + "_" + important_col + "_transform.csv"
        pivot_df.to_csv(output_path)

        # pivot_df

        if DO_PLOT:
            fig, ax = plt.subplots(figsize=(20, 10))
            sns.boxplot(data=pivot_df)
            ax.set_ylabel(important_col)
            ax.set_xlabel("well_number")
            ax.set_title(input_path)
            plt.savefig(output_path[:-4] + ".png")
            # plt.show()


# print(stitched_dict)

for column_name, array_dictionary in stitched_dict.items():
    print(column_name)
    df = pd.DataFrame(array_dictionary)
    print(df.head())
    df.to_csv(column_name + ".csv")
