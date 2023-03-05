"""
Run the premade cellprofiler pipeline in batches over the images provided.

Cellprofiler needs python 3.8 to run.
"""

import cellprofiler_core.pipeline
import cellprofiler_core.preferences
import cellprofiler_core.utilities.java
import cellprofiler.modules
import os
import pathlib

from multiprocessing import Pool
from typing import List

# Set params below
OUTPUT_SETTING_NUMBER = 8
NUM_WORKERS = 8
PIPELINE_PATH = ""
INPUT_FOLDERS: List[str] = ["", ""]


def run_pipeline(pipeline_path: str, input_folder: str, output_folder: str) -> None:
    """Within an env with cellprofiler and java running, run the pipeline.

    Args:
        pipeline_path: path to a cppipe file
        input_folder: a folder containing tif files
        output_folder: a folder to save the output to
    Returns:
        None, because the cellprofiler_core.measurement._measurements.Measurements object
            can't be pickled.
    """
    print(
        f"Running pipeline {pipeline_path} on {input_folder} and saving to {output_folder}"
    )
    # Start cellprofiler
    cellprofiler_core.preferences.set_headless()
    cellprofiler_core.utilities.java.start_java()
    # Read the pipeline
    pipeline = cellprofiler_core.pipeline.Pipeline()
    pipeline.load(pipeline_path)

    # set the output directory (this is quite brittle)
    export_module = pipeline.modules()[-1]
    assert isinstance(
        export_module, cellprofiler.modules.exporttospreadsheet.ExportToSpreadsheet
    )

    export_file_directory = export_module.setting(OUTPUT_SETTING_NUMBER)
    assert (
        export_file_directory.to_dict()["text"] == "Output file location"
    ), export_file_directory.to_dict()["text"]
    export_file_directory.set_value(f"Elsewhere...|{output_folder}")
    # Load the images
    file_list = list(pathlib.Path(".").absolute().glob(f"{input_folder}/*.tif"))
    assert file_list  # check that the glob found some files
    files = [file.as_uri() for file in file_list]
    pipeline.read_file_list(files)
    # Run, discard output
    _ = pipeline.run()
    # Need to stop the Java VM otherwise it will hang
    cellprofiler_core.utilities.java.stop_java()
    return None


if __name__ == "__main__":
    input_folders = [f"input_files/{x}" for x in INPUT_FOLDERS]
    for input_folder in input_folders:
        # for each subfolder in the input folder, run the pipeline.
        print("running analysis for: ", input_folder)
        subfolders = [f.path for f in os.scandir(input_folder) if f.is_dir()]
        with Pool(
            NUM_WORKERS
        ) as p:  # parallelise - pool is limited by number of subfolders
            p.starmap(
                run_pipeline,
                [
                    (
                        PIPELINE_PATH,
                        subfolder,
                        subfolder.replace("input_files", "output_files"),
                    )
                    for subfolder in subfolders
                ],
            )
