"""
Run the premade cellprofiler pipeline in batches over the images provided.
"""

import cellprofiler_core.pipeline
import cellprofiler_core.preferences
import cellprofiler_core.utilities.java
import cellprofiler.modules
import pathlib

OUTPUT_SETTING_NUMBER = 8

# Populate the strings below
PIPELINE_PATH = ""
INPUT_PATH = ""
OUTPUT_PATH = ""

cellprofiler_core.preferences.set_headless()

# Start the Java VM
cellprofiler_core.utilities.java.start_java()

# Read the pipeline
pipeline = cellprofiler_core.pipeline.Pipeline()
pipeline.load(PIPELINE_PATH)

# set the output directory (this is quite brittle)
export_module = pipeline.modules()[-1]
assert isinstance(
    export_module, cellprofiler.modules.exporttospreadsheet.ExportToSpreadsheet
)
export_file_directory = export_module.setting(OUTPUT_SETTING_NUMBER)
assert (
    export_file_directory.to_dict()["text"] == "Output file location"
), export_file_directory.to_dict()["text"]
export_file_directory.set_value(f"Elsewhere...|{OUTPUT_PATH}")

# Load the images
file_list = list(pathlib.Path(".").absolute().glob(f"{INPUT_PATH}/*.tif"))
assert file_list  # check that the glob found some files
files = [file.as_uri() for file in file_list]
pipeline.read_file_list(files)

# Run
output_measurements = pipeline.run()

print(output_measurements.get_measurement_columns())

# Need to stop the Java VM otherwise it will hang
cellprofiler_core.utilities.java.stop_java()
