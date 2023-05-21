"""
nd2_file_analyser.py


For all nd2 files in a given folder:

1. Read the nd2 file. It is a set of 2d images specified by channels, times, z coords and position (in the well)
2. Max-project along the z dimension
3. Adjust intensity with visual inspection ?!
4. Export a set of .tiff files for CellProfiler, or run CP-equivalent analysis in Python.
"""

import napari
import nd2
import os


from glob import glob

FOLDER_PATH = os.getcwd()

print(f"reading files in {FOLDER_PATH}")
for nd2_filepath in glob(f"{FOLDER_PATH}/*.nd2"):
    print("processing", nd2_filepath)
    with nd2.ND2File(nd2_filepath) as nd2_file:
        print(nd2_file.sizes)
        print(nd2_file.shape)

        print(dir(nd2))
        print(dir(nd2_file))
        nd2_array = nd2_file.to_dask()
        # # max project
        zmax = nd2_array.max(
            axis=2
        )  # run some asserts and checks to make sure 2 is the correct axis.
        print(dir(nd2_array))

        print(zmax.shape)

        # viewer = napari.view_image(zmax, channel_axis=2, ndisplay=3)
        # napari.run()
