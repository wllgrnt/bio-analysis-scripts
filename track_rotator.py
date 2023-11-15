"""
- input data structure:
    - folder for each of three conditions: no_TRAK_77 / TRAK1_79 / TRAK2_78
    - each folder has a subfolder for each cell
    - each cell subfolder has ‘tracks.csv’, ‘edges.csv’, ‘spots.csv’
- output
    - the number of backwards and forwards tracks per cell (to be able to plot final stats)
    - a graph plotting the tracks (to check that vector is correct)
    - a sh∫eets of speeds for ‘anterograde’ and ‘retrograde’ movement in each condition
    - raw: TRACK_ID header with speeds in columns for each track
    - stacked: all TRACK_ID speeds stacked in a single column (for histogram)
"""

# %%
import matplotlib.pyplot as plt
import numpy as np
import os
import polars as pl
import pandas as pd
import seaborn as sns


max_distance_travelled_filter = 2.5 # exclude tracks travelling less than this distance
path_to_tracks = './input_folder/tracks/'
path_to_cell_orientations= './input_folder/tracks/cell_orientation_coordinates.xlsx'

max_distance_travelled_filter = 2.5 # exclude tracks travelling less than this distance
# %%
cell_orientations = pl.read_excel(path_to_cell_orientations, sheet_name='Angles_per_cell')


def zero_positions(df):
    df = df.sort('FRAME')
    return df.with_columns(POSITION_X_ZEROED = df['POSITION_X'] - df['POSITION_X'][0],
                    POSITION_Y_ZEROED = df['POSITION_Y'] - df['POSITION_Y'][0],
                    POSITION_T_ZEROED = df['POSITION_T'] - df['POSITION_T'][0])


def get_rotation_angle(condition, cell):
    angle_df = cell_orientations.filter((pl.col('Cell') == cell) & (pl.col('Condition') == condition))
    assert len(angle_df) == 1
    return angle_df['angle (rad)'][0]


if __name__ == '__main__':
    for condition in os.listdir(path_to_tracks):
        condition_path = os.path.join(path_to_tracks, condition)
        if not os.path.isdir(condition_path):
            continue

        for cell in os.listdir(condition_path):
            
            cell_path = os.path.join(condition_path, cell)
            for file in ['tracks.csv', 'edges.csv', 'spots.csv']:
                assert os.path.exists(os.path.join(cell_path, file))

            # From tracks.csv, we need the max distance travelled.
            tracks_df = pl.read_csv(os.path.join(cell_path, 'tracks.csv'), skip_rows_after_header=3)
            # assert tracks_df['TRACK_INDEX'][0] == 1
            assert 'TRACK_ID' in tracks_df.columns
            spots_df = pl.read_csv(os.path.join(cell_path, 'spots.csv'), skip_rows_after_header=3)
            
            spots_df = spots_df[['TRACK_ID', 'POSITION_X', 'POSITION_Y', 'POSITION_Z', 'POSITION_T', 'FRAME']]
            tracks_df = tracks_df[['TRACK_ID', 'MAX_DISTANCE_TRAVELED']]
            spots_with_max_distance = spots_df.join(tracks_df, how='left', on='TRACK_ID', validate='m:1')
            spots_with_max_distance_filtered = spots_with_max_distance.filter(pl.col('MAX_DISTANCE_TRAVELED') > 2.5).sort(by=['TRACK_ID', 'FRAME'])
            if not len(spots_with_max_distance_filtered):
                print(f'warning!, no points found for cell {cell}, condition {condition}')
                continue
            spots_zeroed = spots_with_max_distance_filtered.group_by('TRACK_ID').map_groups(zero_positions)
            # rotation_angle_rads = 1.5

            rotation_angle_rads = get_rotation_angle(condition, cell)
            print(condition, cell, rotation_angle_rads)
            spots_rotated = spots_zeroed.with_columns(POSITION_X_ROTATED = pl.col('POSITION_X_ZEROED') * np.sin(rotation_angle_rads) -  pl.col('POSITION_Y_ZEROED') * np.cos(rotation_angle_rads),
                                    POSITION_Y_ROTATED = pl.col('POSITION_X_ZEROED') * np.cos(rotation_angle_rads) +  pl.col('POSITION_Y_ZEROED') * np.sin(rotation_angle_rads))

            fig, ax = plt.subplots()
            for track_id in spots_rotated['TRACK_ID'].unique():
                single_track = spots_rotated.filter(pl.col('TRACK_ID') == track_id)
                plt.plot(single_track['POSITION_X_ROTATED'], single_track['POSITION_Y_ROTATED']) #, c='TRACK_ID')
            plt.axis('equal')
            plt.show()

# %%
