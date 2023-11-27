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


def zero_positions(df):
    df = df.sort('FRAME')
    return df.with_columns(POSITION_X_ZEROED = df['POSITION_X'] - df['POSITION_X'][0],
                    POSITION_Y_ZEROED = df['POSITION_Y'] - df['POSITION_Y'][0],
                    POSITION_T_ZEROED = df['POSITION_T'] - df['POSITION_T'][0])

def get_rotation_angle(condition, cell, path_to_cell_orientations):
    cell_orientations = pl.read_excel(path_to_cell_orientations, sheet_name='angles_per_cell')
    angle_df = cell_orientations.filter((pl.col('Cell') == cell) & (pl.col('Condition') == condition))
    if len(angle_df) != 1:
        return None
    else:
       return angle_df['angle (rad)'][0]

def generate_rotated_spot_positions(cell_path: str, path_to_cell_orientations) -> pd.DataFrame | None:
    for file in ['tracks.csv', 'edges.csv', 'spots.csv']:
        assert os.path.exists(os.path.join(cell_path, file))

    # From tracks.csv, we need the max distance travelled.
    tracks_df = pl.read_csv(os.path.join(cell_path, 'tracks.csv'), skip_rows_after_header=3)
    # assert tracks_df['TRACK_INDEX'][0] == 1
    assert 'TRACK_ID' in tracks_df.columns
    spots_df = pl.read_csv(os.path.join(cell_path, 'spots.csv'), skip_rows_after_header=3)
    
    # edges_df = pl.read_csv(os.path.join(cell_path, 'edges.csv'), skip_rows_after_header=3)
    # assert 'TRACK_ID' in edges_df.columns


    spots_df = spots_df[['TRACK_ID', 'POSITION_X', 'POSITION_Y', 'POSITION_Z', 'POSITION_T', 'FRAME']]
    tracks_df = tracks_df[['TRACK_ID', 'MAX_DISTANCE_TRAVELED',
                            'TRACK_MEAN_SPEED',
                            'TRACK_MAX_SPEED',
                            'TRACK_MIN_SPEED',
                            'TRACK_MEDIAN_SPEED',
                            'TRACK_STD_SPEED',
                            'MEAN_STRAIGHT_LINE_SPEED',
                            'MEAN_DIRECTIONAL_CHANGE_RATE']]
    spots_with_max_distance = spots_df.join(tracks_df, how='left', on='TRACK_ID', validate='m:1')
    spots_with_max_distance_filtered = spots_with_max_distance.sort(by=['TRACK_ID', 'FRAME'])
    if not len(spots_with_max_distance_filtered):
        print(f'warning!, no points found for cell {cell}, condition {condition}')
        return None
    spots_zeroed = spots_with_max_distance_filtered.group_by('TRACK_ID').map_groups(zero_positions)
    rotation_angle_rads = get_rotation_angle(condition, cell, path_to_cell_orientations)
    if rotation_angle_rads is None:
        print(f'warning!, no angle found for cell {cell}, condition {condition}')
        return None
    print(condition, cell, rotation_angle_rads)
    spots_rotated = spots_zeroed.with_columns(POSITION_X_ROTATED = pl.col('POSITION_X_ZEROED') * np.sin(rotation_angle_rads) -  pl.col('POSITION_Y_ZEROED') * np.cos(rotation_angle_rads),
                            POSITION_Y_ROTATED = pl.col('POSITION_X_ZEROED') * np.cos(rotation_angle_rads) +  pl.col('POSITION_Y_ZEROED') * np.sin(rotation_angle_rads))

    return spots_rotated  # ignore polars for now

def zero_positions(df):
    df = df.sort('FRAME')
    return df.with_columns(POSITION_X_ZEROED = df['POSITION_X'] - df['POSITION_X'][0],
                    POSITION_Y_ZEROED = df['POSITION_Y'] - df['POSITION_Y'][0],
                    POSITION_T_ZEROED = df['POSITION_T'] - df['POSITION_T'][0])


def get_rotation_angle(condition, cell, path_to_cell_orientations, sheet_name='Sheet1'):
    cell_orientations = pl.read_excel(path_to_cell_orientations, sheet_name=sheet_name)
    angle_df = cell_orientations.filter((pl.col('cell') == cell) & (pl.col('condition') == condition))
    if len(angle_df) != 1:
        return None
    else:
       return angle_df['angle'][0]



# if __name__ == '__main__':

dfs = {}
edge_dfs = {}

# generate the positions over time for each condition, cell, track. plot.
for date in ['231027', '231102', '231103', '231117']:
    path_to_tracks = f'./input_folder/all_tracks/{date}'
    path_to_cell_orientations= f'{path_to_tracks}/{date}_cell_orientation_coordinates.xlsx'
    path_to_plots = f'plots/{date}'
    if not os.path.exists(path_to_plots):
        os.mkdir(path_to_plots)
    for condition in os.listdir(path_to_tracks):
        condition_path = os.path.join(path_to_tracks, condition)
        if not os.path.isdir(condition_path):
            continue

        for cell in os.listdir(condition_path):
            cell_path = os.path.join(condition_path, cell)
            spots_rotated =  generate_rotated_spot_positions(cell_path, path_to_cell_orientations)
            # spots_rotated = spots_rotated.filter(pl.col('MAX_DISTANCE_TRAVELED') > distance_filter)
            dfs[(date, condition, cell)]  = spots_rotated.to_pandas()  # my polars is insufficient for what follows.    
            edge_dfs[(date, condition, cell)] = pl.read_csv(os.path.join(cell_path, 'edges.csv'), skip_rows_after_header=3).to_pandas() 
            fig, ax = plt.subplots()
            for track_id in spots_rotated['TRACK_ID'].unique():
                single_track = spots_rotated.filter(pl.col('TRACK_ID') == track_id)
                plt.plot(single_track['POSITION_X_ROTATED'], single_track['POSITION_Y_ROTATED']) #, c='TRACK_ID')
            plt.axis('equal')
            plt.title(f'{date}: {condition}_{cell}')

            plt.savefig(f'{path_to_plots}/{condition}_{cell}.png')
            plt.close()
            # plt.show()
            # if condition == 'no_TRAK_77':
            #     plt.show()
            
# %%
# collate the tracks for each cell into one big dataframe.s
all_points_df = pd.concat(dfs, keys=dfs.keys()).reset_index(names=['date', 'condition', 'cell', 'index']).drop(columns=['index'])
# information on the tracks, one row per track 
per_track = all_points_df.groupby(['date', 'condition', 'cell', 'TRACK_ID']).agg(
    max_distance_traveled=('MAX_DISTANCE_TRAVELED','first'),
    track_mean_speed=('TRACK_MEAN_SPEED','first'),
    track_max_speed=('TRACK_MAX_SPEED','first'),
    track_min_speed=('TRACK_MIN_SPEED','first'),
    track_median_speed=('TRACK_MEDIAN_SPEED','first'),
    track_std_speed=('TRACK_STD_SPEED','first'),
    mean_straight_line_speed=('MEAN_STRAIGHT_LINE_SPEED','first'),
    mean_directional_change_rate=('MEAN_DIRECTIONAL_CHANGE_RATE','first'),
    final_y_position=('POSITION_Y_ROTATED', 'last'),
).reset_index()


# %%

# do as above but with edges

all_edges_df = pd.concat(edge_dfs, keys=edge_dfs.keys()).reset_index(names=['date', 'condition', 'cell', 'index']).drop(columns=['index'])

all_edges_df_with_track_info = all_edges_df.merge(per_track, how='left', on=['date', 'condition', 'cell', 'TRACK_ID'], validate='m:1')

# %%

distance_filter = 2.8
all_edges_df_with_track_info_filtered = all_edges_df_with_track_info[all_edges_df_with_track_info['max_distance_traveled'] > distance_filter]

all_edges_df_with_track_info_filtered['final_y_is_above_zero'] = all_edges_df_with_track_info_filtered['final_y_position'] > 0


sns.violinplot(hue='condition', y='SPEED', x='final_y_is_above_zero', data=all_edges_df_with_track_info_filtered[all_edges_df_with_track_info_filtered['condition'] != 'no_TRAK_77'], cut=0)
# fig, ax = plt.subplots(figsize=(10, 10))

# sns.violinplot(x='condition', y='final_y_position', hue='final_y_is_above_zero', data=per_track_filtered, cut=0)
# plt.axhline(y=0)


# for col in ['track_mean_speed', 'track_max_speed', 'track_min_speed', 'track_median_speed', 'track_std_speed', 'mean_straight_line_speed', 'mean_directional_change_rate']:
#     fig, ax = plt.subplots(figsize=(10, 10))
#     sns.violinplot(x='condition', y=col, hue='final_y_is_above_zero', data=all_edges_df_with_track_info_filtered[all_edges_df_with_track_info_filtered['condition'] != 'no_TRAK_77'], cut=0)
#     plt.title(col)
#     plt.show()
#     # plt.savefig(f'plots/{col}.png')
#     # plt.close()







# %%




# # %%
distance_filter = 2.8
per_track_filtered = per_track[per_track['max_distance_traveled'] > distance_filter]

per_track_filtered['final_y_is_above_zero'] = per_track_filtered['final_y_position'] > 0

# fig, ax = plt.subplots(figsize=(10, 10))

# sns.violinplot(x='condition', y='final_y_position', hue='final_y_is_above_zero', data=per_track_filtered, cut=0)
# plt.axhline(y=0)


# for col in ['track_mean_speed', 'track_max_speed', 'track_min_speed', 'track_median_speed', 'track_std_speed', 'mean_straight_line_speed', 'mean_directional_change_rate']:
for col in ['track_median_speed']:
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.violinplot(hue='condition', y=col, x='final_y_is_above_zero', data=per_track_filtered[per_track_filtered['condition'] != 'no_TRAK_77'], cut=0)
    plt.title(col)
    plt.show()
    # plt.savefig(f'plots/{col}.png')
    # plt.close()
# %%




# # %%
# # Add instantaneous speed to the dataframe
# unique_id_cols = ['date', 'condition', 'cell', 'TRACK_ID']
# all_points_df['dx'] = all_points_df.groupby(unique_id_cols)['POSITION_X_ZEROED'].diff()
# all_points_df['dy'] = all_points_df.groupby(unique_id_cols)['POSITION_Y_ZEROED'].diff()
# all_points_df['dt'] = all_points_df.groupby(unique_id_cols)['POSITION_T_ZEROED'].diff()

# all_points_df['distance'] = np.sqrt(all_points_df['dx']**2 + all_points_df['dy']**2)
# all_points_df['speed'] = all_points_df['distance'] / all_points_df['dt']
# # %%
# all_points_df['angle'] = np.arctan2(all_points_df['dy'], all_points_df['dx'])
# all_points_df['angle_difference'] = all_points_df.groupby(unique_id_cols)['angle'].diff()

# # %%
# all_points_df['angle_degrees'] = np.degrees(all_points_df['angle'])
# all_points_df['angle_difference_degrees'] = all_points_df.groupby(unique_id_cols)['angle_degrees'].diff()

# # %%
# g = sns.displot(all_points_df, x='angle_degrees', hue='condition', kind='kde', fill=True)
# g.set(xticks=range(-180, 181, 45))
# plt.show()
# # %%
# cell_path
edge_df = pl.read_csv('./input_folder/all_tracks/231117/no_TRAK_77/231117_03/edges.csv', skip_rows_after_header=3)

# # %%
# sns.displot(edge_df, x='SPEED', kind='kde', fill=True, cut=0)
# # %%
# sns.displot(all_points_df, x='speed', kind='kde', fill=True, cut=0)

# %%
