""" This file contains a function that merges the subtile output of the shade calculation into one map"""
import os
import glob
import rasterio
import numpy as np
from tqdm import tqdm
import time

def merge_tif_files_by_time(main_folder, output_folder, merged_name, start_time=900, end_time=2000, delete_input_files=False,
                            nodata_value=-9999):
    """
    Merge multiple TIF files from specified subfolders into a single mosaic file for each specified time step.
    -------
    Input:
    - main_folder (str): Path to the main folder containing subfolders with TIF files to be merged.
    - output_folder (str): Path to the folder where the merged TIF files will be saved.
    - merged_name (str): Base name for the merged output files. Each file will be suffixed with the date and time.
    - nodata_value (int, optional): Value to represent 'no data' in the output TIF files. Default is -9999.
    - start_time (int, optional): Starting time for filtering TIF files based on their filenames. Default is 900.
    - end_time (int, optional): Ending time for filtering TIF files based on their filenames. Default is 2000.
    - delete_input_files (bool, optional): If True, delete the original TIFF files after merging. Default is False.

    Output:
    - none: The function saves merged TIFF files directly to the specified `output_folder` with the naming format
      `{merged_name}_{date_str}_{time}.TIF`.
    """
    os.makedirs(output_folder, exist_ok=True)
    tile_folders = os.listdir(main_folder)
    # Store all files for each time
    files_by_time = {}
    merged_files = []

    # Iterate through all tile folders
    for tile_folder in tqdm(tile_folders, desc="Processing Tile Folders"):
        tile_path = os.path.join(main_folder, tile_folder)

        if not os.path.isdir(tile_path):
            continue

        tif_files = glob.glob(os.path.join(tile_path, '*.TIF'))

        for tif_file in tif_files:
            file_name = os.path.basename(tif_file)
            parts = file_name.split('_')

            # Skip the agglomerated file
            if len(parts) < 5:
                continue

            time_str = parts[4]
            time = int(time_str)
            date_str = parts[3]

            if start_time <= time <= end_time:
                if time not in files_by_time:
                    files_by_time[time] = []
                files_by_time[time].append(tif_file)

    # Merge files for each time step
    for time, files in tqdm(files_by_time.items(), desc="Merging by Time", leave=False):
        if not files:
            continue

        src_files_to_mosaic = []
        bounds = []
        resolutions = []

        # Open all files and save bounds
        for file in files:
            src = rasterio.open(file)
            src_files_to_mosaic.append(src)
            bounds.append(src.bounds)
            resolutions.append(src.res[0])  # pixels should be square

        # Calculating merged bounds
        min_x = min(b[0] for b in bounds)
        min_y = min(b[1] for b in bounds)
        max_x = max(b[2] for b in bounds)
        max_y = max(b[3] for b in bounds)

        # Output shape based on the merged bounds and resolution
        out_shape = (
        1, int((max_y - min_y) / resolutions[0]), int((max_x - min_x) / resolutions[0]))  # (bands, height, width)
        out_transform = rasterio.transform.from_bounds(min_x, min_y, max_x, max_y, out_shape[2], out_shape[1])

        # Initialize the mosaic
        mosaic = np.full(out_shape, nodata_value, dtype=np.float32)

        for src in src_files_to_mosaic:
            data = src.read(1).astype(np.float32)  # Read the first band

            # Getting start position of columns and rows
            col_start, row_start = ~out_transform * (src.bounds.left, src.bounds.top)
            row_start, col_start = int(round(row_start)), int(round(col_start))

            # Ensure indices are within bounds
            if (0 <= row_start < mosaic.shape[1] and 0 <= col_start < mosaic.shape[2]):

                window_height, window_width = data.shape
                end_row = min(row_start + window_height, mosaic.shape[1])
                end_col = min(col_start + window_width, mosaic.shape[2])

                data_slice = data[:end_row - row_start, :end_col - col_start]

                # Update the mosaic based on conditions: replace if no data, replace if incoming value is lower & not nodata
                mosaic[0, row_start:end_row, col_start:end_col] = np.where(
                    (mosaic[0, row_start:end_row, col_start:end_col] == nodata_value) |
                    ((data_slice < mosaic[0, row_start:end_row, col_start:end_col]) &
                     (data_slice != nodata_value)),
                    data_slice,
                    mosaic[0, row_start:end_row, col_start:end_col]
                )
            else:
                print(f"Warning: Data for {src.name} is outside the bounds of the mosaic, skipping.")

        # Replace any remaining NaN values with nodata_value
        mosaic = np.nan_to_num(mosaic, nan=nodata_value)

        # Update metadata for the output file
        out_meta = src_files_to_mosaic[0].meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": mosaic.shape[1],
            "width": mosaic.shape[2],
            "transform": out_transform,
            "count": mosaic.shape[0],
            "nodata": nodata_value
        })

        # Create output file name
        output_file = os.path.join(output_folder, f"{merged_name}_{date_str}_{time}.TIF")
        with rasterio.open(output_file, "w", **out_meta) as dest:
            dest.write(mosaic)

        merged_files.append(os.path.basename(output_file))
        print(f"Merged {len(files)} TIF files for time {time} into {output_file}")

        for src in src_files_to_mosaic:
            src.close()

    # Optional: delete input files
    if delete_input_files:
        print(f"Deleting original TIF files in {main_folder}, excluding merged files")
        for tile_folder in tile_folders:
            tile_path = os.path.join(main_folder, tile_folder)
            if not os.path.isdir(tile_path):
                continue

            # Get all .TIF files in the current tile folder
            tif_files_to_delete = glob.glob(os.path.join(tile_path, '*.TIF'))
            for tif_file in tif_files_to_delete:
                if os.path.basename(tif_file) not in merged_files:  # Only delete if not merged
                    os.remove(tif_file)
"""
merge_tif_files_by_time("D:\Geomatics\output", "D:\Geomatics\correct_merged", "amsterdam")
"""