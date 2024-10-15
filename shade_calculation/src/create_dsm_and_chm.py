"""
This script contains multiple functions to create a DSM with ground and building heights, and to create the normalized
CHM. As input, the CHM calculated from the pointcloud and the AHN DSM & DTM are used.
"""

import os
import glob
import re
import geopandas as gpd
import numpy as np
import rasterio
from shapely.geometry import mapping
from rasterio.features import geometry_mask
from scipy.interpolate import NearestNDInterpolator
from functions import write_output
import tqdm
import startinpy
import time
from rasterio import Affine

def get_bbox(raster_paths):
    """
    Compute the overlapping bounding box of the CHM, DSM and DTM file.
    ----
    Input:
    - raster_paths (list of strings): paths to the CHM, DSM and DTM .TIF files

    Output:
    - rasterio.coords.BoundingBox(): bbox of the
    """
    bboxes = []
    for raster_path in raster_paths:
        with rasterio.open(raster_path) as src:
            bbox = src.bounds
            bboxes.append(bbox)
    left = max([bbox.left for bbox in bboxes])
    bottom = max([bbox.bottom for bbox in bboxes])
    right = min([bbox.right for bbox in bboxes])
    top = min([bbox.top for bbox in bboxes])
    return rasterio.coords.BoundingBox(left, bottom, right, top)


def crop_raster(raster_path, bbox, no_data=-9999):
    """
    Crop the input rasters to the size of the overlapping bounding box.
    ----
    Input:
    - raster_path (string): paths to the .TIF file.
    - bbox (4-tuple): overlapping bounding box of the CHM, DSM and DTM file.
    - no_data (int, optional): no_data value to replace source no data value with.

    Output:
    - cropped_data (2d numpy array): cropped raster data.
    - src.window_transform(window): affine transform matrix for given window.
    - src.crs (rasterio src): A PROJ4 dict representation of the CRS of the input raster.
    """
    with rasterio.open(raster_path) as src:
        window = src.window(bbox.left, bbox.bottom, bbox.right, bbox.top)

        # Read the data for the specified window
        cropped_data = src.read(window=window)

        # Ensure window attributes are integers (sometimes bbox in float)
        row_off = int(window.row_off)
        col_off = int(window.col_off)
        height = int(window.height)
        width = int(window.width)

        # Replace no-data values in the data
        if src.nodata is not None:
            cropped_data[cropped_data == src.nodata] = no_data

        return cropped_data, src.window_transform(window), src.crs


def extract_center_cells(cropped_data, no_data=-9999):
    """
    Extract the values of each cell in the input data and save these with the x and y (row and col)
    indices. Thereby, make sure that the corners of the dataset are filled for a full coverage triangulation
    in the next step.
    ----
    Input:
    - cropped_data (2d numpy array): cropped raster data.
    - no_data (int, optional): no_data value to replace source no data value with.

    Output:
    - xyz_filled (list): list containing x, y and z coordinates of the cells.
    """
    # Get coordinates of center cells
    cropped_data = cropped_data[0, :, :]

    # Get the indices of the rows and columns
    rows, cols = np.indices(cropped_data.shape)

    # Identify corner coordinates
    corners = {
        "top_left": (0, 0),
        "top_right": (0, cropped_data.shape[1] - 1),
        "bottom_left": (cropped_data.shape[0] - 1, 0),
        "bottom_right": (cropped_data.shape[0] - 1, cropped_data.shape[1] - 1)
    }

    # Mask for valid center cells (non-no_data)
    valid_center_cells = (cropped_data != no_data)

    # Extract x, y, z values for valid cells
    x_valid = cols[valid_center_cells]
    y_valid = rows[valid_center_cells]
    z_valid = cropped_data[valid_center_cells]

    # Create interpolator from valid points
    interpolator = NearestNDInterpolator(list(zip(x_valid, y_valid)), z_valid)

    # Check each corner for no data and interpolate if necessary
    for corner_name, (row, col) in corners.items():
        if cropped_data[row, col] == no_data:
            # Interpolate the nearest valid value
            cropped_data[row, col] = interpolator((col, row))

    # Extract non-no_data and center cells again after filling corners
    valid_center_cells = (cropped_data != no_data)

    # Extract final x, y, z values after filling corners
    x_filled = cols[valid_center_cells]
    y_filled = rows[valid_center_cells]
    z_filled = cropped_data[valid_center_cells]

    # Prepare final list of [x, y, z]
    xyz_filled = []
    for x_i, y_i, z_i in zip(x_filled, y_filled, z_filled):
        xyz_filled.append([x_i, y_i, z_i])

    return xyz_filled


def fill_raster(cropped_data, nodata_value, transform):
    """
    Fill the no data values of a given raster using Laplace interpolation.
    ----
    Input:
    - cropped_data (2d numpy array): cropped raster data.
    - nodata_value (int): nodata value to replace NAN after interplation with.
    - transform (rasterio transform): affine transform matrix.

    Output:
    - new_data[0, 1:-1, 1:-1] (2d numpy array): filled raster data with first and last rows and columns remove to ensure
                                                there are no nodata values.
    - new_transform (rasterio transform): affine transform matrix reflecting the one column one row removal shift.
    """

    # creating delaunay
    points = extract_center_cells(cropped_data, no_data=nodata_value)
    dt = startinpy.DT()
    dt.insert(points, "BBox")

    # now interpolation
    new_data = np.copy(cropped_data)

    # for interpolation, grid of all column and row positions, excluding the first and last rows/cols
    cols, rows = np.meshgrid(
        np.arange(1, cropped_data.shape[2] - 1),
        np.arange(1, cropped_data.shape[1] - 1)

    )

    # flatten the grid to get a list of all (col, row) locations
    locs = np.column_stack((cols.ravel(), rows.ravel()))

    # laplace interpolation
    interpolated_values = dt.interpolate({"method": "Laplace"}, locs)

    # reshape interpolated grid back to original
    interpolated_grid = np.reshape(interpolated_values, (cropped_data.shape[1] - 2, cropped_data.shape[2] - 2))

    # fill new_data with interpolated values
    new_data[0, 1:-1, 1:-1] = interpolated_grid
    new_data = np.where(np.isnan(new_data), nodata_value, new_data)

    new_transform = transform * Affine.translation(1, 1)

    return new_data[0, 1:-1, 1:-1], new_transform


def chm_finish(chm_array, dtm_array, transform, min_height=2, max_height=40):
    """
    Finish the CHM file by first removing the ground height. Then remove vegetation height
    below and above a certain range to ensure effective shade and remove noise.
    ----
    Input:
    - chm_array (2d numpy array):       cropped raster array of the CHM.
    - dtm_array (2d numpy array):       cropped raster array of the filled DSM.
    - transform (rasterio transform):   affine transform matrix.
    - min_height (float, optional):     minimal height for vegetation to be included.
    - max_height (float, optional):     maximum height for vegetation to be included.

    Output:
    - result_array (2d numpy array):    Array of the CHM with normalized height and min and max heights removed.
    - new_transform (rasterio transform): affine transform matrix reflecting the one column one row removal shift.
    """

    result_array = chm_array[0, 1:-1, 1:-1] - dtm_array
    result_array[(result_array < min_height) | (result_array > max_height)] = 0
    result_array[np.isnan(result_array)] = 0

    new_transform = transform * Affine.translation(1, 1)

    return result_array, new_transform


def replace_buildings(filled_dtm, dsm_buildings, buildings_geometries, transform):
    """
    Replace the values of the filled dtm with the values of the filled dsm, if there is a building.
    ----
    Input:
    - filled_dtm (2d np array):         filled array of the cropped AHN dtm.
    - dsm_buildings (2d np array):      Filled array of the cropped AHN dsm.
    - building_geometries (list):       A list of the building geometries
    - transform (rasterio transform):   affine transform matrix.

    Output:
    - final_dsm (2d numpy array):   a np array representing the final dsm, containing only ground and building
                                    heights.

    """

    building_mask = geometry_mask(buildings_geometries, transform=transform, invert=False, out_shape=filled_dtm.shape)

    # Apply the mask to the filled DTM
    final_dsm = np.where((building_mask), filled_dtm, dsm_buildings)

    return final_dsm


def load_buildings(buildings_path, layer):
    """
    Load in the building shapes from a geopackage file.
    ----
    Input:
    - buildings_path (string):   path to the geopackage file.
    - layer (string):            (Tile) name of the layer of buildings to be used

    Output:
    - List of dictionaries: A list of geometries in GeoJSON-like dictionary format.
      Each dictionary represents a building geometry with its spatial coordinates.
    """
    buildings_gdf = gpd.read_file(buildings_path, layer=layer)
    return [mapping(geom) for geom in buildings_gdf.geometry]


def extract_tilename(filename):
    """
    Extract the name of the AHN tile from the file name
    ----
    Input:
    - filename(string): the name of the input chm.tif.

    output:
    - match.group(1):   the name of the AHN tile.
    """

    match = re.match(r'CHM_(\w+)_\d+\.TIF', filename)
    if match:
        return match.group(1)
    return None


def process_all_folders(input_base_folder, dtm_path, dsm_path, buildings_path, output_base_folder, nodata_value=-9999):
    """
    Iterate through each subfolder in the input base folder and process TIF files.
    ------
    Input:
    - input_base_folder (string): Path to the folder containing subfolders of CHM .tif files.
    - dtm_path (string):          Path to the DTM .tif file.
    - dsm_path (string):          Path to the DSM .tif file.
    - buildings_path (string):    Path to the geopackage file containing building geometries.
    - output_base_folder (string): Base path for saving the output DSM and CHM files.
    - nodata_value (int):         NoData value for raster processing (default: -9999).
    """
    # List all subfolders in the input base folder
    subfolders = [os.path.join(input_base_folder, subfolder) for subfolder in os.listdir(input_base_folder)
                  if os.path.isdir(os.path.join(input_base_folder, subfolder))]

    if not subfolders:
        print(f"No subfolders found in {input_base_folder}.")
        return

    with tqdm.tqdm(total=len(subfolders), desc="Processing Subfolders", unit="subfolder") as pbar:
        # Process each subfolder
        for chm_folder in subfolders:
            print(f"\nProcessing CHM files in folder: {chm_folder}")
            process_files(
                chm_folder=chm_folder,
                dtm_path=dtm_path,
                dsm_path=dsm_path,
                buildings_path=buildings_path,
                output_base_folder=output_base_folder,
                nodata_value=nodata_value
            )
            pbar.update(1)


def process_files(chm_folder, dtm_path, dsm_path, buildings_path, output_base_folder, nodata_value=-9999):
    """
      Function to run the whole process of creating the final DSM and CHM.
      ----
      Input:
      - chm_folder (string):         Path to the folder containing input CHM .tif files.
      - dtm_path (string):           Path to the DTM .tif file.
      - dsm_path (string):           Path to the DSM .tif file.
      - buildings_path (string):     Path to the geopackage file containing building geometries.
      - output_base_folder (string): Base path for saving the output DSM and CHM files.
      - nodata_value (int):          NoData value for raster processing (default: -9999).

      Output:
      - None: The function writes output files directly to the specified `output_base_folder`.
      - Output files:
          - For each input CHM file, the function generates:
              - A DSM file with building data incorporated: `DSM_<tile>_<subtile_number>.tif`
              - A processed CHM file: `CHM_<tile>_<subtile_number>.tif`
          - These files are saved in folders named after the tile in `output_base_folder`.
      """

    chm_files = glob.glob(os.path.join(chm_folder, "*.TIF"))
    if not chm_files:
        print(f"No CHM files found in the folder: {chm_folder}")
        return

    first_chm_filename = os.path.basename(chm_files[0])
    tile = extract_tilename(first_chm_filename)

    # load in the building geometries .
    building_geometries = load_buildings(buildings_path, layer=tile)

    # output folder for tile
    output_folder = os.path.join(output_base_folder, tile)
    os.makedirs(output_folder, exist_ok=True)
    total_start_time = time.time()

    # Use tqdm to wrap the file list for progress tracking
    for chm_path in tqdm.tqdm(chm_files, desc="Processing Files", unit="file"):
        chm_filename = os.path.basename(chm_path)
        file_number = re.search(r'_(\d+)\.TIF', chm_filename).group(1)

        # overlap BBOX dtm, dsm & chm
        raster_paths = [dtm_path, chm_path, dsm_path]
        overlapping_bbox = get_bbox(raster_paths)
        print(f"\nProcessing tile {chm_filename}")
        start_time = time.time()

        # Cropping rasters to bbox
        dtm_cropped, dtm_transform, dtm_crs = crop_raster(dtm_path, overlapping_bbox, no_data=nodata_value, tile=tile,
                                                          file_number=file_number)
        dsm_cropped, dsm_transform, dsm_crs = crop_raster(dsm_path, overlapping_bbox, no_data=nodata_value, tile=tile,
                                                          file_number=file_number)
        chm_cropped, chm_transform, _ = crop_raster(chm_path, overlapping_bbox, no_data=nodata_value, tile=tile,
                                                          file_number=file_number)

        print(f"cropped at {(time.time() - start_time):.2f}")

        # fill no data values dtm & dsm
        filled_dtm, _ = fill_raster(dtm_cropped, nodata_value, dtm_transform)
        print(f"filled dtm at {(time.time() - start_time):.2f}")

        filled_dsm, new_dsm_transform = fill_raster(dsm_cropped, nodata_value, dsm_transform)
        print(f"filled dsm at {(time.time() - start_time):.2f}")

        # normalized CHM calculation
        chm_result, new_chm_transform = chm_finish(chm_cropped, filled_dtm, chm_transform)

        # Insert buildings from DSM in DTM
        final_dsm_with_buildings = replace_buildings(filled_dtm, filled_dsm, building_geometries, new_dsm_transform,
                                                     nodata_value)
        print(f"Inserted buildings at {(time.time() - start_time):.2f}")

        # Define output file paths
        output_dtm_filename = f"DSM_{tile}_{file_number}.tif"
        output_chm_filename = f"CHM_{tile}_{file_number}.tif"
        output_dtm_path = os.path.join(output_folder, output_dtm_filename)
        output_chm_path = os.path.join(output_folder, output_chm_filename)

        # Write output rasters
        write_output(rasterio.open(dtm_path), final_dsm_with_buildings, new_dsm_transform, output_dtm_path, True)
        write_output(rasterio.open(chm_path), chm_result, new_chm_transform, output_chm_path, True)

        elapsed_time = time.time() - start_time
        print(f"Processed {chm_filename} in {elapsed_time:.2f} seconds and saved output to {output_dtm_filename}")

    total_elapsed_time = time.time() - total_start_time
    print(f"\nAll files processed in {total_elapsed_time:.2f} seconds.")


def process_folders(base_chm_folder, dtm_path, dsm_path, buildings_path, output_base_folder, nodata_value=-9999):
    """
    Process each folder containing CHM files concurrently.
    -----------------
    Input:
    - base_chm_folder (string):     Path to the base folder containing subfolders of CHM files.
    - dtm_path (string):            Path to the DTM .tif file.
    - dsm_path (string):            Path to the DSM .tif file.
    - buildings_path (string):      Path to the geopackage file containing building geometries.
    - output_base_folder (string):  Base path for saving the output DSM and CHM files.
    - nodata_value (int):           NoData value for raster processing (default: -9999).

    Output:
    - None: The function writes output files directly to the specified `output_base_folder`.
    """

    # Iterate over each folder in the base folder
    for root, dirs, files in os.walk(base_chm_folder):
        for dir_name in dirs:
            chm_folder = os.path.join(root, dir_name)
            print(f"Processing folder: {chm_folder}")
            process_files(chm_folder, dtm_path, dsm_path, buildings_path, output_base_folder, nodata_value)

"""
geotiff_dtm = "data/DTM_ams.tif"
geotiff_dsm = "data/DSM_ams.tif"
buildings = "data/ams_buildings.gpkg"
chm_folder = "output/25DN2_TEST"
output_base_folder = "final"

process_files(chm_folder, geotiff_dtm, geotiff_dsm, buildings, output_base_folder)
"""