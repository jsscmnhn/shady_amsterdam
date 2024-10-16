""" This file contains the functions to create the 'first CHM', which is the direct output of rasterizing the
point cloud data. This CHM will still contain no data values, and the ground height is not substracted from the
vegetation heights. This file can run the CHM creation process in parallel."""

import concurrent.futures
from shade_calculation.extra import functions
import os
import laspy
import time
from tqdm import tqdm
import numpy as np
import startinpy
from scipy.spatial import cKDTree
from scipy.ndimage import median_filter


def median_filter_chm(chm_array, nodata_value=-9999, size=3):
    """
    Apply a median filter to a CHM, handling NoData values.
    -----
    Parameters:
    - chm_array (np.ndarray): The array representing the height values of the CHM.
    - nodata_value (float): Value representing NoData in the input raster.
    - size (int): Size of the median filter. It defines the footprint of the filter.

    Returns:
    - smoothed_chm (np.ndarray): The smoothed CHM array.
    """
    # Create a mask for valid data
    valid_mask = chm_array != nodata_value

    # Pad the data with nodata_value
    pad_width = size // 2
    padded_chm = np.pad(chm_array, pad_width, mode='constant', constant_values=nodata_value)

    # Apply median filter to padded data
    filtered_padded = median_filter(padded_chm.astype(np.float32), size=size)

    # Remove padding
    smoothed_chm = filtered_padded[pad_width:-pad_width, pad_width:-pad_width]

    # Only keep valid data in smoothed result
    smoothed_chm[~valid_mask] = nodata_value

    return smoothed_chm

def extract_vegetation_points(LasData, ndvi_threshold=0.1, pre_filter=False):
    """
    Extract vegetation points based on classification and NDVI threshold.
    ------
    Input:
    - LasData (laspy.LasData): Input point cloud data in LAS format.
    - ndvi_threshold (float): The NDVI threshold for identifying vegetation points.
                              NDVI values greater than this threshold are considered vegetation.
    - pre_filter (bool): If True, applies an additional filter to remove vegetation points below a certain height
                         threshold (1.5 meters above the lowest vegetation point).
    Output:
    - laspy.LasData: A new LasData object containing only the filtered vegetation points based on the specified criteria.
    """

    # Filter points based on classification (vegetation-related classes), note: vegetation classes are empty in AHN4
    possible_vegetation_points = LasData[(LasData.classification == 1) |  # Unclassified
                                         (LasData.classification == 3) |  # Low vegetation
                                         (LasData.classification == 4) |  # Medium vegetation
                                         (LasData.classification == 5)]  # High vegetation

    # Calculate NDVI

    red = possible_vegetation_points.red
    nir = possible_vegetation_points.nir
    ndvi = (nir.astype(float) - red) / (nir + red)

    # Filter the points whose NDVI is greater than the threshold
    veg_points = possible_vegetation_points[ndvi > ndvi_threshold]

    # Option: already filter away the points with a height below 1.5m from the lowest veg point, introduced because
    # of one very large tile (25GN2_24.LAZ)
    if pre_filter:
        heights =veg_points.z
        min_height = heights.min()

        # Filter out points with heights between the minimum height and 1.5 meters
        filtered_veg_points = veg_points[(heights <= min_height) | (heights > 1.5)]
        return filtered_veg_points

    return veg_points


def chm_creation(LasData, vegetation_data, output_filename, resolution=0.5, smooth=False, nodata_value=-9999, filter_size=3):
    """
    Create a CHM from LiDAR vegetation data and save it as a raster.
    -------
    Input:
    - LasData (laspy.LasData):      Input LiDAR point cloud data used for metadata and output CRS.
    - vegetation_data (tuple):      A tuple containing:
                        - veg_raster (numpy.ndarray): The array representing the height values of vegetation.
                        - grid_centers (tuple of numpy.ndarrays): Contains two arrays (x, y) with the coordinates
                          of the center points of each grid cell.
    - output_filename (str): The name of the output .tif file for saving the CHM.
    - resolution (float, optional): The spatial resolution of the output raster in the same units as the input data
                                    (default: 0.5).
    - smooth (bool, optional): If True, applies a median filter to smooth the CHM.
    - nodata_value (float, optional): The value for NoData pixels (default: -9999).
    - filter_size (int, optional): Size of the median filter (default: 3).

    Output:
    - None: The function saves the CHM as a raster file (.tif) to the specified output path.
    """
    veg_raster = vegetation_data[0]
    grid_centers = vegetation_data[1]
    top_left_x = grid_centers[0][0, 0] - resolution / 2
    top_left_y = grid_centers[1][0, 0] + resolution / 2

    transform = functions.create_affine_transform(top_left_x, top_left_y, resolution)

    if smooth:
        veg_raster = median_filter_chm(veg_raster, nodata_value=nodata_value, size=filter_size)

    functions.write_output(LasData, veg_raster, transform, output_filename, True)


def interpolation_vegetation(LasData, veg_points, resolution, no_data_value=-9999):
    """
    Create a vegetation raster using Laplace interpolation.

    InpurL
    - LasData (laspy.LasData):          Input LiDAR point cloud data.
    - veg_points (laspy.LasData):       Vegetation points to be interpolated.
    - resolution (float):               Resolution of the raster.
    - no_data_value (int, optional):    Value for no data

    Returns:
    - interpolated_grid (np.ndarray): Generated raster for vegetation.
    - grid_center_xy (tuple): Grid of x, y center coordinates for each raster cell.
    """

    # Extents of the pc
    min_x, max_x = round(LasData.x.min()), round(LasData.x.max())
    min_y, max_y = round(LasData.y.min()), round(LasData.y.max())

    # Define size of the region
    x_length = max_x - min_x
    y_length = max_y - min_y

    # Number of rows and columns
    cols = round(x_length / resolution)
    rows = round(y_length / resolution)

    # Initialize raster grid
    vege_raster = np.full((rows, cols), no_data_value, dtype=np.float32)

    # Calculate center coords for each grid cell
    grid_center_xy = functions.raster_center_coords(min_x, max_x, min_y, max_y, resolution)

    if veg_points.x.shape[0] == 0:
        print("There are no vegetation points in the current area.")
        vege_raster = np.full((rows, cols), 0, dtype=np.float32)
        return vege_raster, grid_center_xy

    # create the delaunay triangulation
    dt = startinpy.DT()
    dt.insert(veg_points.xyz, "BBox")

    # Flatten the grid to get a list of all center coords
    locs = np.column_stack((grid_center_xy[0].ravel(), grid_center_xy[1].ravel()))

    vegetation_points = np.column_stack((veg_points.x, veg_points.y))
    tree = cKDTree(vegetation_points)

    # Find the distance to the nearest vegetation point for each grid cell
    distances, _ = tree.query(locs, k=1)

    distance_threshold = 1
    # masking cells that exceed threshold
    within_threshold_mask = distances <= distance_threshold
    # Interpolation only for those near
    valid_locs = locs[within_threshold_mask]

    # laplace interpolation
    interpolated_values = dt.interpolate({"method": "Laplace"}, valid_locs)

    # reshape interpolated grid back to og
    interpolated_grid = np.full_like(vege_raster, no_data_value, dtype=np.float32)  # Start with no_data
    interpolated_grid.ravel()[within_threshold_mask] = interpolated_values

    return interpolated_grid, grid_center_xy


def process_single_laz_file(file_path, output_folder, ndvi_threshold=0.0, resolution=0.5, remove=False,
                            smooth_chm=False, filter_size=3, pre_filter=False):
    """
    Process a LAZ file to extract vegetation points and generate a CHM.
    -------
    Input:
    - file_path (str):          The file_path containing the folders containing the input .LAZ file.
    - output_folder (str):      The folder where the output CHM .tif files will be saved.
    - ndvi_threshold (float):   The NDVI threshold for classifying vegetation points.
    - resolution (float):       The resolution of the output CHM rasters, defining the size of each pixel (default: 0.5).
    - remove (bool):            If True, deletes the original .LAZ files after processing (default: False).
    - smooth_chm (bool):        If True, applies smoothing to the CHM using a median filter (default: False).
    - filter_size (int):        Size of the median filter to use if smoothing (default: 3).

    Output:
    - None: The function process a .LAZ file, creates a corresponding CHM .tif file, and saves it to the output folder.
            Optionally deletes the original .LAZ files if `remove` is set to True.
    """
    tile_name = os.path.splitext(os.path.basename(file_path))[0]
    tile = tile_name.split("_")[0]

    # Create a subfolder for this tile's outputs
    tile_output_folder = os.path.join(output_folder, tile)
    os.makedirs(tile_output_folder, exist_ok=True)

    # Define the output filename in the tile-specific subfolder
    output_filename = os.path.join(tile_output_folder, f"CHM_{tile_name}.TIF")

    print(f"Processing tile {tile_name}")
    start_time = time.time()

    # Load LAS data
    with laspy.open(file_path) as las:
        LasData = las.read()
        print(f"Loaded {file_path} with {len(LasData.points)} points.")

    # check if large tile
    if tile_name == '25GZ2_09':
        pre_filter = True
    # Extract vegetation points
    veg_points = extract_vegetation_points(LasData, ndvi_threshold=ndvi_threshold, pre_filter=pre_filter)

    # Perform interpolation for vegetation data
    vegetation_data = interpolation_vegetation(LasData, veg_points, resolution)

    # Create the CHM and save it
    chm_creation(
        LasData,
        vegetation_data,
        output_filename=output_filename,
        resolution=resolution,
        smooth=smooth_chm,
        nodata_value=-9999,
        filter_size=filter_size
    )

    elapsed_time = time.time() - start_time
    print(
        f"Processed {os.path.basename(file_path)} in {elapsed_time:.2f} seconds and saved output to {output_filename}")

    # Optionally remove the original file
    if remove:
        os.remove(file_path)
        print(f"Deleted file: {file_path}")

    return file_path

def process_laz_files(input_folder, output_folder, ndvi_threshold=0.0, resolution=0.5, remove=False, smooth_chm=False,
                      filter_size=3, pre_filter=False, max_workers=4):
    """
    Process a folder of LAZ files in parallel to extract vegetation points and generate CHMs.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    laz_files = []
    for root, dirs, files in os.walk(input_folder):
        for file in files:
            if file.endswith(".LAZ"):
                laz_files.append(os.path.join(root, file))

    if not laz_files:
        print("No LAZ files found in the input folder or its subfolders.")
        return


    total_start_time = time.time()

    # Use ProcessPoolExecutor for parallel processing
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        list(tqdm(
            executor.map(
                process_single_laz_file,
                laz_files,
                [output_folder] * len(laz_files),
                [ndvi_threshold] * len(laz_files),
                [resolution] * len(laz_files),
                [remove] * len(laz_files),
                [smooth_chm] * len(laz_files),
                [filter_size] * len(laz_files),
                [pre_filter] * len(laz_files)
            ),
            total=len(laz_files),
            desc="Processing files",
            unit="file"
        ))
    total_elapsed_time = time.time() - total_start_time
    print(f"\nAll files processed in {total_elapsed_time:.2f} seconds.")



if __name__ == '__main__':
    input_folder = "E:/temporary_jessica/missed_laz_tiles"
    output_folder = "E:/temporary_jessica/CHM_smoothed"
    max_workers = 20
    process_laz_files(input_folder, output_folder, max_workers=max_workers)

