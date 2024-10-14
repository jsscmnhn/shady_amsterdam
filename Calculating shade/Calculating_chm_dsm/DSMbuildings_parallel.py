import os
import glob
import re
import geopandas as gpd
import numpy as np
import rasterio
from shapely.geometry import mapping
from rasterio.features import geometry_mask
from scipy.interpolate import NearestNDInterpolator
from shade_calculation.src.functions import write_output
import tqdm
import startinpy
import time
from rasterio import Affine
from concurrent.futures import ProcessPoolExecutor


def get_bbox(raster_paths):
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


def crop_raster(raster_path, bbox, no_data=-9999, file_number=None, tile=None):
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

        # Check for specific tile conditions for cropping
        if tile == "25DN2" and file_number in ["25", "24"]:
            cropped_data = cropped_data[:, :2150, :]

        return cropped_data, src.window_transform(window), src.crs


def extract_center_cells(cropped_data, no_data=-9999):
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
    # creating delaunay
    points = extract_center_cells(cropped_data, no_data=nodata_value)
    dt = startinpy.DT()
    dt.insert(points, "BBox")
    print("dt")

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

    # reshape interpolated grid back to og
    interpolated_grid = np.reshape(interpolated_values, (cropped_data.shape[1] - 2, cropped_data.shape[2] - 2))

    # fill new_data with interpolated values
    new_data[0, 1:-1, 1:-1] = interpolated_grid
    new_data = np.where(np.isnan(new_data), nodata_value, new_data)

    new_transform = transform * Affine.translation(1, 1)

    return new_data[0, 1:-1, 1:-1], new_transform


def chm_finish(chm_array, dtm_array, transform):
    result_array = chm_array[0, 1:-1, 1:-1] - dtm_array
    result_array[(result_array < 2) | (result_array > 40)] = 0

    new_transform = transform * Affine.translation(1, 1)

    return result_array, new_transform


def replace_buildings(filled_dtm, dsm_buildings, buildings_geometries, transform, nodata_value=-9999):
    building_mask = geometry_mask(buildings_geometries, transform=transform, invert=False, out_shape=filled_dtm.shape)

    # Apply the mask to the filled DTM
    final_dtm = np.where((building_mask), filled_dtm, dsm_buildings)

    return final_dtm


def load_buildings(buildings_path, layer):
    buildings_gdf = gpd.read_file(buildings_path, layer=layer)
    return [mapping(geom) for geom in buildings_gdf.geometry]


def extract_tilename(filename):
    match = re.match(r'CHM_(\w+)_\d+\.TIF', filename)
    if match:
        return match.group(1)
    return None


def process_single_file(chm_path, dtm_path, dsm_path, building_geometries, output_base_folder, nodata_value=-9999):
    chm_filename = os.path.basename(chm_path)
    file_number = re.search(r'_(\d+)\.TIF', chm_filename).group(1)
    tile = extract_tilename(chm_filename)

    if not tile:
        print(f"Skipping {chm_filename}, couldn't extract common part.")
        return

    # Output folder for tile
    output_folder = os.path.join(output_base_folder, tile)
    os.makedirs(output_folder, exist_ok=True)

    # Overlap BBOX dtm, dsm & chm
    raster_paths = [dtm_path, chm_path, dsm_path]
    overlapping_bbox = get_bbox(raster_paths)

    # Cropping rasters to bbox
    dtm_cropped, dtm_transform, dtm_crs = crop_raster(dtm_path, overlapping_bbox, no_data=nodata_value, tile=tile,
                                                      file_number=file_number)
    dsm_cropped, dsm_transform, dsm_crs = crop_raster(dsm_path, overlapping_bbox, no_data=nodata_value, tile=tile,
                                                      file_number=file_number)
    chm_cropped, chm_transform, _ = crop_raster(chm_path, overlapping_bbox, no_data=nodata_value, tile=tile,
                                                file_number=file_number)

    # Fill raster with DTM and DSM data
    filled_dtm, _ = fill_raster(dtm_cropped, nodata_value, dtm_transform)
    filled_dsm, new_dsm_transform = fill_raster(dsm_cropped, nodata_value, dsm_transform)

    # Norm CHM calculation
    chm_result, new_chm_transform = chm_finish(chm_cropped, filled_dtm, chm_transform)

    # Insert buildings from DSM in DTM
    final_dsm_with_buildings = replace_buildings(filled_dtm, filled_dsm, building_geometries, new_dsm_transform,
                                                 nodata_value)

    # Output file names
    output_dtm_filename = f"DSM_{tile}_{file_number}.tif"
    output_chm_filename = f"CHM_{tile}_{file_number}.tif"

    # Saving final DTM + buildings
    output_dtm_path = os.path.join(output_folder, output_dtm_filename)
    write_output(rasterio.open(dtm_path), final_dsm_with_buildings, new_dsm_transform, output_dtm_path, True)

    # Saving CHM result
    output_chm_path = os.path.join(output_folder, output_chm_filename)
    write_output(rasterio.open(chm_path), chm_result, new_chm_transform, output_chm_path, True)

    print(f"Processed {chm_filename} and saved output to {output_dtm_filename}")


def process_files(chm_files, dtm_path, dsm_path, buildings_path, output_base_folder, nodata_value=-9999):
    total_start_time = time.time()

    # Load building geometries once
    first_chm_filename = os.path.basename(chm_files[0])
    tile = extract_tilename(first_chm_filename)

    if not tile:
        print(f"Skipping {first_chm_filename}, couldn't extract tile information.")
        return

    # Load building geometries once for the tile
    building_geometries = load_buildings(buildings_path, layer=tile)

    # Set up a ProcessPoolExecutor to parallelize the file processing
    with ProcessPoolExecutor() as executor:
        # Use tqdm for progress tracking and parallel submission of tasks
        list(tqdm.tqdm(executor.map(
            process_single_file,
            chm_files,  # CHM files are iterated over
            [dtm_path] * len(chm_files),  # All share the same DTM path
            [dsm_path] * len(chm_files),  # All share the same DSM path
            [building_geometries] * len(chm_files),  # Pass the building geometries loaded once
            [output_base_folder] * len(chm_files),  # All share the same output base folder
            [nodata_value] * len(chm_files)  # All use the same nodata value
        ), total=len(chm_files), desc="Processing Files", unit="file"))

    total_elapsed_time = time.time() - total_start_time
    print(f"\nAll files processed in {total_elapsed_time:.2f} seconds.")

geotiff_dtm = "data/DTM_ams.tif"
geotiff_dsm = "data/DSM_ams.tif"
buildings = "data/ams_buildings.gpkg"

chm_folder = "D:/Geomatics/CHM/25DN2"
output_base_folder = "D:/Geomatics/final"

chm_files = glob.glob(os.path.join(chm_folder, "*.TIF"))
process_files(chm_files, geotiff_dtm, geotiff_dsm, buildings, output_base_folder, nodata_value=-9999)
